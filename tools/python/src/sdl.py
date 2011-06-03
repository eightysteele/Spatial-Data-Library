#!/usr/bin/env python

# Copyright 2011 Jante LLC and University of Kansas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__author__ = "Aaron Steele"

"""
This module includes classes and a command line interface for bulkloading 
WorldClim environment variables to CouchDB.
"""

import csv
import couchdb
import logging
import math
from optparse import OptionParser
import os
import random
import shapefile
import shlex
import subprocess

WORLDCLIM_TILE_RESOLUTION_DEGREES = 30.0

class Variable(object):
    """An environmental variable backed by a .bil and a .hdr file."""

    def __init__(self, bilfile, hdrfile):
        """Constructs a Variable.

        Arguments:
            bilfile - The .bil file path.
            hdrfile - The .hdr file path.
        """
        self.bilfile = bilfile
        self.hdrfile = hdrfile
        
        # Loads xmin, xmax, ymin, and ymax values from the .hdr file:
        for line in open(hdrfile, 'r'):
            if line.startswith('MaxX'):
                self.xmax = int(line.split()[1].strip())
            elif line.startswith('MinX'):
                self.xmin = int(line.split()[1].strip())
            elif line.startswith('MaxY'):
                self.ymax = int(line.split()[1].strip())
            elif line.startswith('MinY'):
                self.ymin = int(line.split()[1].strip())


class TileCell(object):
    """A tile cell described by a polygon with geographic coordinates."""

    def __init__(self, resolution, key, polygon):
        """Constructs a TileCell.

        Arguments:
        """
        self.resolution = resolution
        self.key = key
        self.polygon = polygon
        
    def __str__(self):
        return str(self.__dict__)
    
class Tile(object):
    """A 30 arc-second geographic tile (http://www.worldclim.org/tiles.php)."""

    def __init__(self, row, col, filename=None):
        """Tile constructor.

        Arguments:
            row - The global tile row number.
            col - The global tile column number.
            filename - Shapefile file name.
        """
        self.resolution = WORLDCLIM_TILE_RESOLUTION_DEGREES
        self.row = row
        self.col = col
        self.filename = filename
        self.north = 90.0 - (self.resolution * row)
        self.south = self.north - self.resolution
        self.west = -180.0 + (self.resolution * col)
        self.east = self.west + 30.0        

    def _getcellkey(self, tilerow, tilecol, cellres):
        """Gets the global cell key."""
        multiplier = self.resolution / cellres
        x = int(math.floor((self.col * multiplier) + tilecol))
        y = int(math.floor((self.row * multiplier) + tilerow))
        return '%s-%s' % (y, x)

    def _getcellpolygon(self, lat, lng, cellres):
        """Gets the cell polygon for a lat, lng, and cell resolution."""
        return [
                 [lng, lat],                       
                 [lng + cellres, lat],
                 [lng + cellres, lat - cellres],
                 [lng, lat - cellres],
                 [lng, lat]
               ]

    def __str__(self):
        return str(self.__dict__)

    def _clip2intersect2couchdb(self, cells, options, batchnum):
        logging.info(self.filename)
        filename = os.path.join(os.path.splitext(self.filename)[0], '%s' % batchnum)
        logging.info(filename)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        for cell in cells:
            w.poly(parts=[cell.polygon])
            w.record(CellKey=cell.key)
        w.save(filename)        
        clippedfile = Tile.clip2cell('%s.shp' % filename, self.filename)
        csvfile = Tile.intersect(clippedfile, options)
        server = couchdb.Server(options.couchurl)
        cdb = server['sdl-dev']    
        Tile.csv2couch(csvfile, cdb)

    @classmethod
    def csv2couch(cls, csvfile, cdb):
        logging.info('Bulkloading csv file ' + csvfile)
        dr = csv.DictReader(open(csvfile, 'r'))
        cells = {}
        dr.next() # Skip header
        for row in dr:
            cellkey = row.get('CellKey')
            if not cells.has_key(cellkey):
                cells[cellkey] = {
                    '_id': cellkey, 
                    'coords': [[random.uniform(x, 90),random.uniform(-180, x)] for x in range(-2,3)], # TODO
                    'vars': {}
                    }            
            varname = row.get('RID').split('_')[0]
            varval = row.get('Band1')
            cells.get(cellkey).get('vars')[varname] = varval
        logging.info('Bulkloading %s documents' % len(cells))
        cdb.update(cells.values())

    @classmethod
    def intersect(cls, shapefile, options):      
        """Intesects features in a shapefile with variables via starspan."""
        variables = [os.path.join(options.vardir, x) \
                         for x in os.listdir(options.vardir) \
                         if x.endswith('.bil')]
        variables = reduce(lambda x,y: '%s %s' % (x, y), variables)
        csvfile = shapefile.replace('.shp', '.csv')
        command = 'starspan --vector %s --raster %s --csv %s' \
            % (shapefile, variables, csvfile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        return csvfile
        
    @classmethod
    def clip2cell(cls, src, shapefile):
        """Clips src by shapefile and returns clipped shapefile name."""
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        clipped = src.replace('.shp', '-clipped.shp')
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, shapefile, clipped, src)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        return clipped
                
    def bulkload2couchdb(self, options):
        """Bulkloads the tile to CouchDB using command line options."""
        batchsize = int(options.batchsize)
        batchnum = 0
        cellres = float(options.cellres)
        cells = []
        count = 0
        for cell in self.getcells(cellres):
            if count >= batchsize:
                self._clip2intersect2couchdb(cells, options, batchnum)
                count = 0
                cells = []
                batchnum += 1
                continue
            cells.append(cell)
            count += 1
        if count > 0:
            self._clip2intersect2couchdb(cells, options, batchnum)
    
    def clip(self, shapefile, workspace):
        """Clips shapefile against tile and returns clipped Tile object."""
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        this = self.writeshapefile(workspace)
        clipped = this.replace('.shp', '-clipped.shp')
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, this, clipped, shapefile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        return Tile(self.row, self.col, clipped)

    def writeshapefile(self, workspace):
        """Writes tile shapefile in workspace directory and returns filename."""
        cell = self.getcells(30.0).next()            
        fout = os.path.join(workspace, cell.key)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        w.poly(parts=[cell.polygon])
        w.record(CellKey=cell.key)
        w.save(fout)        
        return '%s.shp' % fout

    def writemultishapefiles(self):
        pass
                                            
    def getcells(self, cellres):
        """Iterates over all cells in the tile at the given cellres.

        Arguments:
            cellres - The cell resolution.
        """
        lng = self.west
        lat = self.north        
        row = 0
        col = 0
        while lng < self.east:
            row = 0
            while lat > self.south:
                yield TileCell(
                    cellres,
                    self._getcellkey(row, col, cellres), 
                    self._getcellpolygon(lat, lng, cellres))
                lat -= cellres
                row += 1
            lat = self.north
            lng += cellres
            col += 1

def _getoptions():
    """Parses command line options and returns them."""
    parser = OptionParser()    
    parser.add_option("-v", 
                      "--vardir", 
                      dest="vardir",
                      help="The directory of variable files.",
                      default=None)
    parser.add_option("-w", 
                      "--workspace", 
                      dest="workspace",
                      help="The workspace directory for temporary files.",
                      default=None)
    parser.add_option("-c", 
                      "--couchurl", 
                      dest="couchurl",
                      help="The CouchDB URL.",
                      default=None)
    parser.add_option("-g", 
                      "--gadm", 
                      dest="gadm",
                      help="The GADM shapefile.",
                      default=None)
    parser.add_option("-t", 
                      "--tile", 
                      dest="tile",
                      help="The Worldclim tile number (rowcol).",
                      default=None)
    parser.add_option("-b", 
                      "--batchsize", 
                      dest="batchsize",
                      help="The batch size (default 25,000)",
                      default=None)
    parser.add_option("-r", 
                      "--cell-resolution", 
                      dest="cellres",
                      help="The cell resolution",
                      default=None)

    return parser.parse_args()[0]

def load(options):    
    row = int(options.tile.split(',')[0])
    col = int(options.tile.split(',')[1])
    tile = Tile(row, col)
    clipped = tile.clip(options.gadm, options.workspace)
    logging.info('Clipped: %s' % str(clipped))
    clipped.bulkload2couchdb(options)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    options = _getoptions()
    load(options)    
    
    
