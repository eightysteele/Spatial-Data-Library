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

__author__ = "Aaron Steele, Dave Vieglais, and John Wieczorek"

"""
This module includes classes and a command line interface for bulkloading 
WorldClim environment variables to CouchDB using the Rectangular Mesh Grid (RMG).
"""

import csv
import time
import couchdb
import logging
import math
from optparse import OptionParser
import os
import random
import shapefile
import shlex
import subprocess
from rmg import *

def maketile(options):
    key = options.key
    nw = map(float, options.nwcorner.split(','))
    se = map(float, options.secorner.split(','))
    nwcorner = Point(nw[0], nw[1])
    secorner = Point(se[0], se[1])
    cells_per_degree = float(options.cells_per_degree)
    tile = Tile(key, nwcorner, secorner, cells_per_degree)
    return tile

def clip(options):
    key = options.key
    nw = map(float, options.nwcorner.split(','))
    se = map(float, options.secorner.split(','))
    nwcorner = Point(nw[0], nw[1])
    secorner = Point(se[0], se[1])
    cells_per_degree = float(options.cells_per_degree)
    tile = Tile(key, nwcorner, secorner, cells_per_degree)
    clipped = tile.clip(options.gadm, options.workspace)
    return clipped

def load(options, clipped):
    clipped.bulkload2couchdb(options)

def getpolygon(key, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
    return RMGCell.polygon(key, cells_per_degree, digits, a, inverse_flattening)

def translatevariable(varname, varval):
    newval = float(varval)
    if newval > 55537: # Actual value is a negative number greater than the nodata value of -9999
        newval = newval - 65536
    newval = truncate(newval,0)
    return newval

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
    """A cell for a Tile described by a polygon with geographic coordinates."""

    def __init__(self, key, polygon, cells_per_degree):
        """Constructs a TileCell.

        Arguments:
        """
        self.key = key
        self.polygon = polygon
        self.cells_per_degree = cells_per_degree
        
    def __str__(self):
        return str(self.__dict__)
    
class Tile(object):
    """A geographic tile defined by a geographic coordinate bounding box."""

    def __init__(self, key, nwcorner, secorner, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING, filename=None):
        """Tile constructor.

        Arguments:
            key - identifier for the tile
            nwcorner - the starting Point in the northwest corner of the Tile.
            secorner - the ending Point in the southeast corner of the Tile.
            cells_per_degree - the desired resolution of the grid
            digits - the number of digits of precision to retain in the coordinates
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
            filename - The name of the input Shapefile for the Tile.
        """
        self.key = key
        self.nwcorner = nwcorner
        self.secorner = secorner
        self.cells_per_degree = cells_per_degree
        self.digits = digits
        self.a = a
        self.inverse_flattening = inverse_flattening
        self.filename = filename

    def __str__(self):
        return str(self.__dict__)

    def _clip2intersect2couchdb(self, cells, options, batchnum):
        t0 = time.time()
        filename = os.path.join(os.path.splitext(self.filename)[0], '%s' % batchnum)
        logging.info('Preparing shapefile %s in _clip2intersect2couchdb().' % (filename) )
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        for cell in cells:
            w.poly(parts=[cell.polygon])
            w.record(CellKey=cell.key)
        w.save(filename)        
        t0 = time.time()
        logging.info('Shapefile %s prepared in %s' % (batchnum, filename, t1-t0))
        clippedfile = Tile.clip2cell('%s.shp' % filename, self.filename)
        csvfile = Tile.intersect(clippedfile, options)
        Tile.csv2couch(csvfile, options)

    def polygon(self):
        """Returns a polygon (list of Points) for the Tile."""
        n = float(truncate(self.nwcorner.lat, self.digits))
        s = float(truncate(self.secorner.lat, self.digits))
        w = float(truncate(self.nwcorner.lng, self.digits))
        e = float(truncate(self.secorner.lng, self.digits))
        return [(w, n), (w, s), (e, s), (e, n), (w, n)]

    @classmethod
    def csv2couch(cls, csvfile, options):
        """Loads values from csv file to couchdb."""
        t0 = time.time()
        logging.info('Beginning csv2couch(), preparing cells for bulkloading from %s.' % (csvfile) )
        server = couchdb.Server(options.couchurl)
        cbd = server[options.couchdbname]
        cells_per_degree = float(options.cells_per_degree)
        dr = csv.DictReader(open(csvfile, 'r'))
        cells = {}
        for row in dr:
            cellkey = row.get('CellKey')
            if not cells.has_key(cellkey):
                cells[cellkey] = {
                    '_id': cellkey, 
                    'coords': getpolygon(cellkey, cells_per_degree),
                    'vars': {}
                    }            
            varname = row.get('RID').split('_')[0]
            # the following dependent on running starspan with --stats %s avg
            varval = row.get('avg_Band1')
            # the following need for Worldclim because of the -9999 NODATA value.
            cells.get(cellkey).get('vars')[varname] = translatevariable(varname, varval)
        t1 = time.time()
        logging.info('%s cells prepared for upload in %s' % (len(cells), t1-t0))
        cdb.update(cells.values())
        t2 = time.time()
        logging.info('%s documents uploaded in %s' % (len(cells), t2-t1))

    @classmethod
    def intersect(cls, shapefile, options):      
        """Intersects features in a shapefile with variables via starspan."""
        t0 = time.time()
        logging.info('Beginning starspan statistics on %s.' % (shapefile) )
        variables = [os.path.join(options.vardir, x) \
                         for x in os.listdir(options.vardir) \
                         if x.endswith('.bil')]
        variables = reduce(lambda x,y: '%s %s' % (x, y), variables)
        csvfile = shapefile.replace('.shp', '.csv')
        # Call starspan requesting mean of variable, exclusing nodata values (-9999 in the file is the same as 55537)
        command = 'starspan --vector %s --raster %s --stats %s avg --nodata 55537' \
            % (shapefile, variables, csvfile)
#        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('starspan statistics finished in %s.' % (t1-t0) )
        return csvfile
        
    @classmethod
    def clip2cell(cls, src, shapefile):
        """Clips src by shapefile and returns clipped shapefile name."""
        t0 = time.time()
        logging.info('Beginning clipping of %s by %s.' % (src, shapefile) )
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        clipped = src.replace('.shp', '-clipped.shp')
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, shapefile, clipped, src)
#        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('%s clipped by %s in %s' % (src, shapefile, t1-t0))
        return clipped
                
    def bulkload2couchdb(self, options):
        """Bulkloads the tile to CouchDB using command line options."""
        t0 = time.time()
        logging.info('Beginning bulkload2couchdb().')
        batchsize = int(options.batchsize)
        batchnum = 0
        cells_per_degree = float(options.cells_per_degree)
        cells = []
        count = 0
        for cell in self.getcells():
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
        t1 = time.time()
        logging.info('Total elapsed time to bulkload2couchdb(): %s' % (t1-t0))
    
    def clip(self, shapefile, workspace):
        """Returns a Tile clipped by shapefile.
        
        Arguments:
            shapefile - the shape file with which to clip the Tile
            workspace - the directory to store the clipped shape file
        """
        t0 = time.time()
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        this = self.writetileshapefile(workspace)
        clipped = this.replace('.shp', '-clipped.shp')
        logging.info('Beginning clipping of %s by %s.' % (shapefile, this))
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, this, clipped, shapefile)
#        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('%s clipped by %s in %s' % (shapefile, this, t1-t0))
        return Tile(self.key, self.nwcorner, self.secorner, self.cells_per_degree, self.digits, self.a, self.inverse_flattening, clipped)

    def writetileshapefile(self, workspace):
        """Writes a shapefile for the Tile the filename.
        
        Arguments:
            workspace - the directory to store the clipped shape file
            """
        fout = os.path.join(workspace, self.key)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('TileKey','C','255')
        w.poly(parts=[self.polygon()])
        w.record(TileKey=self.key)
        w.save(fout)        
        return '%s.shp' % fout

    def writeshapefile(self, workspace):
        """Writes tile shapefile in workspace directory and returns filename."""
        cell = self.getcells().next()
        cellinfo = 'Cell: %s' % (cell)
        fout = os.path.join(workspace, cell.key)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        w.poly(parts=[cell.polygon])
        w.record(CellKey=cell.key)
        w.save(fout)        
        return '%s.shp' % fout

    def getcells(self):
        """Iterates over a set of polygons for cells intersecting a bounding box.""" 
        north = self.nwcorner.lat
        west = self.nwcorner.lng
        south = self.secorner.lat
        east = self.secorner.lng
        crosses_180 = False
        if west > east:
            crosses_180 = True
            crossed = False
            west -= 360
        lng = west
        lat = north
        # key for the NW corner of the tile
        key = RMGCell.key(lng180(lng), lat, self.cells_per_degree, self.a, self.inverse_flattening)
        indexes = key.split('-')
        x_index = int(indexes[0])
        y_index = int(indexes[1])

        while lat > south:
            while lng < east:
                key = str(x_index)+'-'+str(y_index)
                polygon = tuple([(float(x[0]), float(x[1])) for x in RMGCell.polygon(key, self.cells_per_degree, self.digits, self.a, self.inverse_flattening)])
#                yieldingthis = 'Yield CellKey=%s: lat: %s lng: %s' % (key, lat, lng)
#                logging.info(yieldingthis)
                yield TileCell(key, polygon, self.cells_per_degree)
                if crosses_180 == True and crossed == False:
                    elng = RMGCell.east(x_index, y_index, self.cells_per_degree, self.a, self.inverse_flattening) - 360
                    wlng = RMGCell.west(x_index, y_index, self.cells_per_degree, self.a, self.inverse_flattening) - 360
                    if elng > east:
                        crossed = True
                        lng = east
                    elif wlng > -180:
                        crossed = True
                        x_index = 0
                        lng = -180
                    else:
                        x_index += 1
                        lng = RMGCell.west(x_index, y_index, self.cells_per_degree, self.a, self.inverse_flattening) - 360
                        
                else:
                    x_index += 1
                    lng = RMGCell.west(x_index, y_index, self.cells_per_degree, self.a, self.inverse_flattening)
            crossed = False
            lng = west
            y_index += 1
            lat = RMGCell.north(y_index, self.cells_per_degree)
            key = RMGCell.key(lng180(lng), lat, self.cells_per_degree, self.a, self.inverse_flattening)
            indexes = key.split('-')
            x_index = int(indexes[0])

def _getoptions():
    """Parses command line options and returns them."""
    parser = OptionParser()
    parser.add_option("-c", "--command", dest="command",
                      help="SDL command",
                      default=None)
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
    parser.add_option("-p", 
                      "--csvfile", 
                      dest="csvfile",
                      help="A clipped variables csv file to intersect and load.",
                      default=None)
    parser.add_option("-u", 
                      "--couchurl", 
                      dest="couchurl",
                      help="The CouchDB URL.",
                      default=None)
    parser.add_option("-d", 
                      "--couchdbname", 
                      dest="couchdbname",
                      help="The CouchDB database name.",
                      default='worldclim-rmg')
    parser.add_option("-g", 
                      "--gadm", 
                      dest="gadm",
                      help="The Global Administrative shapefile to clip to.",
                      default=None)
    parser.add_option("-k", "--key", dest="key",
                      help="Identifier for the Tile",
                      default=None)
    parser.add_option("-f", "--nwcorner", dest="nwcorner",
                      help="NW corner of bounding box",
                      default=None)
    parser.add_option("-t", "--secorner", dest="secorner",
                      help="SW corner of bounding box",
                      default=None)
    parser.add_option("-n", "--cells-per-degree", dest="cells_per_degree",
                      help="Number of cells per degree",
                      default=None)
    parser.add_option("-b", 
                      "--batchsize", 
                      dest="batchsize",
                      help="The batch size (default 25,000)",
                      default=25000)
    return parser.parse_args()[0]

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    options = _getoptions()
    command = options.command.lower()
    
    if command == 'clip':
        clipped = clip(options)
        logging.info('Finished command clip.')

    if command == 'csv2couchdb':
        tile = maketile(options)
        csvfile = os.path.join(options.workspace, options.csvfile)
        tile.csv2couch(csvfile, options)
        logging.info('Finished command csv2couch.')

    if command == 'load':
        clipped = clip(options)
        load(options, clipped)    
        logging.info('Finished command load.')
#Command line:
#            ./sdl.py -c load -w ~/SDL/workspace -g ~/Spatial-Data-Library/data/gadm/Terrestrial15kBuffer.shp  -k 37 -f 30,0 -t 60,-30 -b 25000 -n 120 -v /SDL/worldclim/37 -u http://eighty.berkeley.edu:5984 &

#Before clearing ~SDL/workspace of Tile 37 files
#Filesystem           1K-blocks      Used Available Use% Mounted on
#/dev/sda2            223856640 124558204  87927100  59% /

#After clearing ~SDL/workspace of Tile 37 files
#Filesystem           1K-blocks      Used Available Use% Mounted on
#/dev/sda2            223856640 106658316 105826988  51% /

