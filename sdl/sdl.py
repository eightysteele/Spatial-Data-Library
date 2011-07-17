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
    """ Creates a Tile object from the command line arguments.
    
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                nwcorner - the northwest corner of the Tile
                secorner - the southeast corner of the Tile
                cells_per_degree - the number of cells in a degree of longitude at
                    the equator - defines a resolution for the grid system.
    """
    key = options.key
    nw = map(float, options.nwcorner.split(','))
    se = map(float, options.secorner.split(','))
    nwcorner = Point(nw[0], nw[1])
    secorner = Point(se[0], se[1])
    cells_per_degree = float(options.cells_per_degree)
    tile = Tile(key, nwcorner, secorner, cells_per_degree)
    return tile

def getworldclimtile(options):
    """ Downloads the files for all of the variables for the Tile given by the Tile key
        in the command line, unzips them, and removes the zipfile.
        
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                vardir - the path to the directory in which to store the Worldclim files
    """
    varset = ['tmean','tmin','tmax','prec','alt','bio']
    for var in varset:
        if os.path.exists(options.vardir):
            command = 'rm -r %s' % (options.vardir)
            logging.info(command)
            args = shlex.split(command)
            subprocess.call(args)

        command = 'mkdir %s' % (options.vardir)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)

        varfile = '%s_%s.zip' % (var, options.key)
        varpath = os.path.join(options.vardir, varfile)
        command = 'wget -P %s http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/%s' % (options.vardir, varfile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)

        command = 'unzip %s -d %s' % (varpath, options.vardir)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        
        command = 'rm %s' % (varpath)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)

def clip(options):
    """ Creates a clipped Tile object and associated files from the command line arguments.
        Returns a clipped Tile object.
    
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                nwcorner - the northwest corner of the Tile
                secorner - the southeast corner of the Tile
                cells_per_degree - the number of cells in a degree of longitude at
                    the equator - defines a resolution for the grid system.
                gadm - the shapefile to clip to (polygon of are containing Worldclim data)
    """
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
    """ Loads a clipped Tile object to CouchDB using the command line arguments.
    
        Arguments:
            options - the options parsed from the command line
    """
    clipped.bulkload2couchdb(options)

def getpolygon(key, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
    """Returns a polygon (list of Points) of the cell defined by the given key.

    Arguments:
        key - the unique identifier for a cell
        cells_per_degree - the desired resolution of the grid
        digits - the number of digits of precision to retain in the coordinates
        a - the semi-major axis of the ellipsoid for the coordinate reference system
        inverse_flattening - the inverse of the ellipsoid's flattening parameter
    """
    return RMGCell.polygon(key, cells_per_degree, digits, a, inverse_flattening)

def translatevariable(varval):
    """ Returns a translated value for a starspan-processed Worldclim variable. Translation includes
        returning spurious large integer values (starspan must use an unsigned short int at some point)
        to their correct negative values, then truncating these values to be integers.
        
        Arguments:
            varval - the variable value as returned from starspan
    """
    newval = float(varval)
    if newval > 55537: # Actual value is a negative number greater than the nodata value of -9999
        newval = newval - 65536
    newval = truncate(newval,0)
    return newval

#class Variable(object):
#    """An environmental variable backed by a .bil and a .hdr file."""
#
#    def __init__(self, bilfile, hdrfile):
#        """Constructs a Variable.
#
#        Arguments:
#            bilfile - The .bil file path.
#            hdrfile - The .hdr file path.
#        """
#        self.bilfile = bilfile
#        self.hdrfile = hdrfile
#        
#        # Loads xmin, xmax, ymin, and ymax values from the .hdr file:
#        for line in open(hdrfile, 'r'):
#            if line.startswith('MaxX'):
#                self.xmax = int(line.split()[1].strip())
#            elif line.startswith('MinX'):
#                self.xmin = int(line.split()[1].strip())
#            elif line.startswith('MaxY'):
#                self.ymax = int(line.split()[1].strip())
#            elif line.startswith('MinY'):
#                self.ymin = int(line.split()[1].strip())

class Cell(object):
    """ A cell described by a key, a polygon, and a grid resolution defined by
        cells_per_degree.
    """

    def __init__(self, key, polygon, cells_per_degree):
        """Constructs a Cell.

        Arguments:
        """
        self.key = key
        self.polygon = polygon
        self.cells_per_degree = cells_per_degree
        
    def __str__(self):
        return str(self.__dict__)
    
class Tile(object):
    """ A tile defined by a geographic coordinate bounding box and the parameters of 
        a coordinate reference system.
    """

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
            filename - the name of the constructed shapefile for the boundaries of the Tile.
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

    def cellbatch2clip2csv2couchdb(self, cells, options, batchnum):
        """ Creates a shapefile for a numbered batch of cells, clips the cells shapefile by 
            the shapefile created from the intersection of the file provided in the gadm
            parameter and the Tile boundary, creates the csv file using starspan to get the 
            variable value statistics out of the Worldclim tile, and finally uploads the 
            resulting cells to CouchDB.
            
            Arguments:
                cells - the list of cells in the batch with their keys and polygons
                options - the options parsed from the command line
                batchnum - the sequential number of the batch of cells

        """
        t0 = time.time()
        filename = os.path.join(os.path.splitext(self.filename)[0], '%s' % batchnum)
        logging.info('Preparing shapefile %s in cellbatch2clip2csv2couchdb().' % (filename) )
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        for cell in cells:
            w.poly(parts=[cell.polygon])
            w.record(CellKey=cell.key)
        w.save(filename)        
        t1 = time.time()
        logging.info('Cell batch shapefile %s prepared in %s' % (filename, t1-t0))
        clippedfile = Tile.clip2cell('%s.shp' % filename, self.filename)
        csvfile = Tile.statistics2csv(clippedfile, options)
        Tile.csv2couch(csvfile, options)

    def polygon(self):
        """ Returns a closed polygon (list of Points - nw, sw, se, ne, nw) for the Tile."""
        n = float(truncate(self.nwcorner.lat, self.digits))
        s = float(truncate(self.secorner.lat, self.digits))
        w = float(truncate(self.nwcorner.lng, self.digits))
        e = float(truncate(self.secorner.lng, self.digits))
        return [(w, n), (w, s), (e, s), (e, n), (w, n)]

    @classmethod
    def csv2couch(cls, csvfile, options):
        """ Loads cells from csv file to CouchDB.
            
            Arguments:
                csvfile - the CSV file containing the cells to load
                options - the options parsed from the command line. Uses:
                    couchurl - the URL of the CouchDB server
                    cells_per_degree - the desired resolution of the grid
        """
        t0 = time.time()
        logging.info('Beginning csv2couch(), preparing cells for bulkloading from %s.' % (csvfile) )
        server = couchdb.Server(options.couchurl)
        cdb = server[options.database]
        cells_per_degree = float(options.cells_per_degree)
        dr = csv.DictReader(open(csvfile, 'r'))
        cells = {}
        for row in dr:
            cellkey = row.get('CellKey')
            if not cells.has_key(cellkey):
                cells[cellkey] = {
                    '_id': cellkey,
                    'tile': options.key,
                    'coords': getpolygon(cellkey, cells_per_degree),
                    'vars': {}
                    }            
            varname = row.get('RID').split('_')[0]
            # the following dependent on running starspan with --stats %s avg
            varval = row.get('avg_Band1')
            # the following need for Worldclim because of the -9999 NODATA value.
            cells.get(cellkey).get('vars')[varname] = translatevariable(varval)
        t1 = time.time()
        logging.info('%s cells prepared for upload in %s' % (len(cells), t1-t0))
        cdb.update(cells.values())
        t2 = time.time()
        logging.info('%s documents uploaded in %s' % (len(cells), t2-t1))

    @classmethod
    def statistics2csv(cls, shapefile, options):      
        """ Extracts statistics on variables in the Worldclim tile for the cells
            in the shapefile via starspan.

            Arguments:
                shapefile - the shapefile containing the cells
                options - the options parsed from the command line. Uses:
                    vardir - the path to the directory in which to store the Worldclim files
        """
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
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('starspan statistics finished in %s.' % (t1-t0) )
        return csvfile
        
    @classmethod
    def clip2cell(cls, src, shapefile):
        """ Creates a new shapefile that is the intersection of the src file 
            and the given shapefile. Returns the name of the clipped shapefile.
            
            Arguments:
                src - the file to intersect with the shapefile
                shapefile - the shapefile to intersect with the src
        """
        t0 = time.time()
        logging.info('Beginning clipping of %s by %s.' % (src, shapefile) )
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        clipped = src.replace('.shp', '-clipped.shp')
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, shapefile, clipped, src)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('%s clipped by %s in %s' % (src, shapefile, t1-t0))
        return clipped
                
    def bulkload2couchdb(self, options):
        """ Loads the Tile to CouchDB using the command line arguments.
        
            Arguments:
                options - the options parsed from the command line. Uses:
                    batchsize - the number of cells to include in a batch to avoid memory overflow
                    cells_per_degree - the number of cells in a degree of longitude at
                        the equator - defines a resolution for the grid system.
        """
        t0 = time.time()
        logging.info('Beginning bulkload2couchdb().')
        batchsize = int(options.batchsize)
        batchnum = 0
        cells_per_degree = float(options.cells_per_degree)
        cells = []
        count = 0
        for cell in self.getcells():
            cells.append(cell)
            count += 1
            if count >= batchsize:
                self.cellbatch2clip2csv2couchdb(cells, options, batchnum)
                count = 0
                cells = []
                batchnum += 1
                continue
        if count > 0:
            self.cellbatch2clip2csv2couchdb(cells, options, batchnum)
        t1 = time.time()
        logging.info('Total elapsed time to bulkload2couchdb(): %s' % (t1-t0))
    
    def clip(self, shapefile, workspace):
        """ Creates a new shapefile in the workspace that is the intersection of 
            this Tile and the given shapefile. Returns a Tile clipped by shapefile.
        
            Arguments:
                shapefile - the shapefile to clip by the boundaries of the Tile
                workspace - the directory in which to store the clipped shape file
        """
        t0 = time.time()
        ogr2ogr = '/usr/local/bin/ogr2ogr'
        this = self.writetileshapefile(workspace)
        clipped = this.replace('.shp', '-clipped.shp')
        logging.info('Beginning clipping of %s by %s.' % (shapefile, this))
        command = '%s -clipsrc %s %s %s' % (ogr2ogr, this, clipped, shapefile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('%s clipped by %s in %s' % (shapefile, this, t1-t0))
        return Tile(self.key, self.nwcorner, self.secorner, self.cells_per_degree, self.digits, self.a, self.inverse_flattening, clipped)

    def writetileshapefile(self, workspace):
        """ Writes a shapefile ([self.key].dbf, [self.key].shp, [self.key].shx) for boundaries 
            of the Tile.
        
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

#    def writeshapefile(self, workspace):
#        """Writes tile shapefile in workspace directory and returns filename.
#        
#        Arguments:
#            workspace - the directory to store the clipped shape file
#        """
#        cell = self.getcells().next()
#        cellinfo = 'Cell: %s' % (cell)
#        fout = os.path.join(workspace, cell.key)
#        w = shapefile.Writer(shapefile.POLYGON)
#        w.field('CellKey','C','255')
#        w.poly(parts=[cell.polygon])
#        w.record(CellKey=cell.key)
#        w.save(fout)        
#        return '%s.shp' % fout

    def getcells(self):
        """ Yields Cells by iterating west to east, north to south over RMG cells
            within the bounding box of a Tile.
        """ 
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
                yield Cell(key, polygon, self.cells_per_degree)
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
    """ Parses command line options and returns them."""
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
                      help="A clipped variables csv file to statistics2csv and load.",
                      default=None)
    parser.add_option("-u", 
                      "--couchurl", 
                      dest="couchurl",
                      help="The CouchDB URL.",
                      default=None)
    parser.add_option("-d", 
                      "--database", 
                      dest="database",
                      help="The CouchDB database name.",
                      default=None)
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
    parser.add_option("-l", 
                      "--logfile", 
                      dest="logfile",
                      help="The name of the log file",
                      default=None)
    return parser.parse_args()[0]

if __name__ == '__main__':
    options = _getoptions()
    command = options.command.lower()

    if options.logfile:
        logfile = os.path.join(options.workspace, options.logfile)
    else:
        logfilename = 'sdl-%s-%s-%s.log' % (options.command, options.key, str(int(time.time())))
        logfile = os.path.join(options.workspace, logfilename)

    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=logfile,
                    filemode='w')    

    if command == 'clip':
        logging.info('Beginning command clip.')
        t0 = time.time()
        clipped = clip(options)
        t1 = time.time()
        logging.info('Finished command clip in %ss.' % (t1-t0))

    if command == 'csv2couchdb':
        logging.info('Beginning command csv2couch.')
        t0 = time.time()
        tile = maketile(options)
        csvfile = os.path.join(options.workspace, options.csvfile)
        tile.csv2couch(csvfile, options)
        t1 = time.time()
        logging.info('Finished command csv2couch in %ss.' % (t1-t0))

    if command == 'load':
        logging.info('Beginning command load.')
        t0 = time.time()
        clipped = clip(options)
        load(options, clipped)    
        t1 = time.time()
        logging.info('Finished command load in %ss.' % (t1-t0))

    if command == 'getworldclimtile':
        logging.info('Beginning command getworldclimtile.')
        t0 = time.time()
        getworldclimtile(options)
        t1 = time.time()
        logging.info('Finished command getworldclimtile in %ss.' % (t1-t0))

    if command == 'full':
        logging.info('Beginning command full.')
        t0 = time.time()
        getworldclimtile(options)
        clipped = clip(options)
        load(options, clipped)    
        t1 = time.time()
        logging.info('Finished command full in %ss.' % (t1-t0))