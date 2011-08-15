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

'''
This module includes classes and a command line interface for bulkloading 
WorldClim environment variables to CouchDB using the Rectangular Mesh Grid (RMG).
'''

import csv
import re
import glob
import simplejson
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
import sys
from rmg import *

VARDICT = {'tmean':'t', 'tmin':'m', 'tmax':'x', 'alt':'a', 'bio':'b', 'prec':'p'}
#OGR2OGR='/Library/Frameworks/GDAL.framework/Programs/ogr2ogr'
OGR2OGR = '/usr/local/bin/ogr2ogr'

#class Variable(object):
#    '''An environmental variable backed by a .bil and a .hdr file.'''
#
#    def __init__(self, bilfile, hdrfile):
#        '''Constructs a Variable.
#
#        Arguments:
#            bilfile - The .bil file path.
#            hdrfile - The .hdr file path.
#        '''
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
    ''' A cell described by a key, a polygon, and a grid resolution defined by
        cells_per_degree.
    '''

    def __init__(self, key, polygon, cells_per_degree):
        '''Constructs a Cell.

        Arguments:
        '''
        self.key = key
        self.polygon = polygon
        self.cells_per_degree = cells_per_degree
        
    def __str__(self):
        return str(self.__dict__)
    
class Tile(object):
    ''' A tile defined by a geographic coordinate bounding box and the parameters of 
        a coordinate reference system.
    '''

    def __init__(self, key, nwcorner, secorner, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING, filename=None):
        '''Tile constructor.

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
        '''
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

    def polygon(self):
        ''' Returns a closed polygon (list of Points - nw, sw, se, ne, nw) for the Tile.'''
        n = float(truncate(self.nwcorner.lat, self.digits))
        s = float(truncate(self.secorner.lat, self.digits))
        w = float(truncate(self.nwcorner.lng, self.digits))
        e = float(truncate(self.secorner.lng, self.digits))
        return [(w, n), (w, s), (e, s), (e, n), (w, n)]

    def getcells(self):
        ''' Yields Cells by iterating west to east, north to south over RMG cells
            within the bounding box of a Tile.
        ''' 
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

    def batchprocesscells(self, batchdir, options):
        ''' Iterate through the cells of the Tile, creating shape files of cells in batches.
        
            Arguments:
                batchdir - the directory in which to put the cell batch files
                options - the options parsed from the command line. Uses:
                    batchsize - the number of cells to include in a batch to avoid memory overflow
                    cells_per_degree - the number of cells in a degree of longitude at
                        the equator. Defines a resolution for the grid system.
        '''
        batchsize = int(options.batchsize)
        batchnum = 0
        cells_per_degree = float(options.cells_per_degree)
        cells = []
        count = 0
        for cell in self.getcells():
            cells.append(cell)
            count += 1
            if count >= batchsize:
                self.cellbatch2shapefile(cells, batchnum, batchdir, options)
                count = 0
                cells = []
                batchnum += 1
                continue
        if count > 0:
            self.cellbatch2shapefile(cells, batchnum, batchdir, options)
        return True

    def cellbatch2shapefile(self, cells, batchnum, batchdir, options):
        ''' Creates a shapefile for a numbered batch of cells.
            
            Arguments:
                cells - the list of cells in the batch with their keys and polygons
                batchnum - the sequential number of the batch of cells
                batchdir - the directory in which to put the cell batch files
                options - the options parsed from the command line
        '''
        filename = os.path.join(batchdir, '%s' % batchnum)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('CellKey','C','255')
        for cell in cells:
            w.poly(parts=[cell.polygon])
            w.record(CellKey=cell.key)
        w.save(filename)        
        return filename

    def clip(self, shapefile, workspace):
        ''' Creates a new shapefile in the workspace that is the intersection of 
            this Tile and the given shapefile. Returns a Tile clipped by shapefile.
        
            Arguments:
                shapefile - the shapefile to clip by the boundaries of the Tile
                workspace - the directory in which to store the clipped shape file
        '''
        t0 = time.time()
        this = self.writetileshapefile(workspace)
        clipped = this.replace('.shp', '-clipped.shp')
        logging.info('Beginning clipping of %s by %s.' % (shapefile, this))
        command = '%s -clipsrc %s %s %s' % (OGR2OGR, this, clipped, shapefile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        t1 = time.time()
        logging.info('%s clipped by %s in %s' % (shapefile, this, t1-t0))
        return Tile(self.key, self.nwcorner, self.secorner, self.cells_per_degree, self.digits, self.a, self.inverse_flattening, clipped)

    def writetileshapefile(self, workspace):
        ''' Writes a shapefile ([self.key].dbf, [self.key].shp, [self.key].shx) for boundaries 
            of the Tile.
        
            Arguments:
                workspace - the directory to store the clipped shape file
        '''
        fout = os.path.join(workspace, self.key)
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('TileKey','C','255')
        w.poly(parts=[self.polygon()])
        w.record(TileKey=self.key)
        w.save(fout)        
        return '%s.shp' % fout

#    @classmethod
#    def clip2cell(cls, src, shapefile):
#        ''' Creates a new shapefile that is the intersection of the src file 
#            and the given shapefile. Returns the name of the clipped shapefile.
#            
#            Arguments:
#                src - the file to intersect with the shapefile
#                shapefile - the shapefile to intersect with the src
#        '''
#        t0 = time.time()
#        logging.info('Beginning clipping of %s by %s.' % (src, shapefile) )
#        clipped = src.replace('.shp', '-clipped.shp')
#        command = '%s -clipsrc %s %s %s' % (OGR2OGR, shapefile, clipped, src)
#        args = shlex.split(command)
#        subprocess.call(args)
#        t1 = time.time()
#        logging.info('%s clipped by %s in %s' % (src, shapefile, t1-t0))
#        return clipped
#                
#    def bulkload2couchdb(self, options):
#        ''' Loads the Tile to CouchDB using the command line arguments.
#        
#            Arguments:
#                options - the options parsed from the command line. Uses:
#                    batchsize - the number of cells to include in a batch to avoid memory overflow
#                    cells_per_degree - the number of cells in a degree of longitude at
#                        the equator - defines a resolution for the grid system.
#        '''
#        t0 = time.time()
#        logging.info('Beginning bulkload2couchdb().')
#        batchsize = int(options.batchsize)
#        batchnum = 0
#        cells_per_degree = float(options.cells_per_degree)
#        cells = []
#        count = 0
#        for cell in self.getcells():
#            cells.append(cell)
#            count += 1
#            if count >= batchsize:
#                self.cellbatch2clip2csv2couchdb(cells, options, batchnum)
#                count = 0
#                cells = []
#                batchnum += 1
#                continue
#        if count > 0:
#            self.cellbatch2clip2csv2couchdb(cells, options, batchnum)
#        t1 = time.time()
#        logging.info('Total elapsed time to bulkload2couchdb(): %s' % (t1-t0))
#    def cellbatch2clip2csv2couchdb(self, cells, options, batchnum):
#        ''' Creates a shapefile for a numbered batch of cells, clips the cells shapefile by 
#            the shapefile created from the intersection of the file provided in the gadm
#            parameter and the Tile boundary, creates the csv file using starspan to get the 
#            variable value statistics out of the Worldclim tile, and finally uploads the 
#            resulting cells to CouchDB.
#            
#            Arguments:
#                cells - the list of cells in the batch with their keys and polygons
#                options - the options parsed from the command line
#                batchnum - the sequential number of the batch of cells
#
#        '''
#        t0 = time.time()
#        filename = os.path.join(os.path.splitext(self.filename)[0], '%s' % batchnum)
#        logging.info('Preparing shapefile %s in cellbatch2clip2csv2couchdb().' % (filename) )
#        w = shapefile.Writer(shapefile.POLYGON)
#        w.field('CellKey','C','255')
#        for cell in cells:
#            w.poly(parts=[cell.polygon])
#            w.record(CellKey=cell.key)
#        w.save(filename)        
#        t1 = time.time()
#        logging.info('Cell batch shapefile %s prepared in %s' % (filename, t1-t0))
#        clippedfile = Tile.clip2cell('%s.shp' % filename, self.filename)
#        csvfile = Tile.statistics2csv(clippedfile, options)
#        Tile.starspan2csv2couch(csvfile, options)

#    @classmethod
#    def starspan2csv2couch_stats(cls, csvfile, options): # unused, but save in case stats other than avg are desired at some point
#        ''' Loads cells having avg stdev min and max values of a variable from 
#            starspan2csv file to CouchDB.
#            
#            Arguments:
#                stats - a list of the statistics calculations expected in the starspan file
#                csvfile - the CSV file containing the cells to load
#                options - the options parsed from the command line. Uses:
#                    couchurl - the URL of the CouchDB server
#                    database - the name of the database in CouchDB
#                    cells_per_degree - the desired resolution of the grid
#        '''
#        t0 = time.time()
#        logging.info('Beginning starspan2csv2couch(), preparing cells for bulkloading from %s.' % (csvfile) )
#        server = couchdb.Server(options.couchurl)
#        cdb = server[options.database]
#        stats = ['avg']
##        stats = ['avg','stdev','min','max']
#        cells_per_degree = float(options.cells_per_degree)
#        dr = csv.DictReader(open(csvfile, 'r'))
#        cells = {}
#        for row in dr:
#            cellkey = row.get('CellKey')
#            if not cells.has_key(cellkey):
#                cells[cellkey] = {
#                    '_id': cellkey,
#                    'tile': options.key,
#                    'coords': getpolygon(cellkey, cells_per_degree),
#                    'vars': {}
#                    }            
#            varname = row.get('RID').split('_')[0]
#            for stat in stats:
#                statname = '%s_Band1' % (stat)
#                statval = row.get(statname)
#                varstat = '%s' % (varname)
#                # removed capacity to process multiple stats - just getting avg now.
##                varstat = '%s_%s' % (varname,stat)
#                # the following needed to process Worldclim stats for representation in the datastore
#                cells.get(cellkey).get('vars')[varstat] = translatestatistic(varname, stat, statval)
#        t1 = time.time()
#        logging.info('%s cells prepared for upload in %s' % (len(cells), t1-t0))
#        cdb.update(cells.values())
#        t2 = time.time()
#        logging.info('%s documents uploaded in %s' % (len(cells), t2-t1))

#    @classmethod
#    def starspan2csv2couch(cls, csvfile, options):
#        ''' Loads cells having avg values of a variable from starspan2csv file to CouchDB.
#            
#            Arguments:
#                csvfile - the CSV file containing the cells to load
#                options - the options parsed from the command line. Uses:
#                    couchurl - the URL of the CouchDB server
#                    database - the name of the database in CouchDB
#                    cells_per_degree - the desired resolution of the grid
#        '''
#        t0 = time.time()
#        logging.info('Beginning starspan2csv2couch(), preparing cells for bulkloading from %s.' % (csvfile) )
#        server = couchdb.Server(options.couchurl)
#        cdb = server[options.database]
#        cells_per_degree = float(options.cells_per_degree)
#        drf = csv.DictReader(open(csvfile, 'r'))
#        firstkey = drf.next().get('CellKey')
#        cells = {}
#        dr = csv.DictReader(open(csvfile, 'r'))
#        for row in dr:
#            cellkey = row.get('CellKey')
#            if not cells.has_key(cellkey):
#                cells[cellkey] = {
#                    '_id': cellkey,
#                    't': options.key,  # tile number
#                    'p': getpolygon(cellkey, cells_per_degree), # polygon 
#                    'v': {}  # dictionary of variables (tmean1, ... tmin1, ... tmax1, ...alt, ... prec1, ... bio1, ...
#                    }
##                cells[cellkey] = {
##                    '_id': cellkey,
##                    'tile': options.key,
##                    'coords': getpolygon(cellkey, cells_per_degree),
##                    'vars': {}
##                    }
#            varname = row.get('RID').split('_')[0]
#            varalias=getvaralias(varname)
#            cells.get(cellkey).get('vars')[varalias] = translatevariable(row.get('avg_Band1'))
#            lastkey=cellkey
#        t1 = time.time()
#        logging.info('%s cells prepared for upload in %s' % (len(cells), t1-t0))
#### Write out to workspace tileno_batchno.csv
#        newcells = dict((k, simplejson.dumps(v)) for k,v in cells.iteritems())
#        outfile = '%s_%s-%s.csv' % (options.key,firstkey,lastkey)
#        os.path.join(options.workspace, outfile)
#        dw = csv.DictWriter(open(outcsvfile,'w'))
#        dw.writerow(cells)
#            
#        cdb.update(cells.values())
#        t2 = time.time()
#        logging.info('%s documents uploaded in %s' % (len(cells), t2-t1))

#    @classmethod
#    def statistics2csv(cls, shapefile, options):
#        ''' Extracts statistics on variables in the Worldclim tile for the cells
#            in the shapefile via starspan.
#
#            Arguments:
#                shapefile - the shapefile containing the cells
#                options - the options parsed from the command line. Uses:
#                    vardir - the path to the directory in which to store the Worldclim files
#        '''
#        t0 = time.time()
#        logging.info('Beginning starspan statistics on %s.' % (shapefile) )
#        variables = [os.path.join(options.vardir, x) \
#                         for x in os.listdir(options.vardir) \
#                         if x.endswith('.bil')]
#        variables = reduce(lambda x,y: '%s %s' % (x, y), variables)
#        csvfile = shapefile.replace('.shp', '.csv')
#        # Call starspan requesting mean of variable, excluding nodata values (-9999 in the file is the same as 55537)
#        # starspan -- vector 0-clipped.shp --raster tmean
#        command = 'starspan --vector %s --raster %s --stats %s avg --nodata 55537' \
#            % (shapefile, variables, csvfile)
##        command = 'starspan --vector %s --raster %s --stats %s avg stdev min max --nodata 55537' \
##            % (shapefile, variables, csvfile)
#        args = shlex.split(command)
#        subprocess.call(args)
#        t1 = time.time()
#        logging.info('starspan statistics finished in %s.' % (t1-t0) )
#        return csvfile
        
#    @classmethod
#    def starspan2csv(cls, shapefile, options):      
#        ''' Extracts statistics on variables in the Worldclim tile for the cells
#            in the shapefile via starspan.
#
#            Arguments:
#                shapefile - the shapefile containing the cells
#                options - the options parsed from the command line. Uses:
#                    vardir - the path to the directory in which to store the Worldclim files
#        '''
#        t0 = time.time()
#        logging.info('Beginning starspan statistics on %s.' % (shapefile) )
#        variables = [os.path.join(options.vardir, x) \
#                         for x in os.listdir(options.vardir) \
#                         if x.endswith('.bil')]
#        variables = reduce(lambda x,y: '%s %s' % (x, y), variables)
#        csvfile = shapefile.replace('.shp', '.csv')
#        # Call starspan requesting mean of variable, excluding nodata values (-9999 in the file is the same as 55537)
#        command = 'starspan --vector %s --raster %s --stats %s avg --nodata 55537' \
#            % (shapefile, variables, csvfile)
#        args = shlex.split(command)
#        subprocess.call(args)
#        t1 = time.time()
#        logging.info('starspan statistics finished in %s.' % (t1-t0) )
#        return csvfile
        
    
def getvaralias(varname):
    varcategory = re.split('[0-9]+',varname)[0]
    newvarname=VARDICT[varcategory]
    if varname != 'alt':
        varnum = re.split('[a-z]+',varname)[1]
        newvarname = '%s%s' % (newvarname, varnum)
    return newvarname

def maketile(options):
    ''' Creates a Tile object from the command line arguments.
    
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                nwcorner - the northwest corner of the Tile
                secorner - the southeast corner of the Tile
                cells_per_degree - the number of cells in a degree of longitude at
                    the equator - defines a resolution for the grid system.
    '''
    key = options.key
    nw = map(float, options.nwcorner.split(','))
    se = map(float, options.secorner.split(','))
    nwcorner = Point(nw[0], nw[1])
    secorner = Point(se[0], se[1])
    cells_per_degree = float(options.cells_per_degree)
    tile = Tile(key, nwcorner, secorner, cells_per_degree)
    return tile

def makeclippedtile(options):
    ''' Returns a Tile object constructed from the command line arguments. Resulting
        Tile is clipped to the given bounding box as well as to the clipping layer
        provided in the -gadm argument. Bounding box for the area of interest must 
        be in or on the boundaries of the Tile.
    
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                nwcorner - the northwest corner of the area of interest within the Tile
                secorner - the southeast corner of the area of interest within the Tile
                cells_per_degree - the number of cells in a degree of longitude at
                    the equator. Defines a resolution for the grid system.
                gadm - the shapefile to clip to (polygon of are containing Worldclim data)
    '''
    tile = maketile(options)
    return tile.clip(options.gadm, options.workspace)

def getpolygon(key, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
    '''Returns a polygon (list of Points) of the cell defined by the given key.

    Arguments:
        key - the unique identifier for a cell
        cells_per_degree - the desired resolution of the grid
        digits - the number of digits of precision to retain in the coordinates
        a - the semi-major axis of the ellipsoid for the coordinate reference system
        inverse_flattening - the inverse of the ellipsoid's flattening parameter
    '''
    return RMGCell.polygon(key, cells_per_degree, digits, a, inverse_flattening)

def getboundingbox(key, cells_per_degree, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
    '''Returns a boundingbox (list of Points) of the cell defined by the given key.

    Arguments:
        key - the unique identifier for a cell
        cells_per_degree - the desired resolution of the grid
        digits - the number of digits of precision to retain in the coordinates
        a - the semi-major axis of the ellipsoid for the coordinate reference system
        inverse_flattening - the inverse of the ellipsoid's flattening parameter
    '''
    return RMGCell.boundingbox(key, cells_per_degree, digits, a, inverse_flattening)

def translatestatistic(varname, stat, statval):
    ''' Returns a translated value for a starspan-processed Worldclim statistic (avg, stdev, min, max). 
        Translation includes returning spurious large integer values 
        (starspan must use an unsigned short int at some point)
        to their correct negative values, then truncating these values to be integers.
        
        Arguments:
            varname - the name of the variable to which the statistic applies (tmean, prec, alt...)
            stat - the statistic the value represents (avg, stdev, min, max)
            statval - the value of the statistic as returned from starspan
    '''
#    logging.info('VARNAME: %s STAT: %s STATVAL: %s' % (varname, stat, statval))
    newval = float(statval)
    if stat == 'stdev':
        return truncate(newval,0)
    else:
        if newval > 55537: # Actual value is a negative number greater than the nodata value of -9999
            newval = newval - 65536
        if varname == 'tmean':
            # tmen in Worldclim is degrees Celsius times 10, turn it into degrees Celsius
            newval = newval/10.0
            return truncate(newval,1)
    return truncate(newval,0)

def translatevariable(varval):
    ''' Returns a translated value for a starspan-processed Worldclim variable. Translation includes
        returning spurious large integer values (starspan must use an unsigned short int at some point)
        to their correct negative values, then truncating these values to be integers.
        
        Arguments:
            varval - the variable value as returned from starspan
    '''
    newval = float(varval)
    if newval > 55537: # Actual value is a negative number greater than the nodata value of -9999
        newval = newval - 65536
    newval = truncate(newval,0)
    return newval

#def clip(options):
#    ''' Creates a clipped Tile object and associated files from the command line arguments.
#        Returns a clipped Tile object.
#    
#        Arguments:
#            options - the options parsed from the command line. Uses:
#                key - the tile identifier (e.g., '37' for Worldclim tile 37)
#                nwcorner - the northwest corner of the Tile
#                secorner - the southeast corner of the Tile
#                cells_per_degree - the number of cells in a degree of longitude at
#                    the equator - defines a resolution for the grid system.
#                gadm - the shapefile to clip to (polygon of are containing Worldclim data)
#    '''
#    key = options.key
#    nw = map(float, options.nwcorner.split(','))
#    se = map(float, options.secorner.split(','))
#    nwcorner = Point(nw[0], nw[1])
#    secorner = Point(se[0], se[1])
#    cells_per_degree = float(options.cells_per_degree)
#    tile = Tile(key, nwcorner, secorner, cells_per_degree)
#    clipped = tile.clip(options.gadm, options.workspace)
#    return clipped

#def cleanworkspace(options):
#    if os.path.exists(options.workspace):
#        base = '%s-clipped' % options.key
#        basedir = os.path.join(options.workspace, base)
#        if os.path.exists(basedir):
#            basedirall = os.path.join(options.workspace, options.key) 
#            command = 'rm %s/*' % (basedir)
#            logging.info(command)
#            args = shlex.split(command)
#            subprocess.call(args)
#
#        base = os.path.join(options.workspace, options.key)
#        command = 'rm %s*' % (base)
#        logging.info(command)
#        args = shlex.split(command)
#        subprocess.call(args)
#    else:
#        command = 'mkdir %s-clipped' % (base)
#        logging.info(command)
#        args = shlex.split(command)
#        subprocess.call(args)
#
#def load(options, clipped):
#    ''' Loads a clipped Tile object to CouchDB using the command line arguments.
#    
#        Arguments:
#            options - the options parsed from the command line
#    '''
#    cleanworkspace(options)
#    clipped.bulkload2couchdb(options)

def prepareworkspace(options):
    if not os.path.exists(options.workspace):
        os.mkdir(options.workspace)
        if not os.path.exists(options.workspace):
            loccing.ingo('Unable to make workspace at %s.' % cleanpath)
            return False
    if len(os.listdir(options.workspace))>0:
        logging.info('Workspace %s not empty. Empty before proceeding. Make sure directory is not in use by other processes.' % options.workspace)
        return False
    return True

def getworldclimtile(options):
    ''' Downloads the files for all of the variables for the Tile given by the Tile key
        in the command line, unzips them, and removes the zipfile.
        
        Arguments:
            options - the options parsed from the command line. Uses:
                key - the tile identifier (e.g., '37' for Worldclim tile 37)
                vardir - the path to the directory in which to store the Worldclim files
    '''
    if not os.path.exists(options.vardir):
        os.mkdir(options.vardir)
        if not os.path.exists(options.vardir):
            logging.info('Unable to create variable directory %s.' % options.vardir)
            return False
    ''' Assume that if there are files in this folfer already, it's OK to proceed.'''
    if len(os.listdir(options.vardir))>0:
        return True

    varset = ['tmean','tmin','tmax','prec','alt','bio']
    for var in varset:
        varfile = '%s_%s.zip' % (var, options.key)
        varpath = os.path.join(options.vardir, varfile)
        ''' Get the file from Worldclim site.'''
        command = 'wget -P %s http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/%s' % (options.vardir, varfile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        ''' Unzip the variable zip file.'''
        command = 'unzip %s -d %s' % (varpath, options.vardir)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        ''' Remove the zip file.'''
        os.remove(varpath)
    return True

### starspan outputs to Couch- and GAE-ready CSVs ###
def starspancsvdir2couchcsvs(batchdir, key, cells_per_degree, workspace):
    if not os.path.exists(batchdir):
        logging.info('Path %s does not exist.')
        return
    os.chdir(batchdir)
    for f in glob.glob("*.csv"):
        starspancsv2couchcsv(f,key,cells_per_degree,workspace)

def starspancsv2couchcsv(csvfile, key, cells_per_degree, workspace):
    ''' Processes the avg statistic from a starspan file to a csv file
        having a cellkey and a document - suitable for loading to CouchDB
        or GAE. Returns the full path to the output file.
        
        Arguments:
            csvfile - the CSV file containing the cells to load
            cells_per_degree - the desired resolution of the grid as float
            key - the tile identifier (e.g., '37' for Worldclim tile 37)
            workspace - the directory in which to store the clipped shape file
    '''
    cells_per_degree = float(cells_per_degree)
    drf = csv.DictReader(open(csvfile, 'r'))
    firstkey = drf.next().get('CellKey')
    cells = {}
    dr = csv.DictReader(open(csvfile, 'r'))
    for row in dr:
        cellkey = row.get('CellKey')
        if not cells.has_key(cellkey):
            cells[cellkey] = {
                '_id': cellkey,
                't': key,  # tile number
                'b': getboundingbox(cellkey, cells_per_degree), # boundingbox 
                'v': {}  # dictionary of variables (tmean1,...tmin1,...tmax1,...alt,...prec1,...bio1,...)
                }
        varname = row.get('RID').split('_')[0]
        varalias=getvaralias(varname)
        cells.get(cellkey).get('v')[varalias] = translatevariable(row.get('avg_Band1'))
        lastkey=cellkey
    ''' Remove cells having all zero values. This is true if the value of bio7 (temperature annual range) is 0.'''
    for k in cells.keys():
        if cells[k].get('v')['b7']=='0':
            del(cells[k])
    newcells = [dict(cellkey=k, doc=simplejson.dumps(cells[k])) for k in cells.keys()]
    outfilename = '%s_%s_%s.csv' % (key,firstkey,lastkey)
    couchdir = os.path.join(workspace,'forcouch')
    if not os.path.exists(couchdir):
        os.mkdir(couchdir)
    outfile = os.path.join(couchdir,outfilename)
    fieldnames=['cellkey','doc']
    with open(outfile,'w') as f:
        dw = csv.DictWriter(f, fieldnames, quotechar="'")
        dw.writer.writerow(dw.fieldnames)
        dw.writerows(newcells)
    return outfile

### CSVs to CouchDB ###
def couchcsvdir2couchdb(couchdir, couchurl, database):
    if not os.path.exists(couchdir):
        os.mkdir(couchdir)
    os.chdir(couchdir)
    for f in glob.glob("*.csv"):
        if not couchcsv2couchdb(f,couchurl,database):
            logging.info('Failed to load %s to %s database %s.' % (f, couchurl, database))

def couchcsv2couchdb(csvfile, couchurl, database):
    ''' Loads cells from a csv file having cellkey and doc properties to CouchDB.
        
        Arguments:
            csvfile - the CSV file containing the cells to load
            couchurl - the URL of the CouchDB server
            database - the name of the database in CouchDB
    '''
    server = couchdb.Server(couchurl)
    cdb = server[database]
    drf = csv.DictReader(open(csvfile, 'r'), quotechar="'")
    cells={}
    for row in drf:
        cellkey=row.get('cellkey')
        doc = row.get('doc')
        cells[cellkey]=simplejson.loads(doc)
    cdb.update(cells.values())
    return True

### Clip cell batch shape files ###
def clipcellbatchfiles(batchdir, cliptoshapefile):
    if not os.path.exists(batchdir):
        logging.info('Unable to find cell batch directory %s.' % batchdir)
        return False
    os.chdir(batchdir)
    for f in glob.glob("*.shp"):
        clipcellbatchshapefile(f,cliptoshapefile)
    return True
    
def clipcellbatchshapefile(cellbatchshapefile, cliptoshapefile):
    ''' Creates a new shapefile that is the intersection of the src file 
        and the given shapefile. Returns the name of the clipped shapefile.
        
        Arguments:
            cellbatchshapefile - the file to intersect with the shapefile
            cliptoshapefile - the shapefile to intersect with the src
    '''
    clipped = cellbatchshapefile.replace('.shp', '-clipped.shp')
    command = '%s -clipsrc %s %s %s' % (OGR2OGR, cliptoshapefile, clipped, cellbatchshapefile)
    logging.info(command)
    args = shlex.split(command)
    subprocess.call(args)
    return clipped

### starspan to CSVs ###
def starspanclippedcellbatchfiles(batchdir, vardir):
    if not os.path.exists(batchdir):
        logging.info('Unable to find cell batch directory %s.' % batchdir)
        return False
    os.chdir(batchdir)
    for f in glob.glob("*-clipped.shp"):
        cellbatchshapefile2statscsv(f,vardir)
    return True

def cellbatchshapefile2statscsv(shapefile, vardir):      
    ''' Extracts statistics on variables in the Worldclim tile for the cells
        in the shapefile via starspan.

        Arguments:
            shapefile - the shapefile containing the cells
            vardir - the path to the directory in which to store the Worldclim files
    '''
    variables = [os.path.join(vardir, x) \
                     for x in os.listdir(vardir) \
                     if x.endswith('.bil')]
    variables = reduce(lambda x,y: '%s %s' % (x, y), variables)
    csvfile = shapefile.replace('.shp', '.csv')
    # Call starspan requesting mean of variable, excluding nodata values (-9999 in the file is the same as 55537)
    # starspan -- vector 0-clipped.shp --raster tmean
    command = 'starspan --vector %s --raster %s --stats %s avg --nodata 55537' \
        % (shapefile, variables, csvfile)
#        command = 'starspan --vector %s --raster %s --stats %s avg stdev min max --nodata 55537' \
#            % (shapefile, variables, csvfile)
    args = shlex.split(command)
    subprocess.call(args)
    return csvfile

def _getoptions():
    ''' Parses command line options and returns them.'''
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
    parser.add_option("-o", 
                      "--outputpath", 
                      dest="outputpath",
                      help="The directory into which to place files ready for uploading to data stores",
                      default=None)
    parser.add_option("--csv_dir", 
                      type='string',
                      dest="csv_dir",
                      help="Directory of CSV file outputs from StartSpan.",
                      default=None)

    parser.add_option('--config_file', 
                      type='string', 
                      dest='config_file',
                      metavar='FILE', 
                      help='Bulkload YAML config file.')
    
    parser.add_option('--filename', 
                      type='string', 
                      dest='filename',
                      metavar='FILE', 
                      help='CSV file with data to bulkload.')                      
    
    parser.add_option('--url', 
                      type='string', 
                      dest='url',
                      help='URL endpoint to /remote_api to bulkload to.')                          

    return parser.parse_args()[0]

def main():
    options = _getoptions()
    command = options.command.lower()
    if options.logfile:
        if options.logfile == 'none':
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=options.logfile,
                    filemode='w')    
    else:
        logfile = 'sdl-%s-%s-%s.log' % (options.command, options.key, str(int(time.time())))
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=logfile,
                    filemode='w')    

    if command == 'csv2appengine':
        '''Bulkloads all CSV files in a directory to App Engine datastore.'''
        cmd = 'appcfg.py upload_data --batch_size=%s --num_threads=%s --config_file=%s --filename=%s --kind CellIndex --url=%s'
        os.chdir(options.csv_dir)
        for csvfile in os.listdir('.'):
            if not csvfile.endswith('.csv'):
                continue            
            cmd_line = cmd % (
                1,
                5,
                options.config_file, 
                os.path.abspath(csvfile),
                options.url)            
            logging.info(cmd_line)
            args = shlex.split(cmd_line)
            subprocess.call(args)            
        sys.exit(1)

    if command=='prepareworkspace':
        logging.info('Beginning prepareworkspace()...%s' % options.workspace)
        t0 = time.time()
        if not prepareworkspace(options):
            logging.info('Unable to prepare workspace %s.' % options.workspace)
        t1 = time.time()
        logging.info('Total elapsed time to prepareworkspace(): %s' % (t1-t0))
        sys.exit(1)

    if command=='getworldclimtile':
        ''' Make sure Worlclim layers are in place for the Tile.'''
        logging.info('Beginning getworldclimtile()...get Worldclim layers for Tile %s' % options.key)
        t0 = time.time()
        if not getworldclimtile(options):
            logging.ingo('Unable to prepare Worldclim tile layers.')
            sys.exit(0)
        t1 = time.time()
        logging.info('Total elapsed time to getworldclimtile(): %s' % (t1-t0))
        sys.exit(1)

    if command=='batchcellstoshapes':
        batchdir = os.path.join(options.workspace,'batches')
        logging.info('Beginning prepareworkspace()...%s' % options.workspace)
        t0 = time.time()
        if not prepareworkspace(options):
            logging.info('Unable to prepare workspace %s.' % options.workspace)
            sys.exit(0)
        t1 = time.time()
        logging.info('Total elapsed time to prepareworkspace(): %s' % (t1-t0))

        logging.info('Beginning makeclippedtile()...create clipped Tile of area to process')
        t0 = time.time()
        clippedtile = makeclippedtile(options)
        t1 = time.time()
        logging.info('Total elapsed time to makeclippedtile(): %s' % (t1-t0))

        if not os.path.exists(batchdir):
            os.mkdir(batchdir)
            if not os.path.exists(batchdir):
                logging.info('Unable to make batch directory %s.' % batchdir) 
                sys.exit(0)

        logging.info('Beginning batchprocesscells()...')
        t0 = time.time()
        clippedtile.batchprocesscells(batchdir, options)
        t1 = time.time()
        logging.info('Total elapsed time to batchprocesscells(): %s' % (t1-t0))
        
        logging.info('Beginning clipcellbatchfiles()...')
        t0 = time.time()
        clippingshape = os.path.join(options.workspace,clippedtile.filename)
        clipcellbatchfiles(batchdir, clippingshape)
        t1 = time.time()
        logging.info('Total elapsed time to clipcellbatchfiles(): %s' % (t1-t0))
        sys.exit(1)

    if command=='starspan':
        batchdir = os.path.join(options.workspace,'batches')
        
        logging.info('Beginning starspanclippedcellbatchfiles()...')
        t0 = time.time()
        starspanclippedcellbatchfiles(batchdir, options.vardir)
        t1 = time.time()
        logging.info('Total elapsed time to starspanclippedcellbatchfiles(): %s' % (t1-t0))
        sys.exit(1)

    if command=='starspan2couch':
        batchdir = os.path.join(options.workspace,'batches')

        logging.info('Beginning starspancsvdir2couchcsvs()...')
        t0 = time.time()
        starspancsvdir2couchcsvs(batchdir, options.key, options.cells_per_degree, options.workspace)
        t1 = time.time()
        logging.info('Total elapsed time to starspancsvdir2couchcsvs(): %s' % (t1-t0))
        sys.exit(1)
                
    if command=='couchfromcsvs':
        couchdir = os.path.join(options.workspace,'forcouch')
        
        logging.info('Beginning couchcsvdir2couchdb()...')
        t0 = time.time()
        couchcsvdir2couchdb(couchdir, options.couchurl, options.database)
        t1 = time.time()
        logging.info('Total elapsed time to couchcsvdir2couchdb(): %s' % (t1-t0))
        sys.exit(1)

    if command=='full':
#./sdl.py -n 120 -u http://eighty.berkeley.edu:5984 -d test -g /Users/tuco/Projects/Spatial-Data-Library/data/gadm/Terrestrial-10min-unbuffered-dissolved.shp -b 50000 -k 37 -v /Users/tuco/Data/SDL/worldclim/37 -w /Users/tuco/Data/SDL/workspace/Tile37-1 -f 30,0 -t 31,-1 -c full -l none
#./sdl.py -n 120 -u http://eighty.berkeley.edu:5984 -d test -g /home/tuco/SDL/Spatial-Data-Library/data/gadm/Terrestrial-10min-unbuffered-dissolved.shp -b 50000 -k 37 -v /home/tuco/Data/SDL/worldclim/37 -w /home/tuco/Data/SDL/workspace/Tile37-1 -f 30,0 -t 31,-1 -c full -l none
        batchdir = os.path.join(options.workspace,'batches')
        couchdir = os.path.join(options.workspace,'forcouch')

        if not prepareworkspace(options):
            logging.info('Unable to prepare workspace.')
            sys.exit(0)

        ''' Make sure Worlclim layers are in place for the Tile.'''
        logging.info('Beginning getworldclimtile()...get Worldclim layers for Tile %s' % options.key)
        t0 = time.time()
        if not getworldclimtile(options):
            loggging.info('Unable to prepare Worldclim tile layers.')
            sys.exit(0)
        t1 = time.time()
        logging.info('Total elapsed time to getworldclimtile(): %s' % (t1-t0))

        logging.info('Beginning makeclippedtile()...create clipped Tile of area to process')
        t0 = time.time()
        clippedtile = makeclippedtile(options)
        t1 = time.time()
        logging.info('Total elapsed time to makeclippedtile(): %s' % (t1-t0))

        if not os.path.exists(batchdir):
            os.mkdir(batchdir)
            if not os.path.exists(batchdir):
                logging.info('Unable to make batch directory %s.' % batchdir) 
                sys.exit(0)

        logging.info('Beginning batchprocesscells()...')
        t0 = time.time()
        clippedtile.batchprocesscells(batchdir, options)
        t1 = time.time()
        logging.info('Total elapsed time to batchprocesscells(): %s' % (t1-t0))
        
        logging.info('Beginning clipcellbatchfiles()...')
        t0 = time.time()
        clippingshape = os.path.join(options.workspace,clippedtile.filename)
        clipcellbatchfiles(batchdir, clippingshape)
        t1 = time.time()
        logging.info('Total elapsed time to clipcellbatchfiles(): %s' % (t1-t0))

        logging.info('Beginning starspanclippedcellbatchfiles()...')
        t0 = time.time()
        starspanclippedcellbatchfiles(batchdir, options.vardir)
        t1 = time.time()
        logging.info('Total elapsed time to starspanclippedcellbatchfiles(): %s' % (t1-t0))

        logging.info('Beginning starspancsvdir2couchcsvs()...')
        t0 = time.time()
        starspancsvdir2couchcsvs(batchdir, options.key, options.cells_per_degree, options.workspace)
        t1 = time.time()
        logging.info('Total elapsed time to starspancsvdir2couchcsvs(): %s' % (t1-t0))
        
        logging.info('Beginning couchcsvdir2couchdb()...')
        t0 = time.time()
        couchcsvdir2couchdb(couchdir, options.couchurl, options.database)
        t1 = time.time()
        logging.info('Total elapsed time to couchcsvdir2couchdb(): %s' % (t1-t0))
        sys.exit(1)

    logging.info('Command %s not understood.' % command)
#    if command == 'clip':
#        logging.info('Beginning command clip.')
#        t0 = time.time()
#        clipped = clip(options)
#        t1 = time.time()
#        logging.info('Finished command clip in %ss.' % (t1-t0))
#        sys.exit(1)

#    if command == 'starspan2csv2couchdb':
#        logging.info('Beginning command starspan2csv2couch.')
#        t0 = time.time()
#        tile = maketile(options)
#        csvfile = os.path.join(options.workspace, options.csvfile)
#        tile.starspan2csv2couch(csvfile, options)
#        t1 = time.time()
#        logging.info('Finished command starspan2csv2couch in %ss.' % (t1-t0))
#        sys.exit(1)

#    if command == 'load':
#        logging.info('Beginning command load.')
#        t0 = time.time()
#        clipped = clip(options)
#        load(options, clipped)    
#        t1 = time.time()
#        logging.info('Finished command load in %ss.' % (t1-t0))
#        sys.exit(1)

#    if command == 'full':
#        logging.info('Beginning command full.')
#        t0 = time.time()
#        getworldclimtile(options)
#        clipped = clip(options)
#        load(options, clipped)    
#        t1 = time.time()
#        logging.info('Finished command full in %ss.' % (t1-t0))
#        sys.exit(1)

#    if command == 'cleanworkspace':
#        logging.info('Beginning command cleanworkspace.')
#        t0 = time.time()
#        cleanworkspace(options)
#        t1 = time.time()
#        logging.info('Finished command cleanworkspace in %ss.' % (t1-t0))
#        sys.exit(1)

#    if command == 'test':
#        logging.info('Beginning command test.')
#        t0 = time.time()
#        # sdl.py -c test -k 32 -d worldclim-rmg-32 -u http://eightysteele.berkeley.edu:5984 -n 120 -w /Users/tuco/Data/SDL/workspace -p 0-clipped.csv
## ./sdl.py -v /home/tuco/Data/SDL/worldclim/32 -w /home/tuco/SDL/workspace -u http://eighty.berkeley.edu:5984 -d worldclim-rmg-32 -g /home/tuco/SDL/Spatial-Data-Library/data/gadm/Terrestrial-10min-buffered_00833.shp -k 32 -f 60,0 -t 30,-30 -n 120 -b 50000 > /home/tuco/SDL/workspace/tile32starspan.log &
#        tile = maketile(options)
#        tile.starspan2csv2couch(os.path.join(options.workspace, options.csvfile), options)
#        t1 = time.time()
#        logging.info('Finished command test in %ss.' % (t1-t0))
#        sys.exit(1)

if __name__ == '__main__':
    main()
