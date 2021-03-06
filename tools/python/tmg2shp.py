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
This module supports the creation of ESRI shp files of the Triangular Mesh Grid
overlapping the rectangular grid given by a bounding box and a rectangular grid
resolution in decimal degrees.

The Triangular Mesh Grid design is specified here:
http://goo.gl/0awGK
"""

import fixpath
fixpath.fix_sys_path()

import logging
import math
from optparse import OptionParser
import os
from sdl.tmg import CellPolygon
from sdl.tmg import gettile
from sdl.tmg import get_rect_tile
import shapefile

def getcellkey(globalrow, globalcol, tilerow, tilecol, res, dim):
    """Gets the global cell key."""
    x = int(math.floor((globalcol * (dim / res)) + tilecol))
    y = int(math.floor((globalrow * (dim /res)) + tilerow))
    return '%s-%s' % (x, y)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-n", "--north", dest="north",
                      help="North latitude of bounding box",
                      default=None)
    parser.add_option("-s", "--south", dest="south",
                      help="South latitude of bounding box",
                      default=None)
    parser.add_option("-w", "--west", dest="west",
                      help="West longitude of bounding box",
                      default=None)
    parser.add_option("-e", "--east", dest="east",
                      help="East longitude of bounding box",
                      default=None)
    parser.add_option("-r", "--resolution", dest="resolution",
                      help="Resolution",
                      default=None)
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="Threshold",
                      default=None)
    parser.add_option("-f", "--filename", dest="filename",
                      help="Output path and file name, .shp, .shx, and .dbf will be created.",
                      default=None)
    parser.add_option("-g", "--globaltile", dest="globaltile",
                      help="Global tile x,y",
                      default=None)


    (options, args) = parser.parse_args()

    north = float(options.north)
    south = float(options.south)
    west = float(options.west)
    east = float(options.east)
    y, x = options.globaltile.split(',')
    globalrow = int(y)
    globalcol = int(x)
    resolution = float(options.resolution)
    # Max number of cells per shapefile:
    threshold = int(options.threshold)
    filename = options.filename
    logging.basicConfig(level=logging.DEBUG)

    # Total number of cells to generate:
    ncells = (int(north - south) / resolution) * (int(east-west) / resolution)
    logging.info('Cell count: %s Threshold: %s' % (ncells, threshold) )

    cells = set()
    cellscount = 0
    cellswritten = 0
    filecount = 0
    lng = west
    lat = north
    xindex = 0
    cellkey = '0-0'

    while lng < east:
        yindex = 0
        while lat > south:

            # If at threshold, write all cells to shapefile and flush:
            if cellscount >= threshold:
                logging.info('Writing shapefile')
                w = shapefile.Writer(shapefile.POLYGON)
                w.field('CellKey','C','255')
                for cell in cells:
                    key = cell.cellkey
                    parts = [list(x) for x in cell.polygon]
                    w.poly(parts=[parts])
                    w.record(CellKey=key)
                w.save(os.path.join(filename, cellkey))
                cellswritten += cellscount
                logging.info('Writing %s cells to shapefile %s (%s total cells written)' % (cellscount, cellkey, cellswritten))
                filecount += 1
                cellscount = 0
                cells = set()

            cellkey = getcellkey(globalrow, globalcol, yindex, xindex, resolution, 30) # 30 is tiles size in degrees
            polygon = tuple([
                (lng, lat),
                (lng + resolution, lat),
                (lng + resolution, lat - resolution),
                (lng, lat - resolution),
                (lng, lat)])
            cells.add(CellPolygon(cellkey, polygon))
            cellscount += 1
            lat -= resolution
            yindex += 1
        lat = north
        lng += resolution
        xindex += 1

# Write any remaining cells to shapefile:
if len(cells) > 0:
    w = shapefile.Writer(shapefile.POLYGON)
    w.field('CellKey','C','255')
    for cell in cells:
        key = cell.cellkey
        parts = [list(x) for x in cell.polygon]
        w.poly(parts=[parts])
        w.record(CellKey=key)
    w.save(os.path.join(filename, cellkey))
    logging.info('Writing shapefile %s' % cellkey)
