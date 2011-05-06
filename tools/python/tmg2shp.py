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

from optparse import OptionParser
from sdl.tmg import gettile
from sdl.tmg import get_rect_tile
import shapefile

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
    parser.add_option("-c", "--cell-count", dest="cell_count",
                      help="Cell count",
                      default=None)
    parser.add_option("-r", "--resolution", dest="resolution",
                      help="Resolution",
                      default=None)
    parser.add_option("-f", "--filename", dest="filename",
                      help="Output path and file name, .shp, .shx, and .dbf will be created.",
                      default=None)

    (options, args) = parser.parse_args()
        
    north = float(options.north)
    south = float(options.south)
    west = float(options.west)
    east = float(options.east)
    resolution = float(options.resolution)
    cell_count = int(options.cell_count)
    filename = options.filename

    w = shapefile.Writer(shapefile.POLYGON)
    w.field('CellKey','C','255')    
#    for x in gettile((west, north), (east, south), resolution, cell_count):
    for x in get_rect_tile((west, north), (east, south), resolution):
        key = x.cellkey
        parts = [list(x) for x in x.polygon]
        w.poly(parts=[parts])
        w.record(CellKey=key)
    w.save(filename)
