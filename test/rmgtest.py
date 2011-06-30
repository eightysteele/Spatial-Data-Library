#!/usr/bin/env python

# Copyright 2011 University of California at Berkeley
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

__author__ = "Aaron Steele and John Wieczorek"

"""This module provides unit testing for Rectangular Mesh Grid classes."""

import logging
import sys
import os
import unittest

#sys.path = [os.path.abspath(os.path.realpath('../'))] + sys.path
sys.path.insert(0, '../')

from sdl.rmg import *
#from sdl.rmg import Point
#from sdl.rmg import RMGCell
#from sdl.rmg import RMGCell
#from sdl.rmg import Point

a = 6378137.0 # WGS84 semi-major axis
inverse_flattening = 298.257223563 # WGS84 inverse flattening

class RMGTest(unittest.TestCase):
    def test_distances_per_degree(self):
        lat = 0
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        self.assertEqual(lngdpd,111319.49079327358)
        self.assertEqual(latdpd,110574.27582159436)

    def test_column_count(self):
        lat = 0
        cells_per_degree = 1
        columns = RMGCell.column_count(lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(columns,358)
        cells_per_degree = 0.01
        columns = RMGCell.column_count(lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(columns,3)

    def test_key0_0(self):
        lng = -180
        lat = 90
        cells_per_degree = 1
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,"0-0")

    def test_key0_1(self):
        lng = -180
        cells_per_degree = 1
        lat = 90 - 1.0/cells_per_degree
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,"0-1")

    def test_key1_1(self):
        lng = 180
        cells_per_degree = 1
        lat = 90 - 1.0/cells_per_degree
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,"6-1")

    def test_center(self):
        lat = 0
        lng = -180
        cells_per_degree = 1
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        center = RMGCell.center(key, cells_per_degree, a, inverse_flattening)
        indexes = key.split('-')
        x_index = int(indexes[0])
        y_index = int(indexes[1])
        mid_lng = RMGCell.mid_lng(x_index, y_index, cells_per_degree, a, inverse_flattening)
        mid_lat = RMGCell.mid_lat(y_index, cells_per_degree)
        self.assertEqual(mid_lng,center.lng)
        self.assertEqual(mid_lat,center.lat)

    def test_tile_count(self):
        # WorldClim Tile 00 cell count
        nwcorner = Point(-180, 90)
        secorner = Point(-150, 60)
        cells_per_degree = 120
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 3334514)
        # WorldClim Tile 37 cell count
        nwcorner = Point(30, 0)
        secorner = Point(60, -30)
        cells_per_degree = 120
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 12308553)
        #Global cell count 120 cells per degree
        nwcorner = Point(-180, 90)
        secorner = Point(180, -90)
        cells_per_degree = 120
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 592726068)
        #Global cell count 1 cells per degree
        nwcorner = Point(-180, 90)
        secorner = Point(180, -90)
        cells_per_degree = 1
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 41246)

    def test_tile(self):
        nwcorner = Point(-180,90)
        secorner = Point(-90,0)
        cells_per_degree = 0.1
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        print tile
        
    def test_tile_kml(self):
        nwcorner = Point(-180,90)
        secorner = Point(-90,0)
        cells_per_degree = 0.1
        digits = 7
        tile = RMGTile(nwcorner, secorner, cells_per_degree, digits)
#        print tile.kml()

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
