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
digits = 7

class RMGTest(unittest.TestCase):
    def test_cells_in_bb(self):
        nwcorner = Point(150, 0)
        secorner = Point(-150, 0)
        cells_per_degree = 0.1
        for cell in RMGCell.cells_in_bb(nwcorner, secorner, startkey=None, cells_per_degree=cells_per_degree, a=a, inverse_flattening=inverse_flattening):
            print cell
        
        startkey = '34-9'
        print 'KEY: %s POLYGON: %s' % (startkey, RMGCell.polygon(startkey, cells_per_degree, digits, a, inverse_flattening))
        startkey = '35-9'
        print 'KEY: %s POLYGON: %s' % (startkey, RMGCell.polygon(startkey, cells_per_degree, digits, a, inverse_flattening))
        startkey = '0-9'
        print 'KEY: %s POLYGON: %s' % (startkey, RMGCell.polygon(startkey, cells_per_degree, digits, a, inverse_flattening))
        
        startkey = RMGCell.key(nwcorner.get_lng(), nwcorner.get_lat(), cells_per_degree, a, inverse_flattening)
        startkey = '34-9'
        print 'startkey = %s' % startkey
        cells=[]
        for cell in RMGCell.cells_in_bb(nwcorner, secorner, startkey=startkey, cells_per_degree=cells_per_degree):
            cells.append(cell)
            if len(cells) == 5:
                print cells
                cells=[]
                startkey = cell
        print cells
    
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

    def test_lat2y(self):
        cells_per_degree = 1
        lat = 90
        y_index = RMGCell.lat2y(lat, cells_per_degree)
        self.assertEqual(y_index, 0)

        lat = 0
        y_index = RMGCell.lat2y(lat, cells_per_degree)
        self.assertEqual(y_index, 90)

        lat = -90
        y_index = RMGCell.lat2y(lat, cells_per_degree)
        self.assertEqual(y_index, 180)

    def test_key0_0(self):
        lng = -180
        lat = 90
        cells_per_degree = 1
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,'0-0')

    def test_key0_1(self):
        lng = -180
        cells_per_degree = 1
        lat = 90 - 1.0/cells_per_degree
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,'0-1')

    def test_key1_1(self):
        lng = 180
        cells_per_degree = 1
        lat = 90 - 1.0/cells_per_degree
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        self.assertEqual(key,'9-1')

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
        # WorldClim Tile 11 cell count
        nwcorner = Point(-150, 60)
        secorner = Point(-120, 30)
        cells_per_degree = 120
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 9058790)
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
        #Cell count within Worldclim bounds at 120 cells per degree
        nwcorner = Point(-180, 90)
        secorner = Point(180, -60)
        cells_per_degree = 120
        tile = RMGTile(nwcorner, secorner, cells_per_degree)
        cellcount = tile.cellcount()
        self.assertEqual(cellcount, 552731769)
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
#        nwcorner = Point(30,0)
#        secorner = Point(60,-30)
        nwcorner = Point(42.883,0)
        secorner = Point(42.892,-0.009)
        cells_per_degree = 120
        digits = 7
        tile = RMGTile(nwcorner, secorner, cells_per_degree, digits)
        print tile.kml()

    def test_cell(self):
        cells_per_degree = 120
        lat = -5.45
        lng = 39.12
        key = RMGCell.key(lng, lat, cells_per_degree)
        print 'lat: '+str(lat) +' lng: '+str(lng) + ' key: '+key

        key='26567-10800'
        polygon = RMGCell.polygon(key, cells_per_degree)
        center = RMGCell.center(key)
        print 'key: %s center: %s polygon: %s' % (key, center, polygon)

        key='18-11'
        polygon = RMGCell.polygon(key, cells_per_degree)
        print polygon

        for corner in RMGCell.polygon(key, cells_per_degree):
            lng = float(corner[0])
            lat = float(corner[1])
            key = RMGCell.key(lng, lat, cells_per_degree)
            polygon = RMGCell.polygon(key, cells_per_degree)
            print key+' '+str(polygon)
            
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()

