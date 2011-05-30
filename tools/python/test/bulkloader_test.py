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

"""This module includes classes for unit testing the bulkloader module."""

from bulkloader import Tile, TileCell
import logging
import unittest

class TileTest(unittest.TestCase):
    def test_getcells(self):
        tile = Tile(3, 7)
        cells = []
        for cell in tile.getcells(30):
            cells.append(cell)
        self.assertEqual(len(cells), 1)
        logging.info(cell)

        tile = Tile(3, 7)
        cells = []
        for cell in tile.getcells(1):
            cells.append(cell)
        self.assertEqual(len(cells), 900)

class TileCellTest(unittest.TestCase):
    polygon = [[1,2],[3,4],[5,6],[7,8],[1,2]]
    cell = TileCell(1, 'key', polygon)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()


