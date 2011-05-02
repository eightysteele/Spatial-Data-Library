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

"""This module provides unit testing for TMG Cell class."""

# Python module imports:
import logging
import sys
import os
import unittest

# Updates sys.path to include geomancer source:
sys.path = [os.path.abspath(os.path.realpath('../'))] + sys.path

# TMG imports
from sdl.tmg import Cell

class CellTest(unittest.TestCase):

    def test_eq_hash_cmp(self):
        c1 = Cell(1, 2, 3)
        c2 = Cell(1, 2, 3)
        self.assertEqual(c1, c2)

        c3 = Cell(11, 22, 33)
        self.assertNotEqual(c1, c3)

        s = set()
        s.add(c1)
        s.add(c2)
        s.add(Cell(1, 2, 3))
        self.assertEqual(1, len(s))
        self.assertEqual(Cell(1, 2, 3), s.pop())

        d = {}
        d[c1] = None
        d[Cell(1, 2, 3)] = None
        self.assertEqual(1, len(d.keys()))

        l = [Cell(1, 2, 3), Cell(1, 1, 2), Cell(11, 22, 33)]
        l.sort()
        self.assertEqual(l[2], Cell(11, 22, 33))
        self.assertEqual(l[1], Cell(1, 2, 3))
        self.assertEqual(l[0], Cell(1, 1, 2))
        
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
