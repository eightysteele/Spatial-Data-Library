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
from sdl.tmg import gettile

class GetTileTest(unittest.TestCase):

    def test_gettile(self):
        north = 0
        west = 30
        south = -30
        east = 60
        resolution = 10
        nwcorner = (west, north)
        secorner = (east, south)
        cells = gettile(nwcorner, secorner, resolution, 3)
        [logging.info(str(x)) for x in cells]
        
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
