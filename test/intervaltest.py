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

"""This module provides unit testing for Bounding Box class."""

import logging
import sys
import unittest

sys.path.insert(0, '../')

from sdl.interval import *

class IntervalTest(unittest.TestCase):
    def test_index_intervals(self):
#        indexes = get_query_intervals(0,11401,1219,1223)
#        print ""
#        print indexes
#        indexes = get_query_intervals(0,11401,1219,1222)
#        print indexes
        print""
        min = -431
        max = 8233
        gte = 1219
        
#        lt = 1220
#        indexes = get_query_intervals(min,max,gte,lt)
#        print 'min: %s max: %s gte: %s lt: %s %s' % (min, max, gte, lt, indexes)
#        lt = 1221
#        indexes = get_query_intervals(min,max,gte,lt)
#        print 'min: %s max: %s gte: %s lt: %s %s' % (min, max, gte, lt, indexes)
#        lt = 1222
#        indexes = get_query_intervals(min,max,gte,lt)
#        print 'min: %s max: %s gte: %s lt: %s %s' % (min, max, gte, lt, indexes)
        lt = 1223

#        intervals = get_indexes(1219,min,max,1)
#        print 'Intervals: %s' % intervals
#        intervals = get_indexes(1221,min,max,1)
#        print 'Intervals: %s' % intervals
        
#        indexes = get_query_intervals2(min,max,gte,lt,1)
#        print 'min: %s max: %s gte: %s lt: %s %s' % (min, max, gte, lt, indexes)
#        sys.exit(1)

        indexes = get_index_intervals(1,0,1,1)
        self.assertEqual(indexes['i0'], 1)
        self.assertEqual(indexes.has_key('i1'), False)

        indexes = get_index_intervals(0,0,1,1)
        self.assertEqual(indexes['i0'], 0)
        self.assertEqual(indexes.has_key('i1'), False)

        indexes = get_index_intervals(-1,0,1,1)
        self.assertEqual(indexes, None)

        indexes = get_index_intervals(2,0,1,1)
        self.assertEqual(indexes, None)

        indexes = get_index_intervals(0,0,2,1)
        self.assertEqual(indexes['i0'], 0)
        self.assertEqual(indexes.has_key('i1'), False)

        indexes = get_index_intervals(1,0,2,1)
        self.assertEqual(indexes['i0'], 1)
        self.assertEqual(indexes.has_key('i1'), False)

        indexes = get_index_intervals(2,0,2,1)
        self.assertEqual(indexes['i0'], 2)
        self.assertEqual(indexes.has_key('i1'), False)

        indexes = get_index_intervals(0,0,3,1)
        self.assertEqual(indexes['i0'], 0)
        self.assertEqual(indexes['i1'], 0)
        self.assertEqual(indexes.has_key('i2'), False)

        indexes = get_index_intervals(1,0,3,1)
        self.assertEqual(indexes['i0'], 1)
        self.assertEqual(indexes['i1'], 0)
        self.assertEqual(indexes.has_key('i2'), False)

        indexes = get_index_intervals(2,0,3,1)
        self.assertEqual(indexes['i0'], 2)
        self.assertEqual(indexes['i1'], 2)
        self.assertEqual(indexes.has_key('i2'), False)

        indexes = get_index_intervals(3,0,3,1)
        self.assertEqual(indexes['i0'], 3)
        self.assertEqual(indexes['i1'], 2)
        self.assertEqual(indexes.has_key('i2'), False)

        intervals = get_query_intervals(-20, 40, -21, 41, 1)
        self.assertEqual(intervals['i0'], 40)
        self.assertEqual(intervals.has_key('i1'), False)
        self.assertEqual(intervals['i2'], 36)
        self.assertEqual(intervals['i3'], 28)
        self.assertEqual(intervals['i4'], 12)
        self.assertEqual(intervals['i5'], -20)
        self.assertEqual(intervals.has_key('i6'), False)

        indexes = get_index_intervals(1219,-431,8233,1)
        self.assertEqual(indexes['i0'], 1219)
        self.assertEqual(indexes['i1'], 1219)
        self.assertEqual(indexes['i2'], 1217)
        self.assertEqual(indexes['i3'], 1217)
        self.assertEqual(indexes['i4'], 1217)
        self.assertEqual(indexes['i5'], 1201)
        self.assertEqual(indexes['i6'], 1169)
        self.assertEqual(indexes['i7'], 1105)
        self.assertEqual(indexes['i8'], 1105)
        self.assertEqual(indexes['i9'], 1105)
        self.assertEqual(indexes['i10'], 593)
        self.assertEqual(indexes['i11'], -431)
        self.assertEqual(indexes['i12'], -431)
        self.assertEqual(indexes['i13'], -431)
        self.assertEqual(indexes.has_key('i14'), False)

        indexes = get_index_intervals(14,-20,40,1)
        self.assertEqual(indexes['i0'], 14)
        self.assertEqual(indexes['i1'], 14)
        self.assertEqual(indexes['i2'], 12)
        self.assertEqual(indexes['i3'], 12)
        self.assertEqual(indexes['i4'], 12)
        self.assertEqual(indexes['i5'], 12)
        self.assertEqual(indexes.has_key('i6'), False)
        
        indexes = get_index_intervals(17.3,-20,40,1)
        self.assertEqual(indexes['i0'], 17)
        self.assertEqual(indexes['i1'], 16)
        self.assertEqual(indexes['i2'], 16)
        self.assertEqual(indexes['i3'], 12)
        self.assertEqual(indexes['i4'], 12)
        self.assertEqual(indexes['i5'], 12)
        self.assertEqual(indexes.has_key('i6'), False)

        indexes = get_index_intervals(320,-454,8850,1)
        self.assertEqual(indexes['i0'], 320)
        self.assertEqual(indexes['i1'], 320)
        self.assertEqual(indexes['i2'], 318)
        self.assertEqual(indexes['i3'], 314)
        self.assertEqual(indexes['i4'], 314)
        self.assertEqual(indexes['i5'], 314)
        self.assertEqual(indexes['i6'], 314)
        self.assertEqual(indexes['i7'], 314)
        self.assertEqual(indexes['i8'], 314)
        self.assertEqual(indexes['i9'], 58)
        self.assertEqual(indexes['i10'], -454)
        self.assertEqual(indexes['i11'], -454)
        self.assertEqual(indexes['i12'], -454)
        self.assertEqual(indexes['i13'], -454)
        self.assertEqual(indexes.has_key('i14'), False)

        indexes = get_index_intervals(320,-454,8850,10)
        self.assertEqual(indexes['i0'], 32)
        self.assertEqual(indexes['i1'], 32)
        self.assertEqual(indexes['i2'], 30)
        self.assertEqual(indexes['i3'], 26)
        self.assertEqual(indexes['i4'], 18)
        self.assertEqual(indexes['i5'], 18)
        self.assertEqual(indexes['i6'], 18)
        self.assertEqual(indexes['i7'], -46)
        self.assertEqual(indexes['i8'], -46)
        self.assertEqual(indexes['i9'], -46)
        self.assertEqual(indexes.has_key('i10'), False)

#    def test_query_intervals(self):
#        intervals = get_query_intervals(-454, 8850, 315, 330, 1)
#        self.assertEqual(intervals.has_key('e'), False)
#        self.assertEqual(intervals['i0'], 315)
#        self.assertEqual(intervals['i1'], 316)
#        self.assertEqual(intervals['i2'], 318)
#        self.assertEqual(intervals['i3'], 322)
#        self.assertEqual(intervals.has_key('i4'), False)
#
#        intervals = get_query_intervals(-454, 8850, 315, 329, 1)
#        self.assertEqual(intervals['e'], 1)
#        self.assertEqual(intervals['i0'], 315)
#        self.assertEqual(intervals['i1'], 316)
#        self.assertEqual(intervals['i2'], 318)
#        self.assertEqual(intervals['i3'], 322)
#        self.assertEqual(intervals.has_key('i4'), False)
#
#        intervals = get_query_intervals(-454, 8850, 315, 328, 1)
#        self.assertEqual(intervals['e'], 2)
#        self.assertEqual(intervals['i0'], 315)
#        self.assertEqual(intervals['i1'], 316)
#        self.assertEqual(intervals['i2'], 318)
#        self.assertEqual(intervals['i3'], 322)
#        self.assertEqual(intervals.has_key('i4'), False)
#
#        intervals = get_query_intervals(-20, 40, 10, 20, 1)
#        self.assertEqual(intervals.has_key('e'), False)
#        self.assertEqual(intervals.has_key('i0'), False)
#        self.assertEqual(intervals['i1'], 10)
#        self.assertEqual(intervals.has_key('i2'), False)
#        self.assertEqual(intervals['i3'], 12)
#        self.assertEqual(intervals.has_key('i4'), False)
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
