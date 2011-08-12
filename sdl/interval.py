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

import math

def get_n(min,max):
    return int(math.ceil(math.log(math.ceil(max) - math.floor(min),2))+1)

def width(n):
    return int(math.pow(2, n))

def start(val,min,max,n):
    a = (val-min)
    b = (val-min)/(math.pow(2,n))
    c = int(math.floor((val-min)/(math.pow(2,n))))
    d = c + min
    while d <= val:
        d+=int(math.pow(2,n))
    if d==val:
        return d
    return d-int(math.pow(2,n))
 
def get_indexes(val,min,max,res):
    newmin = int(math.floor(min/res))
    newmax = int(math.ceil(max/res))
    category_count = get_n(newmin,newmax)
    intervals = []
    for n in range(category_count):
        index = start(val,newmin,newmax,n)
        print 'start %s end %s' % (index, index + int(pow(2,n)))
        intervals.append(index)
    return intervals

def get_index_intervals(val, min, max, res=1):
    intervals = get_indexes(val, min, max, res)
    indexes = dict()
    a = 0
    for i in intervals:
        index = 'i%s' % str(a) 
        indexes[index]=i
        a += 1
    return indexes

def main():
    indexes = get_index_intervals(320,-454,8550,1)
    print indexes

if __name__ == "__main__":
    main()

    
