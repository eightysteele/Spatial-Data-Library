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
from collections import defaultdict
 
def get_indexes(val,min,max,res):
    newval = int(math.floor(val/res))
    newmin = int(math.floor(min/res))
    newmax = int(math.ceil(max/res))
    category_count = int(math.ceil(math.log(math.ceil(newmax) - math.floor(newmin),2)))
    intervals = []
    for n in range(category_count):
        index = newval-(newval-newmin)%(pow(2,n))
#        print 'n %s start %s end %s' % (n, index, index + int(pow(2,n)))
        intervals.append(index)
    return intervals

def get_index_intervals(val,min,max,res=1):
    ''' Returns a dictionary of interval:values where the values are the starting integers of a range containing val.

        Arguments:
            val - the value for which to create the interval dictionary
            min - the minimum value of the variable range, starting value for the intervals for a variable (e.g., -454m for altitude) 
            max - the maximum value of the variable range (e.g., 8850m for altitude)
            res - the resolution of the desired intervals (e.g., 1 means every meter, 10 means every 10 meters)
            
        Example:
            Altitudes between 315 and 330 (greater than or equal to 315, less than 330) using resolution=1.
            i0=315 [315,316} 1
            i1=316 [316,318} 2
            i2=318 [318,322} 4
            i3=322 [322,330} 8          

            Altitudes between 315 and 330 (greater than or equal to 315, less than 330) using resolution=1.
            i0=32 [320,330} 1
            i1=32 [320,340} 2
            i2=30 [300,340} 4
            i3=26 [260,340} 8
            i4=18 [180,340} 16
            i5=18 [180,500} 32
            i6=18 [180,820} 64
            i7=-46 [-460,820} 128
            i8=-46 [-460,2100} 256
            i9=-46 [-460,4660} 512
    '''
    intervals = get_indexes(val,min,max,res)
    indexes = dict()
    j = 0
    for i in intervals:
        index = 'i%s' % str(j) 
        indexes[index]=i
        j += 1
    return indexes

def get_query_intervals(min, max, gte, lt):
    ''' Returns a dictionary of interval:value-lists that need to be queried to get cells having a variable value within a given range.

        Arguments:
            min - the minimum value of the variable range, starting value for the intervals for a variable (e.g., -454m for altitude) 
            max - the maximum value of the variable range (e.g., 8850m for altitude)
            gte - the bottom end of the range of values to include in the results
            lt - an integer one greater than the top end of the range of values include in results
            
        Example:
            Altitudes between 315 and 330 (greater than or equal to 315, less than 329).
            i0=315,328 [315,316} and [328,329} 1
            i1=316,326 [316,318} and [326,328} 2
            i2=318,322 [318,322} and [322,326} 4
    '''
    ''' Clamp gte and lt to min and max values of the variable.'''
    if gte<min:
        gte=min
    if lt>max:
        lt=max
        
    indexes = defaultdict(list)
    ''' Power of 2 for the difference between max and min.'''
    n = int(math.floor(math.log( (lt-gte), 2)))
    ''' m is the magic starting value between min and max.''' 
    m = min+int(math.pow(2,n)*int(math.floor((lt-min)/math.pow(2,n))))
    ''' top is the difference between m and the top of the range.'''
    top = lt-m
    ''' bottom is the difference between m and the bottom of the range.'''
    bottom = m-gte
    diff = bottom
    ''' Iterate through the differences from m to the bottom of the range.'''
    while diff > 0:
        n = int(math.floor(math.log( diff, 2)))
        m=m-int(math.pow(2,n))
        ''' Add the value to the list of values for the interval.'''
        indexes['i%s' % n].append(m)
        diff = m-gte
    m=lt-top
    diff = top
    ''' Iterate through the differences from m to the top of the range.'''
    while diff > 0:
        n = int(math.floor(math.log( diff, 2)))
        ''' Add the value to the list of values for the interval.'''
        indexes['i%s' % n].append(m)
        m=m+int(math.pow(2,n))
        diff = lt-m
    return indexes    
    
def main():
    indexes = get_index_intervals(14,-20,40,1)
    print indexes

    indexes = get_index_intervals(320,-454,8550,1)
    print indexes

    indexes = get_index_intervals(320,-454,8550,10)
    print indexes

    intervals = get_query_intervals(-454, 8850, 315, 329)
    print intervals

    intervals = get_query_intervals(-20, 40, 10, 20)
    print intervals
    
    intervals = get_query_intervals(-20, 40, -21, 41)
    print intervals
    
if __name__ == "__main__":
    main()
