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

# Fix sys.path
from setup_env import fix_sys_path
fix_sys_path()

# Python imports
import logging

# SDL imports
from sdl import interval

# Goole App Engine imports
from django.utils import simplejson
from google.appengine.ext.bulkload import transform

# Datastore Plus imports
from ndb import query, model

def create_key():
    def wrapper(value, bulkload_state):
        '''Returns a CellIndex key with Cell as it's parent.'''
        d = bulkload_state.current_dictionary
        d['varname'] = d['RID'].split('_')[0]
        return transform.create_deep_key(
            ('Cell', 'CellKey'),
            ('CellIndex', 'varname'))(value, bulkload_state)
    return wrapper

def create_cell_key():
    def wrapper(value, bulkload_state):
        '''Returns a CellIndex key with Cell as it's parent.'''
        key_name = bulkload_state.current_dictionary['CellKey']
        return key_name
    return wrapper

def get_varname():
    def wrapper(value, bulkload_state):
        d = bulkload_state.current_dictionary
        return d['RID'].split('_')[0]
    return wrapper

def get_varval():
    def wrapper(value, bulkload_state):
        d = bulkload_state.current_dictionary
        val = d['avg_Band1'].split('.')[0]
        return int(val)
    return wrapper

def within_list(within):
    def wrapper(scaledval, bulkload_state):
        '''Returns a list of values before and after the scaledvalue.'''
        val = scaledval.split('.')[0]
        x = int(val)
        results = [str(x) for x in range(x - within, x + within + 1)]
        logging.info(results)
        return results
    return wrapper

def get_list(within, val):
    x = int(val.split('.')[0])
    return range(x - within, x + within + 1)

def create_index(i):
    """
    bio1 in [-269,314]
    bio12 in [0,9916]
    """
    def wrapper(varval, bulkload_state):
        iname = 'i%s' % i
        varval = int(varval.split('.')[0])
        varname = bulkload_state.current_dictionary['RID'].split('_')[0].lower()        
        if varname == 'alt':
            var_min = -454
            var_max = 8550
        elif varname == 'bio1':
            var_min = -269
            var_max = 314
        elif varname == 'bio12':
            var_min = 0
            var_max = 9916
        else:
            return None
        intervals = interval.get_index_intervals(varval, var_min, var_max)
        for index,value in intervals.iteritems():
            if index == iname:
                return value
        return None
    return wrapper
        
def add_dynamic_properties(input_dict, instance, bulkload_state_copy):    
    """Adds dynamic properties from the CSV input_dict to the entity instance."""
    val = int(input_dict['avg_Band1'].split('.')[0])
    var_min = -454
    var_max = 8550
    for k,v in interval.get_index_intervals(val, var_min, var_max).iteritems():
        instance[k] = v
    return instance


