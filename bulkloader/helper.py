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
import csv
import logging
from optparse import OptionParser
import os
import simplejson

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
        return transform.create_deep_key(
            ('Cell', 'cellkey'),
            ('CellIndex', 'cellkey'))(value, bulkload_state)
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

cells = {}
cells_loaded = {}

def write_cellkey_csv(filename):
    dr = csv.DictReader(open(filename, 'r'))
    cells = {}

    for row in dr:
        cellkey = row['CellKey']
        if not cells.has_key(cellkey):
            cells[cellkey] = {}
        varname = row['RID'].split('_')[0]
        varval = row['avg_Band1'].split('.')[0]
        cells[cellkey][varname] = varval

    name, ext = os.path.splitext(filename)
    dw = csv.DictWriter(open('%s.keys.csv' % name, 'w'), ['cellkey', 'varvals'])
    dw.writeheader()
    for k,v in cells.iteritems():
        dw.writerow(dict(cellkey=k,varvals=simplejson.dumps(v)))
    

def create_index(i):
    """
    bio1 in [-269,314]
    bio12 in [0,9916]
    """
    def wrapper(varvals, bulkload_state):
        varvals = simplejson.loads(varvals)
        varname, index = i.split('-')
        iname = 'i%s' % index        
        if varname == 'a':
            var_min = -454
            var_max = 8550
            varval = int(varvals['alt'])
        elif varname == 'b1':
            var_min = -269
            var_max = 314
            varval = int(varvals['bio1'])
        elif varname == 'b12':
            var_min = 0
            var_max = 9916
            varval = int(varvals['bio12'])
        else:
            logging.info('Unsupported range variable %s' % varname)
            return None
        intervals = interval.get_index_intervals(varval, var_min, var_max)
        for index,value in intervals.iteritems():
            if index == iname:
                return value
        logging.info('No %s range value for variable %s' % (iname, varname))
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


def _getoptions():
    ''' Parses command line options and returns them.'''
    parser = OptionParser()
    parser.add_option('--filename', 
                      type='string', 
                      dest='filename',
                      metavar='FILE', 
                      help='CSV file.')    
    return parser.parse_args()[0]

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    options = _getoptions()
    write_cellkey_csv(options.filename)
