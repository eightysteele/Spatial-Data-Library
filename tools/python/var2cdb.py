#!/usr/bin/env python

import csv
import hashlib
import logging
from optparse import OptionParser
import os
import random
import re
import simplejson
from uuid import uuid4

# CouchDB imports
import couchdb
from couchdb import Server

DB_FILE = 'load.sqlite3.db'
DB_CACHE_TABLE = 'cache'
DB_TMP_TABLE = 'tmp'
VAR_MIN = -50
VAR_MAX = 293
VAR_NODATA = -9999

def scaleval(val, varmin, varmax):
    try:        
        return varmin + float(val) * (varmax - varmin) / 255 
    except:
        return None

def within_list(val, within):
    x = int(val)
    return range(x - within, x + within + 1)

def mock_vars(val):
    vars = {'bio1':val}
    for x in range(2, 25):
        vars['bio%s' % x] = random.uniform(-200, 500)
    return vars

# FID,CellKey,RID,col,row,x,y,Band1
# 0,27178-13595,bio1_37.bil,1979,2795,46.483,-23.292,203

def load(csvdir, couchdb_url):
    csvdir = os.path.abspath(csvdir)
    server = Server(couchdb_url)
    sdl = server['sdl']    
    csvfiles = [x for x in os.listdir(csvdir) if x.endswith('.csv')]
    for csvfile in csvfiles:
        dr = csv.DictReader(open(os.path.join(csvdir, csvfile), 'r'))
        cells = {}
        for row in dr:
            cellkey = row.get('CellKey')
            if not cells.has_key(cellkey):
                cells[cellkey] = {
                    '_id': cellkey, 
                    'coords': [[random.uniform(x, 90),random.uniform(-180, x)] for x in range(-2,3)], # TODO
                    'vars': {}
                    }            
            varname = row.get('RID').split('_')[0]
            varval = row.get('Band1')
            cells.get(cellkey).get('vars')[varname] = varval
        logging.info('Bulkloading %s documents' % len(cells))
        sdl.update(cells.values())
            
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)    
    
    # Parses command line parameters:
    parser = OptionParser()
    parser.add_option("-d", "--csvdir", dest="csvdir",
                      help="The directory of CSV files",
                      default=None)
    parser.add_option("-u", "--url", dest="url",
                      help="The CouchDB URL",
                      default=None)

    (options, args) = parser.parse_args()
    csvdir = options.csvdir
    url = options.url
    
    # Writes the CSV file:
    load(csvdir, url)
