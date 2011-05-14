#!/usr/bin/env python

import csv
import hashlib
import logging
from optparse import OptionParser
import os
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
    return varmin + val * (varmax - varmin) / 255 

def within_list(val, within):
    x = int(val)
    return range(x - within, x + within + 1)

def load(csvdir, couchdb_url):
    csvdir = os.path.abspath(csvdir)
    server = Server(couchdb_url)
    sdl = server['sdl']    
    csvfiles = [x for x in os.listdir(csvdir) if x.endswith('.csv')]
    for csvfile in csvfiles:
        dr = csv.DictReader(open(os.path.join(csvdir, csvfile), 'r'))
        docs = []
        for row in dr:
            cell_value = scaleval(float(row['avg_Band1']), VAR_MIN, VAR_MAX)
            if cell_value > VAR_MAX or cell_value < VAR_MIN or cell_value is VAR_NODATA:
                continue
            docs.append({                
                    '_id': row['CellKey'],
                    'coords': [], # TODO
                    'vars': {row['RID'].split('_')[0]: float(row['avg_Band1'])}
                    })
        logging.info('Bulkloading %s documents' % len(docs))
        sdl.update(docs)

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
