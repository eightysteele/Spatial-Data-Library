#!/usr/bin/env python

import csv
import hashlib
import logging
import math
from optparse import OptionParser
import os
import random
import re
import simplejson
from uuid import uuid4

def getcellkey(globalrow, globalcol, tilerow, tilecol, res, dim):
    """Gets a global cell key."""
    x = int(math.floor((globalcol * (dim / res)) + tilecol))
    y = int(math.floor((globalrow * (dim /res)) + tilerow))
    return '%s-%s' % (x, y)
    
def splitcsv(csvfile):
    """Writes a CSV per cell key containing rows for each variable value."""
    data = {}
    reader = csv.DictReader(
        open(os.path.abspath(csvfile), 'r'))
    reader.next() # Skip header
    for line in reader:
        cellkey = getcellkey(3, 7, int(line.get('col')), int(line.get('row')), .008333, 30)
        line['CellKey'] = cellkey
        line.pop('FID')
        line.pop('x')
        line.pop('y')
        if not data.has_key(cellkey):
            data[cellkey] = []
        rows = data.get(cellkey)
        if len(rows) <= 68:
            data.get(cellkey).append(line)

    for cellkey in data.keys():
        writer = csv.DictWriter(
            open('%s.csv' % os.path.abspath(cellkey), 'w+'), 
            ['CellKey', 'RID', 'row', 'col', 'Band1'],
            extrasaction='ignore')   
        rows = data.get(cellkey)
        logging.info('%s rows for %s.csv' % (len(rows), cellkey))
        writer.writerows(rows)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)    
    
    # Parses command line parameters:
    parser = OptionParser()
    parser.add_option("-c", "--csvfile", dest="csvfile",
                      help="Input CSV file",
                      default=None)

    (options, args) = parser.parse_args()
    csvfile = options.csvfile
    
    splitcsv(csvfile)

