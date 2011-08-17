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

'''
This module includes classes and a command line interface for bulkloading 
WorldClim environment variables to CouchDB using the Rectangular Mesh Grid (RMG).
'''

import csv
import sys
import shlex
import subprocess
import glob
import logging
from optparse import OptionParser
import os

VARDICT = {'tmean':'t', 'tmin':'m', 'tmax':'x', 'alt':'a', 'bio':'b', 'prec':'p'}
VARSET = ['tmean', 'tmin', 'tmax', 'alt', 'prec', 'bio']

def getminmax(vardir, varlimits):
    os.chdir(vardir)
    for f in glob.glob("*.hdr"):
        var = f.split('_')[0]
        file = open(f, 'r')
        for line in file:
            if line.find(' ')>0:
                if line.split()[0]=='MinValue':
                    min = int(line.split()[1])
                    if varlimits.get(var)==None:
                        varlimits[var]={'min':99999999, 'max':-99999999, 'file': '-1'}
                    if varlimits[var]['min']>min:
                        varlimits[var]['min']=min
                        varlimits[var]['file']=f
                if line.split()[0]=='MaxValue':
                    max = int(line.split()[1])
                    if varlimits.get(var)==None:
                        varlimits[var]={'min':99999999, 'max':-99999999, 'file': '-1'}
                    if varlimits[var]['max']<max:
                        varlimits[var]['max']=max
                        varlimits[var]['file']=f
        file.close()

def getworldclimtile(key, vardir):
    ''' Downloads the files for all of the variables for the Tile given by the Tile key
        in the command line, unzips them, and removes the zipfile.
        
        Arguments:
            key - the tile identifier (e.g., '37' for Worldclim tile 37)
            vardir - the path to the directory in which to store the Worldclim files
    '''
    if not os.path.exists(vardir):
        os.mkdir(vardir)
        if not os.path.exists(vardir):
            logging.info('Unable to create variable directory %s.' % vardir)
            return False

    for var in VARSET:
#        varfile = '%s_10m_bil.zip' % (var)
        varfile = '%s_%s.zip' % (var, key)
        varpath = os.path.join(vardir, varfile)
        ''' Get the file from Worldclim site.'''
#        command = 'wget -P %s http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/%s' % (vardir, varfile)
        command = 'wget -P %s http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/%s' % (vardir, varfile)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        ''' Unzip the variable zip file.'''
        command = 'unzip %s -d %s' % (varpath, vardir)
        logging.info(command)
        args = shlex.split(command)
        subprocess.call(args)
        ''' Remove the zip file.'''
        os.remove(varpath)
    return True

def _getoptions():
    ''' Parses command line options and returns them.'''
    parser = OptionParser()
    parser.add_option("-d", 
                      "--vardir", 
                      dest="vardir",
                      help="The directory of variable files.",
                      default=None)
    return parser.parse_args()[0]

def main():
    options = _getoptions()
    varlimits = {}

#    getworldclimtile('', options.vardir)
#    getminmax(options.vardir, varlimits)
    
    for i in range(0,5):
        for j in range (0,12):
            tile = '%s%s' %(i,j)
            savetodir = os.path.join(options.vardir,tile)
            getworldclimtile(tile, vardir)
            getminmax(options.vardir, varlimits)
        print 'Row %s Worldclim variable limits: %s' % varlimits

if __name__ == '__main__':
    main()

'''
[
  {
    "bio1": {
      "key": "bio1",
      "name": "Annual Mean Temperature",
      "minval": -269,
      "maxval": 314,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio2": {
      "key": "bio2",
      "name": "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
      "minval": 9,
      "maxval": 211,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio3": {
      "key": "bio3",
      "name": "Isothermality (BIO2/BIO7) (* 100)",
      "minval": 8,
      "maxval": 95,
      "unit": "%",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio4": {
      "key": "bio4",
      "name": "Temperature Seasonality (standard deviation *100)",
      "minval": 72,
      "maxval": 22673,
      "unit": "deg C * 100",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio5": {
      "key": "bio5",
      "name": "Max Temperature of Warmest Month",
      "minval": -59,
      "maxval": 489,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio6": {
      "key": "bio6",
      "name": "Min Temperature of Coldest Month",
      "minval": -547,
      "maxval": 258,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio7": {
      "key": "bio7",
      "name": "Temperature Annual Range (BIO5-BIO6)",
      "minval": 53,
      "maxval": 725,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio8": {
      "key": "bio8",
      "name": "Mean Temperature of Wettest Quarter",
      "minval": -251,
      "maxval": 375,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio9": {
      "key": "bio9",
      "name": "Mean Temperature of Driest Quarter",
      "minval": -450,
      "maxval": 364,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio10": {
      "key": "bio10",
      "name": "Mean Temperature of Warmest Quarter",
      "minval": -97,
      "maxval": 380,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio11": {
      "key": "bio11",
      "name": "Mean Temperature of Coldest Quarter",
      "minval": -488,
      "maxval": 289,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio12": {
      "key": "bio12",
      "name": "Annual Precipitation",
      "minval": 0,
      "maxval": 9916,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio13": {
      "key": "bio13",
      "name": "Precipitation of Wettest Month",
      "minval": 0,
      "maxval": 2088,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio14": {
      "key": "bio14",
      "name": "Precipitation of Driest Month",
      "minval": 0,
      "maxval": 652,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio15": {
      "key": "bio15",
      "name": "Precipitation Seasonality (Coefficient of Variation)",
      "minval": 0,
      "maxval": 261,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio16": {
      "key": "bio16",
      "name": "Precipitation of Wettest Quarter",
      "minval": 0,
      "maxval": 5043,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio17": {
      "key": "bio17",
      "name": "Precipitation of Driest Quarter",
      "minval": 0,
      "maxval": 2159,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio18": {
      "key": "bio18",
      "name": "Precipitation of Warmest Quarter",
      "minval": 0,
      "maxval": 4001,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "bio19": {
      "key": "bio19",
      "name": "Precipitation of Coldest Quarter",
      "minval": 0,
      "maxval": 3985,
      "unit": "mm",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "alt": {
      "key": "alt",
      "name": "Altitude (m)",
      "minval": -431,
      "maxval": 8233,
      "unit": "m",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin1": {
      "key": "tmin1",
      "name": "Minimum Temperature, January",
      "minval": -547,
      "maxval": 266,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin2": {
      "key": "tmin2",
      "name": "Minimum Temperature, February",
      "minval": -525,
      "maxval": 273,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin3": {
      "key": "tmin3",
      "name": "Minimum Temperature, March",
      "minval": -468,
      "maxval": 277,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin4": {
      "key": "tmin4",
      "name": "Minimum Temperature, April",
      "minval": -379,
      "maxval": 283,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin5": {
      "key": "tmin5",
      "name": "Minimum Temperature, May",
      "minval": -225,
      "maxval": 295,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin6": {
      "key": "tmin6",
      "name": "Minimum Temperature, June",
      "minval": -170,
      "maxval": 312,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin7": {
      "key": "tmin7",
      "name": "Minimum Temperature, July",
      "minval": -171,
      "maxval": 311,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin8": {
      "key": "tmin8",
      "name": "Minimum Temperature, August",
      "minval": -178,
      "maxval": 312,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin9": {
      "key": "tmin9",
      "name": "Minimum Temperature, September",
      "minval": -192,
      "maxval": 300,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin10": {
      "key": "tmin10",
      "name": "Minimum Temperature, October",
      "minval": -302,
      "maxval": 268,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin11": {
      "key": "tmin11",
      "name": "Minimum Temperature, November",
      "minval": -449,
      "maxval": 267,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin12": {
      "key": "tmin12",
      "name": "Minimum Temperature, December",
      "minval": -522,
      "maxval": 268,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmean1": {
      "key": "tmean1",
      "name": "Mean Temperature, January",
      "minval": -513,
      "maxval": 338,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmean2": {
      "key": "tmean2",
      "name": "Mean Temperature, February",
      "minval": -473,
      "maxval": 333,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin3": {
      "key": "tmin3",
      "name": "Minimum Temperature, March",
      "minval": -468,
      "maxval": 277,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin4": {
      "key": "tmin4",
      "name": "Minimum Temperature, April",
      "minval": -379,
      "maxval": 283,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin5": {
      "key": "tmin5",
      "name": "Minimum Temperature, May",
      "minval": -225,
      "maxval": 295,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin6": {
      "key": "tmin6",
      "name": "Minimum Temperature, June",
      "minval": -170,
      "maxval": 312,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin7": {
      "key": "tmin7",
      "name": "Minimum Temperature, July",
      "minval": -171,
      "maxval": 311,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin8": {
      "key": "tmin8",
      "name": "Minimum Temperature, August",
      "minval": -178,
      "maxval": 312,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin9": {
      "key": "tmin9",
      "name": "Minimum Temperature, September",
      "minval": -192,
      "maxval": 300,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin10": {
      "key": "tmin10",
      "name": "Minimum Temperature, October",
      "minval": -302,
      "maxval": 268,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin11": {
      "key": "tmin11",
      "name": "Minimum Temperature, November",
      "minval": -449,
      "maxval": 267,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
    "tmin12": {
      "key": "tmin12",
      "name": "Minimum Temperature, December",
      "minval": -522,
      "maxval": 268,
      "unit": "deg C * 10",
      "database": "WorldClim",
      "version": "1.4",
      "release": 3,
      "created": "2006-01-04",
      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    },
  }
]
'''
