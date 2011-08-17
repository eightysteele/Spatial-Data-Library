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
import simplejson
import logging
from optparse import OptionParser
import os

knownvarlimits = {'bio19': {'max': 5162, 'file': 'bio19_28.hdr', 'min': 0}, 'tmax12': {'max': 417, 'file': 'tmax12_39.hdr', 'min': -484}, 'tmax11': {'max': 404, 'file': 'tmax11_39.hdr', 'min': -403}, 'tmax10': {'max': 402, 'file': 'tmax10_27.hdr', 'min': -286}, 'prec9': {'max': 1231, 'file': 'prec9_23.hdr', 'min': 0}, 'bio10': {'max': 383, 'file': 'bio10_25.hdr', 'min': -143}, 'bio11': {'max': 289, 'file': 'bio11_27.hdr', 'min': -521}, 'bio12': {'max': 11401, 'file': 'bio12_29.hdr', 'min': 0}, 'bio13': {'max': 2949, 'file': 'bio13_29.hdr', 'min': 0}, 'bio14': {'max': 752, 'file': 'bio14_23.hdr', 'min': 0}, 'bio4': {'max': 22721, 'file': 'bio4_21.hdr', 'min': 0}, 'bio16': {'max': 8019, 'file': 'bio16_29.hdr', 'min': 0}, 'bio17': {'max': 2495, 'file': 'bio17_23.hdr', 'min': 0}, 'bio18': {'max': 6090, 'file': 'bio18_29.hdr', 'min': 0}, 'prec12': {'max': 1041, 'file': 'prec12_23.hdr', 'min': 0}, 'alt': {'max': 8233, 'file': 'alt_28.hdr', 'min': -431}, 'prec10': {'max': 1238, 'file': 'prec10_23.hdr', 'min': 0}, 'prec11': {'max': 1109, 'file': 'prec11_23.hdr', 'min': 0}, 'tmean6': {'max': 386, 'file': 'tmean6_43.hdr', 'min': -195}, 'bio2': {'max': 214, 'file': 'bio2_21.hdr', 'min': 0}, 'tmin11': {'max': 269, 'file': 'tmin11_310.hdr', 'min': -483}, 'tmin10': {'max': 278, 'file': 'tmin10_25.hdr', 'min': -322}, 'bio3': {'max': 96, 'file': 'bio3_23.hdr', 'min': 0}, 'bio9': {'max': 366, 'file': 'bio9_17.hdr', 'min': -521}, 'tmin6': {'max': 325, 'file': 'tmin6_43.hdr', 'min': -253}, 'bio6': {'max': 258, 'file': 'bio6_30.hdr', 'min': -573}, 'prec1': {'max': 1003, 'file': 'prec1_33.hdr', 'min': 0}, 'tmean3': {'max': 337, 'file': 'tmean3_25.hdr', 'min': -465}, 'tmean2': {'max': 335, 'file': 'tmean2_39.hdr', 'min': -506}, 'tmean1': {'max': 340, 'file': 'tmean1_39.hdr', 'min': -536}, 'prec8': {'max': 2179, 'file': 'prec8_29.hdr', 'min': 0}, 'tmean7': {'max': 394, 'file': 'tmean7_43.hdr', 'min': -193}, 'tmin3': {'max': 280, 'file': 'tmin3_25.hdr', 'min': -495}, 'tmean5': {'max': 363, 'file': 'tmean5_28.hdr', 'min': -264}, 'tmean4': {'max': 345, 'file': 'tmean4_25.hdr', 'min': -381}, 'prec3': {'max': 827, 'file': 'prec3_311.hdr', 'min': 0}, 'prec2': {'max': 851, 'file': 'prec2_33.hdr', 'min': 0}, 'tmean9': {'max': 360, 'file': 'tmean9_25.hdr', 'min': -192}, 'tmean8': {'max': 384, 'file': 'tmean8_43.hdr', 'min': -190}, 'prec7': {'max': 2949, 'file': 'prec7_29.hdr', 'min': 0}, 'prec6': {'max': 2891, 'file': 'prec6_29.hdr', 'min': 0}, 'prec5': {'max': 1312, 'file': 'prec5_29.hdr', 'min': 0}, 'prec4': {'max': 924, 'file': 'prec4_23.hdr', 'min': 0}, 'tmean12': {'max': 333, 'file': 'tmean12_39.hdr', 'min': -519}, 'bio5': {'max': 490, 'file': 'bio5_25.hdr', 'min': -96}, 'tmean11': {'max': 330, 'file': 'tmean11_310.hdr', 'min': -443}, 'tmax8': {'max': 475, 'file': 'tmax8_43.hdr', 'min': -120}, 'tmean10': {'max': 330, 'file': 'tmean10_25.hdr', 'min': -304}, 'tmax5': {'max': 442, 'file': 'tmax5_25.hdr', 'min': -236}, 'tmax3': {'max': 423, 'file': 'tmax3_27.hdr', 'min': -449}, 'tmin12': {'max': 270, 'file': 'tmin12_310.hdr', 'min': -553}, 'tmax2': {'max': 417, 'file': 'tmax2_27.hdr', 'min': -454}, 'tmin5': {'max': 300, 'file': 'tmin5_27.hdr', 'min': -292}, 'tmin4': {'max': 286, 'file': 'tmin4_25.hdr', 'min': -398}, 'tmin7': {'max': 316, 'file': 'tmin7_43.hdr', 'min': -251}, 'bio1': {'max': 320, 'file': 'bio1_27.hdr', 'min': -290}, 'tmin1': {'max': 266, 'file': 'tmin1_35.hdr', 'min': -573}, 'bio7': {'max': 725, 'file': 'bio7_21.hdr', 'min': 0}, 'tmax9': {'max': 443, 'file': 'tmax9_25.hdr', 'min': -179}, 'tmin2': {'max': 273, 'file': 'tmin2_35.hdr', 'min': -558}, 'tmax7': {'max': 490, 'file': 'tmax7_43.hdr', 'min': -134}, 'tmax6': {'max': 469, 'file': 'tmax6_43.hdr', 'min': -136}, 'bio8': {'max': 378, 'file': 'bio8_25.hdr', 'min': -285}, 'tmax4': {'max': 430, 'file': 'tmax4_25.hdr', 'min': -364}, 'tmin9': {'max': 316, 'file': 'tmin9_43.hdr', 'min': -258}, 'tmin8': {'max': 313, 'file': 'tmin8_43.hdr', 'min': -260}, 'tmax1': {'max': 419, 'file': 'tmax1_39.hdr', 'min': -500}, 'bio15': {'max': 265, 'file': 'bio15_27.hdr', 'min': 0}}

VARDICT = {'tmean':'t', 'tmin':'m', 'tmax':'x', 'alt':'a', 'bio':'b', 'prec':'p'}
VARSET = ['tmean', 'tmin', 'tmax', 'alt', 'prec', 'bio']
MONTHS = {
          '1':'January', '2':'February', '3':'March','4':'April',
          '5':'May','6':'June','7':'July','8':'August',
          '9':'September','10':'October','11':'November','12':'December'
          }
VARCOMMONMETADATA = {
                      "database": "WorldClim",
                      "version": "1.4",
                      "release": 3,
                      "created": "2006-01-04",
                      "srs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    }
TEMPMETADATA = {
                'tmin':'Minimum Temperature', 'tmean':'Mean Temperature', 'tmax':'Maximum Temperature', 'prec':'Precipitation'
                }
TEMPUNITS = "deg C * 10"
PRECUNITS = "mm"
ALTUNITS = "m"
VARMETADATA = {
               "bio1": {
                        "name": "Annual Mean Temperature",
                        "unit": "deg C * 10"
               },
               "bio2": {
                      "name": "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                      "unit": "deg C * 10"
               },
               "bio3": {
                      "name": "Isothermality (BIO2/BIO7) (* 100)",
                      "unit": "%"
               },
               "bio4": {
                      "name": "Temperature Seasonality (standard deviation *100)",
                      "unit": "deg C * 100"
               },
               "bio5": {
                      "name": "Max Temperature of Warmest Month",
                      "unit": "deg C * 10"
               },
               "bio6": {
                      "name": "Min Temperature of Coldest Month",
                      "unit": "deg C * 10"
               },
               "bio7": {
                      "name": "Temperature Annual Range (BIO5-BIO6)",
                      "unit": "deg C * 10"
               },
               "bio8": {
                      "name": "Mean Temperature of Wettest Quarter",
                      "unit": "deg C * 10"
               },
               "bio9": {
                      "name": "Mean Temperature of Driest Quarter",
                      "unit": "deg C * 10"
               },
               "bio10": {
                      "name": "Mean Temperature of Warmest Quarter",
                      "unit": "deg C * 10"
               },
               "bio11": {
                      "name": "Mean Temperature of Coldest Quarter",
                      "unit": "deg C * 10"
               },
               "bio12": {
                      "name": "Annual Precipitation",
                      "unit": "mm"
               },
               "bio13": {
                      "name": "Precipitation of Wettest Month",
                      "unit": "mm"
               },
               "bio14": {
                      "name": "Precipitation of Driest Month",
                      "unit": "mm"
               },
               "bio15": {
                      "name": "Precipitation Seasonality (Coefficient of Variation)",
                      "unit": "mm"
               },
               "bio16": {
                      "name": "Precipitation of Wettest Quarter",
                      "unit": "mm"
               },
               "bio17": {
                      "name": "Precipitation of Driest Quarter",
                      "unit": "mm"
               },
               "bio18": {
                      "name": "Precipitation of Warmest Quarter",
                      "unit": "mm"
               },
               "bio19": {
                      "name": "Precipitation of Coldest Quarter",
                      "unit": "mm"
               },
               "alt": {
                      "name": "Altitude",
                      "unit": "m"
               }
            }

def getmetadata(varlimits):
    metadata = {}

    var = 'alt'
    metadata[var]={}
    metadata[var]['key']=var
    metadata[var]['name']='Altitude'
    metadata[var]['minval']=varlimits[var]['min']
    metadata[var]['maxval']=varlimits[var]['max']
    metadata[var]['unit']=ALTUNITS
    for f in VARCOMMONMETADATA.keys():
        metadata[var][f]=VARCOMMONMETADATA[f]
    
    for i in range(1,2):
        var = 'tmean%s'%(i)
        metadata[var]={}
        metadata[var]['key']=var
        metadata[var]['name']='%s, %s' % (TEMPMETADATA['tmean'],MONTHS[str(i)])
        metadata[var]['minval']=varlimits[var]['min']
        metadata[var]['maxval']=varlimits[var]['max']
        metadata[var]['unit']=TEMPUNITS
        for f in VARCOMMONMETADATA.keys():
            metadata[var][f]=VARCOMMONMETADATA[f]
            
        var = 'tmin%s'%(i)
        metadata[var]={}
        metadata[var]['key']=var
        metadata[var]['name']='%s, %s' % (TEMPMETADATA['tmean'],MONTHS[str(i)])
        metadata[var]['minval']=varlimits[var]['min']
        metadata[var]['maxval']=varlimits[var]['max']
        metadata[var]['unit']=TEMPUNITS
        for f in VARCOMMONMETADATA.keys():
            metadata[var][f]=VARCOMMONMETADATA[f]
            
        var = 'tmax%s'%(i)
        metadata[var]={}
        metadata[var]['key']=var
        metadata[var]['name']='%s, %s' % (TEMPMETADATA['tmean'],MONTHS[str(i)])
        metadata[var]['minval']=varlimits[var]['min']
        metadata[var]['maxval']=varlimits[var]['max']
        metadata[var]['unit']=TEMPUNITS
        for f in VARCOMMONMETADATA.keys():
            metadata[var][f]=VARCOMMONMETADATA[f]

        var = 'prec%s'%(i)
        metadata[var]={}
        metadata[var]['key']=var
        metadata[var]['name']='%s, %s' % (TEMPMETADATA['prec'],MONTHS[str(i)])
        metadata[var]['minval']=varlimits[var]['min']
        metadata[var]['maxval']=varlimits[var]['max']
        metadata[var]['unit']=PRECUNITS
        for f in VARCOMMONMETADATA.keys():
            metadata[var][f]=VARCOMMONMETADATA[f]
   
    for i in range(1,20):
        var = 'bio%s'%(i)
        metadata[var]={}
        metadata[var]['key']=var
        metadata[var]['name']=VARMETADATA[var]['name']
        metadata[var]['minval']=varlimits[var]['min']
        metadata[var]['maxval']=varlimits[var]['max']
        metadata[var]['unit']=VARMETADATA[var]['unit']
        for f in VARCOMMONMETADATA.keys():
            metadata[var][f]=VARCOMMONMETADATA[f]
    return metadata

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
    
    metadata = getmetadata(knownvarlimits)
    print simplejson.dumps(metadata)
    sys.exit(1)
    
    varlimits = {}

    for i in range(0,5):
        for j in range (0,12):
            tile = '%s%s' %(i,j)
            savetodir = os.path.join(options.vardir,tile)
#            getworldclimtile(tile, savetodir)
            getminmax(savetodir, varlimits)
    print 'Worldclim variable limits: %s' % (varlimits)

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
    }
  }
]
'''
'''
'tmean3': {
  'key': 'tmean3', 
  'name': 'Mean Temperature', 
  'minval': -465, 
  'maxval': 337, 
  'unit': 'deg C * 10',
  'database': 'WorldClim', 
  'version': '1.4', 
  'release': 3, 
  'created': '2006-01-04', 
  'srs': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' 
}, 

'''
