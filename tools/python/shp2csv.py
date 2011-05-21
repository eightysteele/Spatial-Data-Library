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

"""
This module uses Starspan to calculate statistics for the intersection of 
shapefile vectors with worldclim raster layers.
"""

import logging
from optparse import OptionParser
import os
import shlex
import subprocess



def intersect(indir, raster, outdir):
    """Intersects a directory of shapefiles with a raster file using Starspan.

    Args:
        indir - Directory containing multiple shapefiles
        raster - A raster data file
        outdir - Directory to write resulting CSV files
    """
    indir = os.path.abspath(indir)
    raster = os.path.abspath(raster)
    outdir = os.path.abspath(outdir)

    os.chdir(indir)
    rasterfiles = [x for x in os.listdir(raster) if x.endswith('.bil')]
    vectorfiles = [x for x in os.listdir('.') if x.endswith('.shp')]
    for v in vectorfiles:
        csvfile = os.path.join(outdir, v.replace('.shp', '-starspan-stats-ave.csv'))
        commandline = 'starspan --vector %s --raster %s --csv %s' % (v, raster, csvfile)
        logging.info(commandline)
        args = shlex.split(commandline)
        subprocess.call(args)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)    
    
    # Parses command line parameters:
    parser = OptionParser()
    parser.add_option("-i", "--input-dir", dest="indir",
                      help="The input directory of shapefiles",
                      default=None)
    parser.add_option("-o", "--output-dir", dest="outdir",
                      help="The output directory for CSV files",
                      default=None)
    parser.add_option("-r", "--raster-file", dest="rasterfile",
                      help="The raster file",
                      default=None)

    (options, args) = parser.parse_args()
    indir = options.indir
    rasterfile = options.rasterfile
    outdir = options.outdir
    
    intersect(indir, rasterfile, outdir)
