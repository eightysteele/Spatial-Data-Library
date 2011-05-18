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
This module uses GDAL 1.8 ogr2ogr to clip multiple shapefiles against a single
shapefile.
"""

import logging
from optparse import OptionParser
import os
import shlex
import subprocess

def clip(indir, clipsrc, outdir):
    """Clips a directory of shapefiles with another shapefile file using ogr2ogr.

    Args:
        indir - Directory containing multiple shapefiles
        clipsrc - Source shapefile to clip from
        outdir - Directory to write resulting shapefiles
    """
    indir = os.path.abspath(indir)
    clipsrc = os.path.abspath(clipsrc)
    outdir = os.path.abspath(outdir)
    os.chdir(indir)
    shapefiles = [x for x in os.listdir('.') if x.endswith('.shp')]
    for v in shapefiles:
        clipped = os.path.join(outdir, v.replace('.shp', '-clipped.shp'))
        commandline = '/usr/local/bin/ogr2ogr -clipsrc %s %s %s' % (clipsrc, clipped, v)
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
                      help="The output directory for clipped shapefiles",
                      default=None)
    parser.add_option("-c", "--clipsrc", dest="clipsrc",
                      help="The clip src shapefile",
                      default=None)

    (options, args) = parser.parse_args()
    indir = options.indir
    clipsrc = options.clipsrc
    outdir = options.outdir
    
    clip(indir, clipsrc, outdir)
