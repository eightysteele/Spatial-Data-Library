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

"""This module contains some shit."""

# Standard Python imports:
from datetime import datetime
import logging
import os
import simplejson

# Google App Engine imports:
from google.appengine.api import mail, memcache, urlfetch
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app

# SDL imports:
from sdl import rmg

# Datastore Plus imports
from ndb import query, model

# CouchDb connection parameters:
COUCHDB_HOST = 'http://eighty.berkeley.edu'
COUCHDB_PORT = 5984
COUCHDB_DATABASE = 'worldclim-rmg'
COUCHDB_DESIGN = 'api'
COUCHDB_VIEW = 'cells'
COUCHDB_URL = '%s:%s/%s/_design/%s/_view/%s' % \
    (COUCHDB_HOST, 
     COUCHDB_PORT, 
     COUCHDB_DATABASE, 
     COUCHDB_DESIGN, 
     COUCHDB_VIEW)

# SI conversion functions for variable values
SI_CONVERSIONS = dict(
    alt=lambda x: int(x),
    bio1=lambda x: int(x)/10.0,
    bio2=lambda x: int(x)/10.0,
    bio3=lambda x: int(x)/100.0,
    bio4=lambda x: int(x)/100.0,
    bio5=lambda x: int(x)/10.0,
    bio6=lambda x: int(x)/10.0,
    bio7=lambda x: int(x)/10.0,
    bio8=lambda x: int(x)/10.0,
    bio9=lambda x: int(x)/10.0,
    bio10=lambda x: int(x)/10.0,
    bio11=lambda x: int(x)/10.0,
    bio12=lambda x: int(x),
    bio13=lambda x: int(x),
    bio14=lambda x: int(x),
    bio15=lambda x: int(x),
    bio16=lambda x: int(x),
    bio17=lambda x: int(x),
    bio18=lambda x: int(x),
    bio19=lambda x: int(x),
    prec1=lambda x: int(x),
    prec10=lambda x: int(x),
    prec11=lambda x: int(x),
    prec12=lambda x: int(x),
    prec2=lambda x: int(x),
    prec3=lambda x: int(x),
    prec4=lambda x: int(x),
    prec5=lambda x: int(x),
    prec6=lambda x: int(x),
    prec7=lambda x: int(x),
    prec8=lambda x: int(x),
    prec9=lambda x: int(x),
    tmax1=lambda x: int(x)/10.0,
    tmax10=lambda x: int(x)/10.0,
    tmax11=lambda x: int(x)/10.0,
    tmax12=lambda x: int(x)/10.0,
    tmax2=lambda x: int(x)/10.0,
    tmax3=lambda x: int(x)/10.0,
    tmax4=lambda x: int(x)/10.0,
    tmax5=lambda x: int(x)/10.0,
    tmax6=lambda x: int(x)/10.0,
    tmax7=lambda x: int(x)/10.0,
    tmax8=lambda x: int(x)/10.0,
    tmax9=lambda x: int(x)/10.0,
    tmean1=lambda x: int(x)/10.0,
    tmean10=lambda x: int(x)/10.0,
    tmean11=lambda x: int(x)/10.0,
    tmean12=lambda x: int(x)/10.0,
    tmean2=lambda x: int(x)/10.0,
    tmean3=lambda x: int(x)/10.0,
    tmean4=lambda x: int(x)/10.0,
    tmean5=lambda x: int(x)/10.0,
    tmean6=lambda x: int(x)/10.0,
    tmean7=lambda x: int(x)/10.0,
    tmean8=lambda x: int(x)/10.0,
    tmean9=lambda x: int(x)/10.0,
    tmin1=lambda x: int(x)/10.0,
    tmin10=lambda x: int(x)/10.0,
    tmin11=lambda x: int(x)/10.0,
    tmin12=lambda x: int(x)/10.0,
    tmin2=lambda x: int(x)/10.0,
    tmin3=lambda x: int(x)/10.0,
    tmin4=lambda x: int(x)/10.0,
    tmin5=lambda x: int(x)/10.0,
    tmin6=lambda x: int(x)/10.0,
    tmin7=lambda x: int(x)/10.0,
    tmin8=lambda x: int(x)/10.0,
    tmin9=lambda x: int(x)/10.0)

class Cell(model.Model):
    """Models a CouchDB cell document.

    key_name - The cell key (e.g., 1-2).
    """
    rev = model.StringProperty('r')
    coords = model.StringProperty('c')
    varvals = model.TextProperty('v')

    def __eq__(self, other):
        if isinstance(other, Cell):
            return self.key() == other.key()
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        return hash(self.key().name())

    def __cmp__(self, other):
        return self.key().__cmp__(other.key())

class CellIndex(model.Expando): # parent=Cell, key_name=varname    
    n = model.StringProperty('n', required=True)
    v = model.IntegerProperty('v', required=True)

    @classmethod
    def search(cls, varname, within, pivot, limit, offset):
        prop = 'within_%s' % within
        gql = "SELECT * FROM CellIndex WHERE n = '%s'" % varname
        gql = "%s AND %s = %d" % (gql, prop, int(pivot))
        logging.info(gql)
        qry = query.parse_gql(gql)[0]
        logging.info('QUERY='+str(qry))
        results = qry.fetch(limit, offset=offset, keys_only=True)
        cell_keys = [key.parent().id() for key in results]
        return set(cell_keys)

class CellValuesHandler(webapp.RequestHandler):
    """Handler for cell value requests."""

    @classmethod
    def fromds(cls, cell_keys):
        """Returns CouchDBCell entities from a datastore query on cell keys.

        Arguments:
            cell_keys - A set of cell key strings (e.g., 9-15).

        Returns:
            A dictionary of cell key to CouchDBCell.
        """

        # test
        keys = [model.Key('Cell', x) for x in cell_keys]
        entities = model.get_multi(keys)

        #entities = Cell.get_by_key_name(cell_keys)
        cells = {}
        for x in entities:
            if x:
                cells[x.key.id()] = x
        return cells

    @classmethod
    def fromcouchdb(cls, cell_keys):
        """Returns a list of CouchDBCell entities from a CouchDB query on cell keys.
        
        Arguments:
            cell_keys - A set of cell key strings (e.g., 9-15).

        Returns:
            A dictionary of cell key to CouchDBCell instance.
        """
        response = urlfetch.fetch(
            url=COUCHDB_URL,
            payload=simplejson.dumps({'keys': list(cell_keys)}),
            method=urlfetch.POST,
            headers={"Content-Type":"application/json"})
        if response.status_code != 200:
            return {}
        cells = {}        
        for row in simplejson.loads(response.content).get('rows'):            
            key = row.get('key')
            value = row.get('value')
            cells[key] = Cell(
                id=key,
                rev=value.get('rev'),
                coords=simplejson.dumps(value.get('coords')),
                varvals=simplejson.dumps(value.get('varvals')))
        return cells

    @classmethod
    def getcells(cls, cell_keys):
        """Gets CouchDBCell entities corresponding to a set of cell keys.
        
        Arguments:
            cell_keys - A set of cell key strings (e.g., 9-15).

        Returns:
            A dictionary of cell key to CouchDBCell.
        """
        cells = {}
        
        # Get cached cells
        cached = memcache.get_multi(cell_keys)        
        cells.update(cached)
        
        # Calculate any uncached cell keys
        cell_keys = cell_keys.difference(cached.keys())

        cachecells = False

        # Check datastore for any uncached cell keys
        if len(cell_keys) > 0:
            cachecells = True
            stored = cls.fromds(cell_keys)
            cells.update(stored)
            # Calculate any cell keys not in datastore
            cell_keys = cell_keys.difference(stored.keys())
            
        # Check CouchDB for any cells not in datastore
        if len(cell_keys) > 0:
            cachecells = True
            couched = cls.fromcouchdb(cell_keys)
            cells.update(couched)

            # Put cells from CouchDB into datastore
            #db.put(couched.values())
            model.put_multi(couched.values())
            
        # Cache cells
        if cachecells:
            memcache.set_multi(cells)                        

        return cells

    def cell_keys_from_coords(self, coords):
        """Returns a list of cell keys generated from a list of coords.

        Arguments:
            coords - list of lon,lat pairs
        """
        cell_keys = set()
        for coord in coords:
            lon, lat = coord.split(',')
            cell_key = rmg.RMGCell.key(
                float(lon), 
                float(lat), 
                rmg.CELLS_PER_DEGREE, 
                rmg.SEMI_MAJOR_AXIS, 
                rmg.INVERSE_FLATTENING)
            cell_keys.add(cell_key)
        return cell_keys

    def get(self):
        """Handles a cell API request by proxing to post."""
        return self.post()

    def post(self):
        """Handles a cell API request.

        URL parameters:
            xy - coordinates (lon,lat|lon,lat|...)
            k - cell keys (cellkey,cellkey,...)
            v - variable names (varname,varname,...)
            c - return cell coordinates if true
            si - return converted variable values to standard SI units if true
            bb - bounding box (north,west|south,east)
            bbo - bounding box offset cell key
            bbl - bounding box cell limit
            range, pivot, variable
        """                      
        w = self.request.get_range('within', min_value=1, default=0)
        p = self.request.get('pivot', None)
        variable = self.request.get('variable', None)

        xy = self.request.get('xy', None) 
        k = self.request.get('k', None)  
        v = self.request.get('v', None) 
        c = 'true' == self.request.get('c') 
        si = 'true' == self.request.get('si')
        bb = self.request.get('bb', None)
        bb_offset = self.request.get('bb_offset', None)
        limit = self.request.get_range('limit', min_value=1, max_value=100, default=10)
        offset = self.request.get_range('offset', min_value=0, default=0)

        # Invalid request
        if not k and not xy and not bb and not (w and p):
            self.error(404)
            return
        
        # Get cell key unqiues
        cell_keys = set()
        offset_key = None

        if w and p and variable: # If range query ignore bb, k, xy params
            cell_keys = CellIndex.search(variable, w, p, limit, offset)
        elif bb: # If bb then ignore other cell key sources (e.g., k and xy)
            nw,se = bb.split('|')
            w,n = nw.split(',')
            e,s = se.split(',')
            nwpoint = rmg.Point(float(w), float(n)) # lon,lat
            sepoint = rmg.Point(float(e), float(s)) # lon,lat
            count = 0
            for cell_key in rmg.RMGCell.cells_in_bb(nwpoint, sepoint, startkey=bb_offset):
                if count == limit:
                    offset_key = cell_key
                    break
                count += 1
                cell_keys.add(cell_key)
        else: 
            if k: 
                cell_keys.update([x.strip() for x in k.split(',')])
            if xy:
                cell_keys.update(self.cell_keys_from_coords([x.strip() for x in xy.split('|')]))

        # Invalid request
        if not cell_keys:
            logging.error('No cell keys for k=%s, xy=%s' % (k, xy))
            self.error(404)
            return
        
        # Get variable names
        if v:
            variable_names = set([x.strip() for x in v.split(',')])
        elif variable:
            variable_names = [variable]
        else:
            variable_names = []
        
        # Get cells by key
        cells = CellValuesHandler.getcells(cell_keys)        

        # Prepare results
        results = []                    
        for cellkey, cell in cells.iteritems():
            varvals = simplejson.loads(cell.varvals)
            requested_varvals = {}

            # Prepare only requested variables
            if len(variable_names) > 0:
                for name in variable_names:
                    var = varvals.get(name)
                    if not var: # Ignore invalid variable names
                        continue
                    val = None
                    if si:
                        val = apply(SI_CONVERSIONS.get(name), [var])
                    else:
                        val = int(var)
                    requested_varvals[name] = val
            else: # Send back all variables
                if si:
                    requested_varvals = dict(
                        (k, apply(SI_CONVERSIONS.get(k), [v])) for k,v in varvals.iteritems())
                else:
                    requested_varvals = dict((k, int(v)) for k,v in varvals.iteritems())
            
            # Package result as a dictionary
            result = dict(
                cell_key=cellkey, 
                cell_values=requested_varvals)
            
            # Include cell coordinates
            if c:
                result['cell_coords'] = simplejson.loads(cell.coords)
                
            # Add result
            results.append(result)

        if offset_key:
            results = dict(
                offset_cell_key=offset_key,
                cells=results)
            
        # Return all results as JSON
        json = simplejson.dumps(results)
        self.response.headers["Content-Type"] = "application/json"
        self.response.out.write(json)

application = webapp.WSGIApplication(
    [('/api/cells/values', CellValuesHandler),], debug=True)    

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
