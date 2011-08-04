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

class CouchDbCell(db.Model):
    """Models a CouchDB cell document.

    key_name - The cell key (e.g., 1-2).
    """
    rev = db.StringProperty(required=True, indexed=False)
    coords = db.StringProperty(required=True, indexed=False)
    varvals = db.TextProperty(required=True)

    def __eq__(self, other):
        if isinstance(other, CouchDbCell):
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

class Variable(db.Expando):
    """Variable metadata."""
    name = db.StringProperty()

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
        entities = CouchDbCell.get_by_key_name(cell_keys)
        cells = {}
        for x in entities:
            if x:
                cells[x.key().name()] = x
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
            cells[key] = CouchDbCell(
                key_name=key,
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
            db.put(couched.values())
            
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
        """
        xy = self.request.get('xy', None) 
        k = self.request.get('k', None)  
        v = self.request.get('v', None) 
        c = 'true' == self.request.get('c') 

        # Invalid request
        if not k and not xy:
            self.error(404)
            return
        
        # Get cell key unqiues
        cell_keys = set()
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
                    requested_varvals[name] = varvals.get(name)
            else: # Send back all variables
                requested_varvals = varvals
            
            result = dict(
                cell_key=cellkey, 
                cell_values=requested_varvals)
            
            # Send back cell coordinates
            if c:
                result['cell_coords'] = simplejson.loads(cell.coords)
                
            # Add result
            results.append(result)
            
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
