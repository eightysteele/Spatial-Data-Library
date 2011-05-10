#!/usr/bin/env python
#
# Copyright 2011 Jante LLC and University of Kansas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from datetime import datetime
from google.appengine.api import mail, memcache, urlfetch
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
from sdl import tmg
import logging
import os
import simplejson

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
    """Gets cell values.

    Gets cell values for a list of variables that correspond to a list of 
    coordinates (or cell keys).
    
    Example:
        GET /api/cells/values?xy=1,2|7,3&v=bio1,bio2

    Returns:
    [
      {
        "coordinate": [1,2], 
        "cell-key": "1-2-3", 
        "cell-values": 
          [
            {"bio1": -50},
            {"bio2": -89}
          ]
      },
      ...
    ]
    """

    @staticmethod
    def fromcache(cell_keys):
        return set([memcache.get(x) for x in cell_keys if memcache.get(x)])

    @staticmethod
    def fromds(cell_keys):
        return set([x for x in CouchDbCell.get_by_key_name(cell_keys) if x])

    @staticmethod
    def fromcouchdb(cell_keys):
        results = set()
        url = "http://ec2-75-101-194-134.compute-1.amazonaws.com:5984/sdl/_design/api/_view/cell-values"        
        payload = simplejson.dumps({'keys': list(cell_keys)})
        logging.info('PAYLOAD: ' + payload)
        method = urlfetch.POST
        headers = {"Content-Type":"application/json"}
        result = urlfetch.fetch(url, payload=payload, method=method, headers=headers)
        if result.status_code == 200:
            logging.info('Content: ' + result.content)
            rows = simplejson.loads(result.content).get('rows')
            logging.info(rows)
            for row in rows:
                key = row.get('key')
                value = row.get('value')
                rev = value.get('rev')
                coords = simplejson.dumps(value.get('coords'))
                varvals = simplejson.dumps(value.get('varvals'))
                results.add(CouchDbCell(
                        key_name=key, rev=rev, coords=coords, varvals=varvals))                    
        else:
            logging.error(result.status_code)
        return results

    @classmethod
    def getcells(cls, cell_keys):
        """Returns a list of cell value JSON strings.
        
        cell_keys - set of CellValue key_name strings.
        """
        cells = set()
        
        # Checks cache:
        cached = cls.fromcache(cell_keys)
        cells.update(cached)
        cell_keys = cell_keys.difference(set([x.key().name() for x in cached]))

        # Checks datastore:
        if len(cell_keys) > 0:
            stored = cls.fromds(cell_keys)
            cells.update(stored)
            key_names = set()
            for cell in stored:
                key_name = cell.key().name()
                key_names.add(key_name)
                memcache.add(key_name, cell)                
            cell_keys = cell_keys.difference(key_names)
            
        # Checks CouchDB:
        if len(cell_keys) > 0:
            couched = cls.fromcouchdb(cell_keys)
            cells.update(couched)
            db.put(couched)
            for cell in couched:
                memcache.add(cell.key().name(), cell)

        return cells
        

    def get(self):
        return self.post()
    
    def post(self):
        logging.info('POSTING')
        k = self.request.get('k', None) 
        v = self.request.get('v', None)
        c = 'true' == self.request.get('c')
        if not k or not v:
            self.error(404)
            return

        cell_keys = set([x.strip() for x in k.split(',')])
        variable_names = set([x.strip() for x in v.split(',')])
        if not cell_keys or not variable_names:
            self.error(404)
            return
        
        cells = CellValuesHandler.getcells(cell_keys)        
        results = []
        for cell in cells:
            requested_varvals = {}
            varvals = simplejson.loads(cell.varvals)
            for name in variable_names:
                requested_varvals[name] = varvals.get(name)
                result = {'cell-key': cell.key().name(), 
                          'cell-values': requested_varvals}
                if c:
                    result['cell-coords'] = simplejson.loads(cell.coords)
            results.append(result)
        json = simplejson.dumps(results)
        self.response.headers["Content-Type"] = "application/json"
        self.response.out.write(json)

application = webapp.WSGIApplication(
    [('/api/cells/values', CellValuesHandler),], debug=True)    

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
