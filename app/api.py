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

class Cell(db.Expando):
    """Cell model with x-y key_name.

    The coords property is a JSON string:
   
    {
      "cell-key": "1-2-3",
      "coordinates": 
        [
          [-122.2887875,37.8603069],
          [-122.2774275,37.8603069],
          [-122.2729888,37.8686002],
          [-122.2843488,37.8686002],
          [-122.2887875,37.8603069]
        ],
    }

    The values property is a JSON string:

    {
      "cell-key": "1-2-3", 
      "cell-values": 
        [
          {"bio1": -50},
          ...
          {"bioN": -89}          
        ]
    }

    """
    coords = db.TextProperty() # "[[x0,y0],...,[x4,y4]]"
    values = db.TextProperty() # "{"bio1": -50, "bio2": 0, ...}"

class CellValuePage(db.Model):
    """Models a page of at most 5k cell values in JSON. Roughly 250k
    payload since 1 value object is ~50 characters.

    
    13,064,305 cells per tile.

    key_name: variablename-within-pivot-offset

    [
      {"cell-key": "1-2-3", "cell-value": 8},
      ...
    ]
    """

class CellValue(db.Model):
    """Models a cell value for a specific variable name.

    key_name: x-y-variablename
    
    {"bio1": -50}
    """
    json = db.TextProperty(required=True)

class CellCoordinates(db.Model):
    """Models a cell key and it's coordinates as JSON. 

    key_name: x-y

    {
      "cell-key": "1-2",
      "coordinates": 
        [
          [-122.2887875,37.8603069],
          [-122.2774275,37.8603069],
          [-122.2729888,37.8686002],
          [-122.2843488,37.8686002],
          [-122.2887875,37.8603069]
        ],
    }
    """
    json = db.TextProperty(required=True)

class CellIndex(db.Expando):
    """Cell index."""
    value = db.IntegerProperty()
    variable = db.StringProperty()
    metadata = db.TextProperty() # JSON .hdr info

class Variable(db.Expando):
    """Variable metadata."""
    name = db.StringProperty()

class CellValueJson(object):
    def __init__(self, cell_key, source):
        self.cell_key = cell_key
        self.source = source
        self.cell_values = []

    def addval(self, json):
        self.cell_values.append(eval(json))

    def keepvars(self, varnames):
        results = []
        for val in self.cell_values:
            for name in varnames:
                if name in val.keys():
                    results.append(val)
        self.cell_values = results                    

    def __str__(self):
        return str(self.__dict__)

class CellValuesHandler(webapp.RequestHandler):
    """Gets cell values.

    Gets cell values for a list of variables that correspond to a list of 
    coordinates (or cell keys).

    GET /api/cells/values?xy=1,2|7,3&v=bio1,bio2

    Returns:

    [
      {
        "coordinate": [1,2], 
        "cell-key": "1-2-3", 
        "cell-values": 
          [
            {"name": "bio1", "value": -50},
            {"name": "bio2", "value": -89}
          ]
      },
      ...
    ]
    """

    @staticmethod
    def cellkey(key_name):
        return key_name[:key_name.rfind('-')]

    @staticmethod
    def varname(key_name):
        return key_name[key_name.rfind('-'):]

    @staticmethod
    def fromcache(key_names):
        results = {}
        for key_name in key_names:
            val = memcache.get(key_name)
            if val:
                results[key_name] = val
        return results

    @staticmethod
    def fromds(key_names):
        results = {}
        for val in CellValue.get_by_key_name(key_names):
            if val:
                key_name = val.key().name()
                results[key_name] = val
        return results

    @staticmethod
    def fromcouchdb(key_names):
        results = {}
        cell_keys = map(lambda key_name: CellValuesHandler.cellkey(key_name), key_names)
        url = "http://ec2-75-101-194-134.compute-1.amazonaws.com:5984/sdl/_design/api/_view/cell-values"        
        payload = simplejson.dumps({'keys': cell_keys})
        method = urlfetch.POST
        headers = {"Content-Type":"application/json"}
        result = urlfetch.fetch(url, payload=payload, method=method, headers=headers)
        if result.status_code == 200:
            rows = simplejson.loads(result.content).get('rows')
            for row in rows:
                for variable in row.get('value').keys():
                    key_name = '%s-%s' % (row.get('key'), variable)
                    value = row.get('value').get(variable)
                    results[key_name] = CellValue(
                        key_name=key_name, 
                        json=simplejson.dumps({variable: value}))
        return results

    @staticmethod
    def updatecellvals(cellvals, entities, source):
        for key_name in entities.keys():
            variable = CellValuesHandler.varname(key_name)
            cell_key = CellValuesHandler.cellkey(key_name)
            if not cellvals.has_key(cell_key):
                cellvals[cell_key] = CellValueJson(cell_key, source)
            val = simplejson.loads(entities.get(key_name).json)
            cellvals.get(cell_key).addval(entities.get(key_name).json)
        return cellvals

    @staticmethod
    def getvals(key_names):
        """Returns a list of cell value JSON strings.
        
        key_names - set of CellValue key_name strings.
        """
        cellvals = {}
        
        # Checks cache:
        cached = CellValuesHandler.fromcache(key_names)
        if not cached:
            logging.info('CACHE MISS')
        else:
            logging.info('CACHE HIT')
        cellvals = CellValuesHandler.updatecellvals(cellvals, cached, 'memcache')
        diffs = key_names.difference(set(cached.keys()))

        # Checks datastore if cache didn't have everything and updates memcache:
        if len(diffs) > 0:
            stored = CellValuesHandler.fromds(diffs)
            if not stored:
                logging.info('DATASTORE MISS')
            else:
                logging.info('DATASTORE HIT')
            cellvals = CellValuesHandler.updatecellvals(cellvals, stored, 'datastore')
            for key_name in stored.keys():
                memcache.add(key_name, stored.get(key_name))                
            diffs = diffs.difference(set(stored.keys()))
            
        # Downloads from CouchDB if cache and datastore didn't have everything
        # and updates datastore and cache:
        if len(diffs) > 0:
            couched = CellValuesHandler.fromcouchdb(diffs)
            if not couched:
                logging.info('COUCHDB MISS')
            else:
                logging.info('COUCHDB HIT')
            cellvals = CellValuesHandler.updatecellvals(cellvals, couched, 'couchdb')
            db.put(couched.values())
            for entity in couched.values():
                memcache.add(entity.key().name(), entity)

        #return simplejson.dumps(map(lambda x: x.__dict__ if, cellvals.values()))
        return cellvals.values()
        

    def get(self):
        return self.post()
    
    def post(self):
        logging.info('POSTING')
        k = self.request.get('k', None) 
        v = self.request.get('v', None)
        if not k or not v:
            self.error(404)
            return

        cell_keys = [x.strip() for x in k.split(',')]
        variable_names = [v.strip() for x in v.split(',')]
        if not cell_keys or not variable_names:
            self.error(404)
            return
        
        key_names = set()
        for varname in variable_names:
            for cellkey in cell_keys:
                key_names.add('%s-%s' % (cellkey, varname))
        
        vals = CellValuesHandler.getvals(key_names)        
        for val in vals:
            val.keepvars(variable_names)
        json = simplejson.dumps([x.__dict__ for x in vals])
        self.response.headers["Content-Type"] = "application/json"
        self.response.out.write(json)

application = webapp.WSGIApplication(
    [('/api/cells/values', CellValuesHandler),], debug=True)    

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
