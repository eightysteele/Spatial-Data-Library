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
import base64
from datetime import datetime
import logging
import os
import simplejson
import sys
import urllib

# Google App Engine imports:
from google.appengine.api import mail, memcache, urlfetch
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app

# SDL imports:
from sdl import rmg, interval

# Datastore Plus imports
from ndb import query, model
from ndb.query import OR, AND

# CouchDb connection parameters:
COUCHDB_COOKIE = 'void'
COUCHDB_HOST = 'http://spatial.iriscouch.com'
COUCHDB_PORT = 5984
COUCHDB_DATABASE = 'worldclim'
COUCHDB_DESIGN = 'api'
COUCHDB_VIEW = 'cells'
COUCHDB_URL = '%s:%s/%s/_design/%s/_view/%s' % \
    (COUCHDB_HOST, 
     COUCHDB_PORT, 
     COUCHDB_DATABASE, 
     COUCHDB_DESIGN, 
     COUCHDB_VIEW)

WC_ALIAS = dict(
    tmean1='m1',
    tmean2='m2',
    tmean3='m3',
    tmean4='m4',
    tmean5='m5',
    tmean6='m6',
    tmean7='m7',
    tmean8='m8',
    tmean9='m9',
    tmean10='m10',
    tmean11='m11',
    tmean12='m12',
    tmin1='n1',
    tmin2='n2',
    tmin3='n3',
    tmin4='n4',
    tmin5='n5',
    tmin6='n6',
    tmin7='n7',
    tmin8='n8',
    tmin9='n9',
    tmin10='n10',
    tmin11='n11',
    tmin12='n12',
    tmax1='x1',
    tmax2='x2',
    tmax3='x3',
    tmax4='x4',
    tmax5='x5',
    tmax6='x6',
    tmax7='x7',
    tmax8='x8',
    tmax9='x9',
    tmax10='x10',
    tmax11='x11',
    tmax12='x12',
    prec1='p1',
    prec2='p2',
    prec3='p3',
    prec4='p4',
    prec5='p5',
    prec6='p6',
    prec7='p7',
    prec8='p8',
    prec9='p9',
    prec10='p10',
    prec11='p11',
    prec12='p12',
    alt='a',
    bio1='b1',
    bio2='b2',
    bio3='b3',
    bio4='b4',
    bio5='b5',
    bio6='b6',
    bio7='b7',
    bio8='b8',
    bio9='b9',
    bio10='b10',
    bio11='b11',
    bio12='b12',
    bio13='b13',
    bio14='b14',
    bio15='b15',
    bio16='b16',
    bio17='b17',
    bio18='b18',
    bio19='b19' )

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

class CellIndex(model.Model): # parent=Cell, key_name=varname    
    #n = model.StringProperty('n', required=True) # variable name
    #v = model.IntegerProperty('v', required=True) # variable value
    a0 = model.IntegerProperty()
    a1 = model.IntegerProperty()
    a2 = model.IntegerProperty()
    a3 = model.IntegerProperty()
    a4 = model.IntegerProperty()
    a5 = model.IntegerProperty()
    a6 = model.IntegerProperty()
    a7 = model.IntegerProperty()
    a8 = model.IntegerProperty()
    a9 = model.IntegerProperty()
    a10 = model.IntegerProperty()
    a11 = model.IntegerProperty()
    a12 = model.IntegerProperty()
    a13 = model.IntegerProperty()

    b10 = model.IntegerProperty()
    b11 = model.IntegerProperty()
    b12 = model.IntegerProperty()
    b13 = model.IntegerProperty()
    b14 = model.IntegerProperty()
    b15 = model.IntegerProperty()
    b16 = model.IntegerProperty()
    b17 = model.IntegerProperty()
    b18 = model.IntegerProperty()
    b19 = model.IntegerProperty()
    b110 = model.IntegerProperty()
    b111 = model.IntegerProperty()
    b112 = model.IntegerProperty()
    b113 = model.IntegerProperty()

    b120 = model.IntegerProperty()
    b121 = model.IntegerProperty()
    b122 = model.IntegerProperty()
    b123 = model.IntegerProperty()
    b124 = model.IntegerProperty()
    b125 = model.IntegerProperty()
    b126 = model.IntegerProperty()
    b127 = model.IntegerProperty()
    b128 = model.IntegerProperty()
    b129 = model.IntegerProperty()
    b1210 = model.IntegerProperty()
    b1211 = model.IntegerProperty()
    b1212 = model.IntegerProperty()
    b1213 = model.IntegerProperty()


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

    COUCHDB_COOKIE = None

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
    def setcouchcookie(cls):
        # Set new cookie
        logging.info('BAD COOKIE -- getting new one')
        u,p = open('couchdb-creds.txt', 'r').read().split(':')
        payload = urllib.urlencode(dict(name=u.strip(), password=p.strip()))
        url = 'http://spatial.iriscouch.com:5984/_session'
        r = urlfetch.fetch(
            follow_redirects=False,
            url=url,
            payload=payload,
            method=urlfetch.POST,
            headers={'Content-Type': 'application/x-www-form-urlencoded'})
        cls.COUCHDB_COOKIE = r.headers['set-cookie']
        logging.info('New cookie %s' % cls.COUCHDB_COOKIE)

    @classmethod
    def fromcouchdb(cls, cell_keys):
        """Returns a list of CouchDBCell entities from a CouchDB query on cell keys.
        
        Arguments:
            cell_keys - A set of cell key strings (e.g., 9-15).

        Returns:
            A dictionary of cell key to CouchDBCell instance.
        """
        logging.info('COUCHDB=%s' % COUCHDB_URL)

        if not cls.COUCHDB_COOKIE:
            cls.setcouchcookie()

        try:
            response = urlfetch.fetch(
                url=COUCHDB_URL,
                payload=simplejson.dumps({'keys': list(cell_keys)}),
                method=urlfetch.POST,
                headers={"Content-Type":"application/json",
                         "X-CouchDB-WWW-Authenticate": "Cookie",
                         "Cookie": cls.COUCHDB_COOKIE})
            logging.info('CONTENT=%s' % response.content)

            if response.status_code == 401:
                # Retry with new cookie since it expired
                cls.setcouchcookie()            
                response = urlfetch.fetch(
                    url=COUCHDB_URL,
                    payload=simplejson.dumps({'keys': list(cell_keys)}),
                    method=urlfetch.POST,
                    headers={"Content-Type":"application/json",
                             "X-CouchDB-WWW-Authenticate": "Cookie",
                             "Cookie": cls.COUCHDB_COOKIE})
                logging.info('CONTENT=%s' % response.content)

            elif response.status_code != 200:
                return {}

        except Exception as e:
            logging.error(e)
            return {}

        cells = {}        
        for row in simplejson.loads(response.content).get('rows'):            
            key = row.get('key')
            value = row.get('value')
            logging.info('VALUE=%s' % value)
            cells[key] = Cell(
                id=key,
                rev=value.get('rev'),
                coords=simplejson.dumps(value.get('b')),
                varvals=simplejson.dumps(value.get('v')))
        return cells

    @classmethod
    def getcells(cls, cell_keys):
        """Gets CouchDBCell entities corresponding to a set of cell keys.
        
        Arguments:
            cell_keys - A set of cell key strings (e.g., 9-15).

        Returns:
            A dictionary of cell key to CouchDBCell.
        """
        logging.info(cell_keys)
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

    def process_cell_keys(self, cell_keys, si=None, c=None, offset_key=None, variable_names = []):
        if not c:
            c = 'true' == self.request.get('c') 
        if not si:
            si = 'true' == self.request.get('si')

        #logging.info('cell_keys=%s' % cell_keys)

        # Get cells by key
        cells = CellValuesHandler.getcells(cell_keys)        
        
        # Prepare results
        results = []                    
        for cellkey, cell in cells.iteritems():
            logging.info(cell.varvals)
            varvals = simplejson.loads(cell.varvals)
            requested_varvals = {}
            
            # Prepare only requested variables
            if len(variable_names) > 0:
                for name in variable_names:
                    var = varvals.get(WC_ALIAS[name], None)
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

    def range_query(self, ranges, limit, offset): #gte, lt, var, limit, offset):
        variables = []

        if len(ranges) > 1:
            qry = "CellIndex.query(AND"
        else:
            qry = "CellIndex.query"

        for r in ranges:
            var = r[0]
            variables.append(var)
            gte = int(r[1])
            lt = int(r[2])

            if var == 'alt':
                var_min = -454
                var_max = 8550
            elif var == 'bio1':
                var_min = -269
                var_max = 314
            elif var == 'bio12':
                var_min = 0
                var_max = 9916
            else:
                self.error(404)
                return
            
            var = WC_ALIAS.get(r[0])
            intervals = interval.get_query_intervals(var_min, var_max, gte, lt)
            logging.info('var=%s, gte=%s, lt=%s' % (var, gte, lt))
                
            # Build the query
            qry = "%s(AND(OR(" % qry
            #qry = "CellIndex.query(AND(CellIndex.n == '%s', OR(" % var
            for index,value in intervals.iteritems():
                if not value or not index.startswith('i'):
                    continue
                index = index.replace('i', var)
                logging.info('index=%s, value=%s' % (index, value))
                qry = '%sCellIndex.%s == %d,' % (qry, index, value)
        
        
        qry = '%s)))' % qry[:-1]
        logging.info('qry=%s' % qry)
        qry = eval(qry)

        logging.info(qry)
        results = qry.fetch(limit, offset=offset, keys_only=True)
        cell_keys = set([key.parent().id() for key in results])
        self.process_cell_keys(cell_keys, variable_names=variables)

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
            gte = greater than or equal value for range query
            lt = less than value for range query
            variable = single variable name
        """ 
        RANGE_DEFAULT = sys.maxint - 1

        # Get request params
        xy = self.request.get('xy', None) 
        k = self.request.get('k', None)  
        v = self.request.get('v', None) 
        c = 'true' == self.request.get('c') 
        si = 'true' == self.request.get('si')
        bb = self.request.get('bb', None)
        bb_offset = self.request.get('bb_offset', None)
        limit = self.request.get_range('limit', min_value=1, max_value=100, default=10)
        offset = self.request.get_range('offset', min_value=0, default=0)

        ranges = [r.split(',') for r in self.request.get_all('r')]
        logging.info('ranges=%s' % ranges)

        for r in ranges:
            var = r[0]
            if not WC_ALIAS.get(var, None):
                logging.error('Unknown variable names %s' % var)
                self.error(404)
                return
        if v:
            unknowns = [x for x in v.split(',') if not WC_ALIAS.get(x, None)]
            if len(unknowns) > 0:
                logging.error('Unknown variable names %s' % unknowns)
                self.error(404)
                return
        
        # Handle range query and return
        #if gte is not RANGE_DEFAULT and lt is not RANGE_DEFAULT and variable:
        if len(ranges) > 0:
            self.range_query(ranges, limit, offset) #gte, lt, variable, limit, offset)
            return

        # Invalid request
        if not k and not xy and not bb:
            self.error(404)
            logging.error('No k, xy, or bb')
            return
        
        # Get cell key unqiues
        cell_keys = set()
        offset_key = None

        if bb: # If bb then ignore other cell key sources (e.g., k and xy)
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

        # Invalid request since no cell keys
        if not cell_keys:
            logging.error('No cell keys for k=%s, xy=%s' % (k, xy))
            self.error(404)
            return
        
        # Get variable names for response
        if v:
            variable_names = set([x.strip() for x in v.split(',')])
        else:
            variable_names = []
            
        # Process cell keys and return results
        self.process_cell_keys(cell_keys, si=si, c=c, offset_key=offset_key, 
                               variable_names=variable_names)

application = webapp.WSGIApplication(
    [('/api/cells/values', CellValuesHandler),], debug=True)    

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
