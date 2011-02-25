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
from google.appengine.api import mail, memcache as m
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
from sdl import tmg
import logging
import os
import simplejson

memcache = m.Client()

# ==============================================================================
# Datastore model

class Cell(db.Expando):
    """Triangular mesh grid cell. 
    
    The key name is the icosahedron face number along with the orthogonal 
    offsets i and j:    
    
    key_name = face-i-j
    """
    json = db.TextProperty()
    kml = db.TextProperty()

class CellIndex(db.Expando):
    """Cell index."""
    value = db.IntegerProperty()
    variable = db.StringProperty()
    metadata = db.TextProperty() # JSON .hdr info

class Variable(db.Expando):
    """Variable metadata."""
    name = db.StringProperty()

def pretty_date(time=False):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """

    import datetime as dt

    now = datetime.now()
    if type(time) is int:
        diff = now - datetime.fromtimestamp(time)
    elif not time:
        diff = now - now
    else:
        diff = now - time

    if type(diff) is dt.timedelta:
        second_diff = diff.seconds
        day_diff = diff.days
    else:
        second_diff = diff.second
        day_diff = diff.day

    if day_diff < 0:
        return ''

    if day_diff == 0:
        if second_diff < 10:
            return "just now"
        if second_diff < 60:
            return str(second_diff) + " seconds ago"
        if second_diff < 120:
            return  "a minute ago"
        if second_diff < 3600:
            return str(second_diff / 60) + " minutes ago"
        if second_diff < 7200:
            return "an hour ago"
        if second_diff < 86400:
            return str(second_diff / 3600) + " hours ago"
    if day_diff == 1:
        return "Yesterday"
    if day_diff < 7:
        return str(day_diff) + " days ago"
    if day_diff < 31:
        return str(day_diff / 7) + " weeks ago"
    if day_diff < 365:
        return str(day_diff / 30) + " months ago"
    return str(day_diff / 365) + " years ago"

def _getprops(obj):
    """Returns dictionary of all static and dynamic properties of an entity."""
    dict = {}
    for key in obj.properties().keys():
        val = obj.properties()[key].__get__(obj, CellIndex)
        if type(val) is datetime.datetime:
            dict[key] = str(val)
        else:
            dict[key] = val
    for key in obj.dynamic_properties():
        val = obj.__getattribute__(key)
        if type(val) is datetime.datetime:
            dict[key] = str(val)
        else:
            dict[key] = val
    dict['key'] = obj.key().name()
    return dict

PLACEMARK = u"""                                                                                                                                          
    <Placemark>                                                                                                                                               
        <name>Cellz</name>                                                                                                                              
        <visibility>1</visibility>                                                                                                                            
        <styleUrl>#transGreenPoly</styleUrl>                                                                                                                  
        <description><![CDATA[%s]]></description>
                        %s                                                                                                                                  
                        %s                                                                                                                                  
                        %s                                                                                                                                  
                        %s                                                                                                                                  
                    </coordinates>                                                                                                                            
                </LinearRing>                                                                                                                                 
            </outerBoundaryIs>                                                                                                                                
        </Polygon>                                                                                                                                            
    </Placemark>"""
#                        %s,1                                                                                                                                  
#                        %s,1                                                                                                                                  
#                        %s,1                                                                                                                                  
#                        %s,1                                                                                                                                  
#                        %s,1                                                                                                                                  
    
KML = u'''<?xml version="1.0" encoding="UTF-8"?>                                                                                                          
<kml xmlns="http://www.opengis.net/kml/2.2">                                                                                                                  
    <Document>                                                                                                                                                
        <name>KmlFile</name>                                                                                                                                  
        <Style id="transGreenPoly">                                                                                                                           
            <LineStyle>                                                                                                                                       
                <width>1.5</width>                                                                                                                            
                <color>11111111</color>                                                                                                                       
            </LineStyle>                                                                                                                                      
            <PolyStyle>                                                                                                                                       
                <color>7d00ff00</color>                                                                                                                       
            </PolyStyle>                                                                                                                                      
        </Style>                                                                                                                                              
        <Folder>                                                                                                                                              
            <name>Triangular Mesh Grid</name>                                                                                                                 
            <visibility>1</visibility>                                                                                                                        
            <description>Global triangular mesh grid coverage.</description>                                                                                  
            %s                                                                                                                                                
        </Folder>                                                                                                                                             
    </Document>                                                                                                                                               
</kml>'''

def createKmlMesh(cell_count):
    placemarks = []
    for n in range(10):
        for x in range(cell_count):
            for y in range(cell_count):
                polygon = tmg.Cell.polygon(n, x, y, cell_count)                    
#                points = tuple(['%s,%s' % (c[0], c[1]) for c in polygon])
                points = tuple(['%s,%s' % (c[0], c[1]) for c in polygon])
                key = '%s-%s-%s' % (n, x, y)
#                data = (key, key) + points
                data = (points, points) + points
                p = PLACEMARK % data
                placemarks.append(p)
    return KML % ' '.join(placemarks)

# ==============================================================================
# Request handlers

class BaseHandler(webapp.RequestHandler):
    """Base handler for common functions like template rendering."""
    def render_template(self, file, template_args):
        path = os.path.join(os.path.dirname(__file__), "html", file)
        self.response.out.write(template.render(path, template_args))

class MeshHandler(BaseHandler):
    def get(self, cell_count):
        format = self.request.get('format', None)
        if format == 'text':
            memcache_key = 'mesh-%s' % cell_count
            kml = memcache.get(memcache_key)
            if not kml:
                kml = createKmlMesh(cell_count)
                memcache.set(memcache_key, kml)
            self.response.out.write(kml)
        elif format == 'earth':
            self.render_template("meshearth.html", {})
        else:
            self.render_template("meshmap.html", {})
            
class KmlHandler(BaseHandler):
    def get(self, cell_count):
        cell_count = int(cell_count)
        memcache_key = 'mesh-%s' % cell_count
        kml = memcache.get(memcache_key)
        if not kml:
            kml = createKmlMesh(cell_count)
            memcache.set(memcache_key, kml)
        self.response.headers['Content-Type'] = 'application/vnd.google-earth.kml+xml'
        self.response.out.write(kml)

class CellHandler(BaseHandler):    
    def get(self, n, x, y):
        try:
            cell_count = self.request.get('cc', None)
            if cell_count is not None:
                cell_count = int(cell_count)
            polygon = tmg.Cell.polygon(int(n), int(x), int(y), cell_count)
        except (Exception), e:
            logging.error(str(e))
            self.error(404)
            return
        self.response.headers['Content-Type'] = 'application/json'
        self.response.out.write(simplejson.dumps(polygon))
                                
class GitHubPostReceiveHooksHandler(BaseHandler):    
    def post(self):
        payload = self.request.get('payload')
        json = simplejson.loads(payload)
        title = '[%s] New GitHub activity - GIT push' % json['repository']['name']
        body = 'The following commits were just pushed:\n\n'
        for c in json['commits']:
            body += '%s\n' % c['message']
            body += '%s (author)\n' % c['author']['name']
            body += '%s\n' % c['timestamp']
            body += '%s\n\n' % c['url']
        logging.info(body)
        mail.send_mail(sender="Spatial Datastore Library <admin@geo-ds.appspotmail.com>",
              to="Aaron <eightysteele@gmail.com>, John <tuco@berkeley.edu>, Dave <dave.vieglais@gmail.com>",
              subject=title,
              body=body)        
        
class DataListHandler(BaseHandler):
    '''Lists all Variable entities as JSON.'''
    def get(self):
        self.response.out.write(simplejson.dumps([{_getprops(x)['key']: _getprops(x)} for x in Variable.all()]))

class DataHandler(BaseHandler):
    '''Lists metadata for a single Variable entity as JSON.'''
    def get(self, variable):
        entity = Variable.get_by_key_name(variable)
        if entity is None:
            self.error(404)
            return
        self.response.out.write(simplejson.dumps(_getprops(entity)))


class ApiHandler(BaseHandler):
    """API handler."""
    
    def searchLatLng(self, ll):
        '''Searches for Cell associated with lat/lng and returns Cell.json.'''
        cell = memcache.get(ll)
        if cell:
            self.response.out.write(cell.bio1)
            return;
        try:
            lat, lng = ll.replace(' ', '').split(',')
            key = tmg.Rhomboid(float(lat), float(lng)).key
        except (Exception), e:
            logging.error(str(e))
            self.error(400)
            return
        cell = Cell.get_by_key_name(key)
        logging.info('key=%s, cell=%s' % (key, cell))
        if not cell:
            self.error(404)
            return
        memcache.add(ll, cell)
        self.response.out.write(cell.bio1)
    
    def get(self):
        self.post()
             
    def post(self):
        ll = self.request.get('ll', None)
        if ll:
            self.searchLatLng(ll)
            return
        
        within = self.request.get('within', None)
        variable = self.request.get('variable', None)
        value = int(self.request.get('pivot', None))

        within_filter = 'within_%s =' % within

        query = db.Query(CellIndex, keys_only=True)
        query.filter('variable =', variable).filter(within_filter, value)
        index = query.get()
        if not index:
            self.error(404)
            return
        cell = db.get(index.parent())
        self.response.out.write(cell.__getattribute__(variable))

application = webapp.WSGIApplication(
        [('/data/api', ApiHandler),
         ('/data/([\w]*)', DataHandler),
         ('/data', DataListHandler),
         ('/cells/([\d]+)/([\d]+)/([\d]+)', CellHandler),
         ('/cells/mesh/([\d]+)', MeshHandler),
         ('/cells/mesh/kml/([\d]+)', KmlHandler),
         ('/github/post-commit-hook', GitHubPostReceiveHooksHandler),
         ], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
        
