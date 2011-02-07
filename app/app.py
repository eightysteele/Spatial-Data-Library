from django.utils import simplejson
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
import os
import logging
import datetime

# ==============================================================================
# Datastore model

class Cell(db.Model):
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


# ==============================================================================
# Request handlers

class BaseHandler(webapp.RequestHandler):
    """Base handler for common functions like template rendering."""
    def render_template(self, file, template_args):
        path = os.path.join(os.path.dirname(__file__), "templates", file)
        self.response.out.write(template.render(path, template_args))

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

    def get(self):
        action = self.request.get('action', None)
        if not action:
            self.error(404)
            return
        if action == 'search':
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
            self.response.out.write(simplejson.dumps(_getprops(cell)))

application = webapp.WSGIApplication([('/data/api', ApiHandler),
                                      ('/data/([\w]*)', DataHandler),
                                      ('/data', DataListHandler), ], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
