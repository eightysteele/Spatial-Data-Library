from django.utils import simplejson
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
import os
import logging

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

# ==============================================================================
# Request handlers

class BaseHandler(webapp.RequestHandler):
    """Base handler for common functions like template rendering."""
    def render_template(self, file, template_args):
        path = os.path.join(os.path.dirname(__file__), "templates", file)
        self.response.out.write(template.render(path, template_args))

class SandboxHandler(BaseHandler):
    """Handler for the sandbox stuff."""

    def _getprops(self, obj):
        """Returns dictionary of all static and dynamic properties of an entity."""
        dict = {}
        for key in obj.properties().keys():
            dict[key] = obj.properties()[key].__get__(obj, CellIndex)
        for key in obj.dynamic_properties():
            logging.info(key)
            dict[key] = obj.__getattribute__(key)
        dict['key'] = str(obj.key())
        return dict

    def get(self):
        """Prototyping cell API stuff..."""

        action = self.request.get('action', None)
        if not action:
            self.error(404)

        if action == 'search':
            # Handles something like:
            # ?action=search&within=1&value=9&variable=bio1

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
            self.response.out.write(simplejson.dumps(self._getprops(cell)))

        if action == 'create_test_data':
            """Just creates a test entity in the datastore for bio1.
                
            Variable  = BIO1 = Annual Mean Temperature
            MinValue -50
            MaxValue 293
            """
            val = 10
            cell = Cell(key_name='f-i-j', json=simplejson.dumps({'variable':'bio1', 'value':10}))
            cell_key = cell.put()
            cell_index = CellIndex(parent=cell_key, variable='bio1', value=val)
            cell_index.within_1 = [9, 10, 11]
            cell_index.within_10 = range(val - 10, val + 10 + 1)
            cell_index.put()
            self.response.out.write(simplejson.dumps(self._getprops(cell_index)))

application = webapp.WSGIApplication([('/', SandboxHandler), ], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
