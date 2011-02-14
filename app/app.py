from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
from datetime import datetime
import os
import simplejson
from google.appengine.api import mail
import logging


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


# ==============================================================================
# Request handlers

class BaseHandler(webapp.RequestHandler):
    """Base handler for common functions like template rendering."""
    def render_template(self, file, template_args):
        path = os.path.join(os.path.dirname(__file__), "templates", file)
        self.response.out.write(template.render(path, template_args))

class PostReceiveHooksHandler(BaseHandler):    
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
              to="Aaron <eightysteele@gmail.com>", #, John <tuco@berkeley.edu>, Dave <dave.vieglais@gmail.com>",
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
            
    def post(self):
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
            self.response.out.write(cell.__getattribute__(variable))

application = webapp.WSGIApplication([('/data/api', ApiHandler),
                                      ('/data/([\w]*)', DataHandler),
                                      ('/data', DataListHandler),
                                      ('/github/post-commit-hook', PostReceiveHooksHandler),
                                      ], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
