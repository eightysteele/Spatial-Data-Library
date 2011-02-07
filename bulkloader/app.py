from google.appengine.ext import webapp
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app
import os

# ===============================================================================
# Request handlers

class BaseHandler(webapp.RequestHandler):
    def render_template(self, file, template_args):
        path = os.path.join(os.path.dirname(__file__), "templates", file)
        self.response.out.write(template.render(path, template_args))

class IndexHandler(BaseHandler):
    """Handler for the sandbox UI."""
    def get(self):
        self.render_template("index.html", {})

application = webapp.WSGIApplication(
         [('/', IndexHandler), ],
         debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
