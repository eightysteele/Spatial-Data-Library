# Copyright 2011 Aaron Steele
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

__author__ = "Aaron Steele"

"""This module handles post recieve service hooks from GitHub."""

import common

from datetime import datetime
import logging
import os
import simplejson

from google.appengine.api import mail
from google.appengine.ext import webapp, db
from google.appengine.ext.webapp import template
from google.appengine.ext.webapp.util import run_wsgi_app

class PostReceiveHandler(common.BaseHandler):    
    """Handler for GitHub Post Receive service hook.
    https://github.com/VertNet/Software/admin/hooks
    """
    def post(self):
        payload = self.request.get('payload')
        logging.info('GitHub post receive payload: %s' % payload)
        json = simplejson.loads(payload)
        title = '[%s] GitHub push' % json['repository']['name']
        body = 'The following commits were just pushed:\n\n'
        for c in json['commits']:
            body += '%s\n%s' % (c['message'], c['url'])
            body += '\n%s (author)' % c['author']['name']
            body += '\n%s\n\n' % c['timestamp'] # TODO: user friendly date
        mail.send_mail(
            sender="GitHub <commits@geo-ds.appspotmail.com>",
            to="sdl-developers@googlegroups.com",
            subject=title,
            body=body)     

application = webapp.WSGIApplication(
    [('/service-hooks/post-receive', PostReceiveHandler),], debug=True)        

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
