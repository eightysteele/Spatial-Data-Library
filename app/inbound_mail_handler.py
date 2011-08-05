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

"""This module handles listserve type emails."""

import email
import logging
import simplejson
import urllib

from google.appengine.api import mail
from google.appengine.api import urlfetch
from google.appengine.ext import webapp 
from google.appengine.ext.webapp.mail_handlers import InboundMailHandler 
from google.appengine.ext.webapp.util import run_wsgi_app

AUTHORIZED_SENDERS = ['eightysteele@gmail.com', 'gtuco.btuco@gmail.com', 'noreply@googlegroups.com']

class EmailHandler(InboundMailHandler):
            
    def receive(self, msg):
        logging.info('Received mail from %s' % msg.sender)
        if msg.sender not in AUTHORIZED_SENDERS:
            return
        
        # Forwards messages from google groups to Aaron:
        if msg.sender == 'noreply@googlegroups.com':
            logging.info('Forwarding message to eighty')
            mail.send_mail(sender='noreply@geo-ds.appspotmail.com',
                           to='eightysteele@gmail.com',
                           subject=msg.subject,
                           body=msg.body)

application = webapp.WSGIApplication([EmailHandler.mapping()], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
