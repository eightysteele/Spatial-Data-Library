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

import logging
import time
import math
from optparse import OptionParser
import httplib
import json
import simplejson
import couchdb
from rmg import *

def prettyPrint(s):
    """Prettyprints the json response of an HTTPResponse object"""

    # HTTPResponse instance -> Python object -> str
    print simplejson.dumps(json.loads(s.read()), sort_keys=True, indent=4)

def getJson(s):
    """Prettyprints the json response of an HTTPResponse object"""

    # HTTPResponse instance -> Python object -> str
    return json.loads(s.read())

class Couch:
    """Basic wrapper class for operations on a couchDB"""

    def __init__(self, host, port=5984, options=None):
        self.host = host
        self.port = port

    def connect(self):
        return httplib.HTTPConnection(self.host, self.port) # No close()

    # Database operations

    def createDb(self, dbName):
        """Creates a new database on the server"""

        r = self.put(''.join(['/',dbName,'/']), "")
        prettyPrint(r)

    def deleteDb(self, dbName):
        """Deletes the database on the server"""

        r = self.delete(''.join(['/',dbName,'/']))
        prettyPrint(r)

    def listDb(self):
        """List the databases on the server"""

        prettyPrint(self.get('/_all_dbs'))

    def infoDb(self, dbName):
        """Returns info about the couchDB"""
        r = self.get(''.join(['/', dbName, '/']))
        prettyPrint(r)

    # Document operations

    def listDoc(self, dbName):
        """List all documents in a given database"""

        r = self.get(''.join(['/', dbName, '/', '_all_docs']))
        prettyPrint(r)

    def listViewDocs(self, dbName, viewName):
        """List all documents in a given database"""

        r = self.get(''.join(['/', dbName, '/', viewName]))
        return r

    def getDoc(self, dbName, docId):
        """Return a document in a given database"""
        r = self.get(''.join(['/', dbName, '/', docId,]))
        return getJson(r)
    
    def openDoc(self, dbName, docId):
        """Open a document in a given database"""
        r = self.get(''.join(['/', dbName, '/', docId,]))
        prettyPrint(r)

    def saveDoc(self, dbName, body, docId=None):
        """Save/create a document to/in a given database"""
        if docId:
            r = self.put(''.join(['/', dbName, '/', docId]), body)
        else:
            r = self.post(''.join(['/', dbName, '/']), body)
        prettyPrint(r)

    def deleteDoc(self, dbName, docId, rev_id):
        r = self.delete(''.join(['/', dbName, '/', docId, '?rev=', rev_id]))
        prettyPrint(r)

    # Basic http methods

    def get(self, uri):
        c = self.connect()
        headers = {"Accept": "application/json"}
        c.request("GET", uri, None, headers)
        return c.getresponse()

    def post(self, uri, body):
        c = self.connect()
        headers = {"Content-type": "application/json"}
        c.request('POST', uri, body, headers)
        return c.getresponse()

    def put(self, uri, body):
        c = self.connect()
        if len(body) > 0:
            headers = {"Content-type": "application/json"}
            c.request("PUT", uri, body, headers)
        else:
            c.request("PUT", uri, body)
        return c.getresponse()

    def delete(self, uri):
         logging.info('Deleting uri %s.' % (uri))
         c = self.connect()
         c.request("DELETE", uri)
         return c.getresponse()

class DeletedRecords(object):
    def __init__(self, conn, couch):
        self.conn = conn
        self.couch = couch
        self.deletesql = 'delete from %s where recguid=?' % CACHE_TABLE
        self.deltasql = 'SELECT * FROM %s LEFT OUTER JOIN %s USING (recguid) WHERE %s.recguid is null' \
            % (CACHE_TABLE, TMP_TABLE, TMP_TABLE)

    def _deletechunk(self, cursor, recs, docs):
        logging.info('%s deleted' % self.totalcount)
        cursor.executemany(self.deletesql, recs)
        self.conn.commit()
        for doc in docs:
            self.couch.delete(doc)
        
    def execute(self, chunksize):
        logging.info('Checking for deleted records')
        
        cursor = self.conn.cursor()
        deletes = cursor.execute(self.deltasql)
        count = 0
        self.totalcount = 0
        docs = []
        recs = []

        for row in deletes.fetchall():            
            if count >= chunksize:
                self._deletechunk(cursor, recs, docs)
                self.totalcount += count
                count = 0
                docs = []
                recs = []
            count += 1
            recguid = row[0]
            recs.append((recguid,))
            docid = row[3]
            docrev = row[4]
            doc = {'_id': docid, '_rev': docrev}
            docs.append(doc)
        
        if count > 0:
            self.totalcount += count
            self._deletechunk(cursor, recs, docs)

        self.conn.commit()

        logging.info('DELETE: %s records deleted' % self.totalcount)
        
def execute(options):
    chunksize = int(options.chunksize)
    couch = couchdb.Server(options.couchurl)['vertnet']
    # Handles deleted records:
    DeletedRecords(conn, couch).execute(chunksize)
    
def _getoptions():
    """Parses command line options and returns them."""
    parser = OptionParser()
    parser.add_option("-c", "--command", dest="command",
                      help="SDL command",
                      default=None)
    parser.add_option("-k", 
                      "--documentkey", 
                      dest="documentkey",
                      help="The key to the document in couchdb.",
                      default=None)
    parser.add_option("-d", 
                      "--database", 
                      dest="database",
                      help="The database in couchdb.",
                      default=None)
    parser.add_option("-v", 
                      "--view", 
                      dest="view",
                      help="The view in couchdb.",
                      default=None)
    parser.add_option("-u", 
                      "--couchurl", 
                      dest="couchurl",
                      help="The CouchDB URL.",
                      default=None)
    return parser.parse_args()[0]

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    options = _getoptions()
    command = options.command.lower()
    
    if command == 'deldoc':
        server = Couch(options.couchurl)
        database = options.database
        documentkey = options.documentkey
        r = server.getDoc(database, documentkey)
        rev_id = r['_rev']
        logging.info('Deleting document %s in %s on %s.' % (documentkey, database, options.couchurl))
        
        server.deleteDoc(database, documentkey, rev_id)
        logging.info('Finished deleting document %s.' % (documentkey))

    if command == 'cleanworldclim':
        """Deletes the documents in the given server/database/view."""
        logging.info('Beginning cleanworldclim on %s.' % (options.database+'/'+options.view) )
        server = couchdb.Server(options.couchurl)
#        database = server[options.database]
        database = server[options.database]
        count = 0
        for row in database.view(options.view):
            doc = dict(_id=row.id, _rev=row.value)
            database.delete(doc)
            count += 1
            if count%1000 == 0:
                logging.info('Deleted %s documents' % (str(count)) )
        logging.info('Finished deleting documents from view %s.' % (options.database+'/'+options.view))
  
    if command == 'updatetest':
        """ Use to time the bulkloading of documents to a Couchdb database. Run cleanworldclim to remove documents added with the updatetest command."""
        server = couchdb.Server(options.couchurl)
        cdb = server[options.database]
        # Prepare a batch of cells to bulkload.
        cells = {}
        wcvars = ["tmax12","tmax11","tmax10","bio15","prec9""bio10","bio11","bio12","bio13","bio14","bio4","bio5","bio17","bio18","bio19","prec11","prec10","alt","bio16","tmin5","tmean6","tmin11","tmin10","tmin12","tmean10","tmean11","tmin6","bio6","prec1","tmean3","tmean2","tmean1","prec8","tmean7","tmin3","tmean5","tmean4","prec3","prec2","tmean9","tmean8","prec7","prec6","prec5","prec4","tmean12","tmax5","tmax4","tmax3","tmax8","prec12","tmax2","bio2","tmin4","tmin7","bio1","tmin1","bio7","tmax9","tmin2","tmax7","tmax6","bio8","bio9","tmin9","tmin8","tmax1","bio3"]
        for i in range(25000):
            cellkey = str(i) + "-" + str(i)
            cells[cellkey] = {
                '_id': cellkey, 
                'coords': RMGCell.polygon('0-0', 120),
                'vars': {}
                }
            for varname in wcvars:
                cells.get(cellkey).get('vars')[varname] = "0"
        
        # bulkload the batch
        t0 = time.time()
        logging.info('Beginning bulkload.')
        cdb.update(cells.values())        
        t1 = time.time()
        logging.info('%s documents bulkloaded in %s' % (len(cells), t1-t0))
        """Results:
            with views api/cells and sdl/zerovals:
                INFO:root:25000 documents bulkloaded in 82.1581799984
                INFO:root:25000 documents bulkloaded in 86.815792799
                INFO:root:25000 documents bulkloaded in 89.949283123
        """

    if command == 'testsave':
        server = Couch(options.couchurl)
        database = options.database
        documentkey = options.documentkey
        doc = """
        {
            "vars":
            {
                "tmax12":"0",
                "tmin12":"0"
            }
        }
        """
        logging.info('Saving document %s in %s on %s as:\n%s.' % (documentkey, database, server, doc))
        server.saveDoc(database, doc, documentkey)
        logging.info('Finished saving document %s.' % (documentkey))
