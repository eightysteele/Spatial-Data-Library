'''
Helper module containing transformation functions for the bulkloader that
are referenced by config.yaml file.
'''

from django.utils import simplejson
from google.appengine.ext.bulkload import transform

def create_key():
    def wrapper(value, bulkload_state):
        '''Returns a CellIndex key with Cell as it's parent.'''
        d = bulkload_state.current_dictionary
        #d['pkey'] = '%s_%s' % (d['key'], d['variable'])
        d['pkey'] = d['key']
        d['ikey'] = d['variable']
        return transform.create_deep_key(('Cell', 'pkey'),
                                         ('CellIndex', 'ikey'))(value, bulkload_state)
    return wrapper

def cell_to_json():
    def wrapper(value):
        coordinates = []
        for pair in value.split():
            latlng = pair.split(',')
            coordinates.append({"lat":float(latlng[0]), "lng":float(latlng[1])})
        return simplejson.dumps(coordinates)
    return wrapper


def to_json():
    def wrapper(val, bulkload_state):
        d = bulkload_state.current_dictionary.copy()
        obj = {'variable':d['variable'],
               'cell':simplejson.loads(cell_to_json()(d['cell'])),
               'value':d['value']}
        return simplejson.dumps(obj)
    return wrapper


def within_list(within):
    def wrapper(scaledval, bulkload_state):
        '''Returns a list of values before and after the scaledvalue.'''
        x = int(scaledval)
        return range(x - within, x + within + 1)
    return wrapper

def to_list(fn):
    def wrapper(value):
        '''Evualuates the CSV value as a list and returns it. 
        
        Arguments:
            value - a string encoded list of ints like "[1, 2, 3]"
        '''
        out = None
        if value is not None and len(value) > 0 and eval(value):
            out = eval(value)
        return out
    return wrapper
