from setup_env import fix_sys_path
fix_sys_path()

# Python imports
import logging

# SDL imports
from sdl import interval

# Goole App Engine imports
from django.utils import simplejson
from google.appengine.ext.bulkload import transform

# Datastore Plus imports
from ndb import query, model

def create_key():
    def wrapper(value, bulkload_state):
        '''Returns a CellIndex key with Cell as it's parent.'''
        d = bulkload_state.current_dictionary
        d['varname'] = d['RID'].split('_')[0]
        return transform.create_deep_key(
            ('Cell', 'CellKey'),
            ('CellIndex', 'varname'))(value, bulkload_state)
    return wrapper

def create_cell_key():
    def wrapper(value, bulkload_state):
        '''Returns a CellIndex key with Cell as it's parent.'''
        key_name = bulkload_state.current_dictionary['CellKey']
        return key_name
    return wrapper


def get_varname():
    def wrapper(value, bulkload_state):
        d = bulkload_state.current_dictionary
        return d['RID'].split('_')[0]
    return wrapper

def get_varval():
    def wrapper(value, bulkload_state):
        d = bulkload_state.current_dictionary
        val = d['avg_Band1'].split('.')[0]
        return int(val)
    return wrapper

def within_list(within):
    def wrapper(scaledval, bulkload_state):
        '''Returns a list of values before and after the scaledvalue.'''
        val = scaledval.split('.')[0]
        x = int(val)
        results = [str(x) for x in range(x - within, x + within + 1)]
        logging.info(results)
        return results
    return wrapper

def get_list(within, val):
    x = int(val.split('.')[0])
    return range(x - within, x + within + 1)

def add_dynamic_properties(input_dict, instance, bulkload_state_copy):    
    """Adds dynamic properties from the CSV input_dict to the entity instance."""
    val = input_dict['avg_Band1']
    ranges = dict(
        within_1=get_list(1, val),
        within_5=get_list(5, val),
        within_10=get_list(10, val))
    for key,value in ranges.iteritems():
        instance[key] = value
    if input_dict['varname'] not in ['bio16']:
        logging.info('skipping variable %s' %  input_dict['varname'])
        return datastore.Entity('CellIndex')
    return instance


