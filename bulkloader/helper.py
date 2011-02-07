'''
Helper module containing transformation functions for the bulkloader that
are referenced by config.yaml file.
'''

def toList(fn):
    '''Returns a closure that evaluates CSV value as a list.'''
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
