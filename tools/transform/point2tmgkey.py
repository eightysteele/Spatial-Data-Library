#!/usr/bin/env python

def my_helper(lat, lng):
    print 'lat=%f, lng=%f' % (lat, lng)

def convert(lat, lng):
    '''Converts the lat/lng into a TGM key and returns it.
    
    Arguments:
        lat - the decimal latitude
        lng - the decimal longitude
    
    Returns:
        The TGM key as a string
    '''
    
    # You can define and call helper methods if it helps:
    my_helper(lat, lng)
    return 'face-i-j'

    
if __name__ == '__main__':    
    convert(37.323345, -4.0334334)
