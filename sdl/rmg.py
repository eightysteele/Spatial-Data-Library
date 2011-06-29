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

"""
This module supports the conversion of geographic coordinates to a
global, equal-area Rectangular Mesh grid, which conserves cell size
and minimizes the number of cells required to cover the area of
the surface of an ellipsoid.

TODO: The Rectangular Mesh Grid design is specified here:
http://goo.gl/0awGK
"""

import logging
import math
from optparse import OptionParser

"""The number of cells in a one degree of longitude at the equator."""
CELLS_PER_DEGREE = 120

"""
DEGREE_DIGITS is the number of significant digits to the right of the decimal
to use in latitude and longitude equality determination and representation. This 
should be set to 7 to preserve reversible transformations between coordinate systems 
down to a resolution of roughly 1 m.
""" 
DEGREE_DIGITS = 7

"""
SEMI_MAJOR_AXIS is the radius of the sphere at the equator for the WGS84 datum. 
Cell construction and lookup are based on projection of the geographic coordinates 
onto a sphere of this radius. As a result, edge lengths and areas of cells projected
onto the WGS84 ellipsoid will vary slightly by latitude.
"""
SEMI_MAJOR_AXIS = 6378137.0

"""
Flattening is the measure of the oblateness of an ellipsoid. Inverse flattening is 
a parameter describing an ellipsoid for a datum. Here the default value is for WGS84.
"""
INVERSE_FLATTENING = 298.257223563

PLACEMARK_1 = u"""                                                                                                                                          
    <Placemark>                                                                                                                                               
        <name>Cell</name>                                                                                                                              
        <visibility>1</visibility>                                                                                                                            
        <styleUrl>#transGreenPoly</styleUrl>                                                                                                                  
        <description><![CDATA[%s %s]]></description>                                                                                                             
        <Polygon id="%s">                                                                                                                                     
            <outerBoundaryIs>                                                                                                                                 
                <LinearRing>                                                                                                                                  
                    <coordinates>
"""                                                                                                                             

PLACEMARK_2 = u"""
                    </coordinates>                                                                                                                            
                </LinearRing>                                                                                                                                 
            </outerBoundaryIs>                                                                                                                                
        </Polygon>                                                                                                                                            
    </Placemark>"""

KML = u"""<?xml version="1.0" encoding="UTF-8"?>                                                                                                          
<kml xmlns="http://www.opengis.net/kml/2.2">                                                                                                                  
    <Document>                                                                                                                                                
        <name>KmlFile</name>                                                                                                                                  
        <Style id="transGreenPoly">                                                                                                                           
            <LineStyle>                                                                                                                                       
                <width>1.5</width>                                                                                                                            
                <color>11111111</color>                                                                                                                       
            </LineStyle>                                                                                                                                      
            <PolyStyle>                                                                                                                                       
                <color>7d00ff00</color>                                                                                                                       
            </PolyStyle>                                                                                                                                      
        </Style>                                                                                                                                              
        <Folder>                                                                                                                                              
            <name>Triangular Mesh Grid</name>                                                                                                                 
            <visibility>1</visibility>                                                                                                                        
            <description>Global triangular mesh grid coverage.</description>                                                                                  
            %s                                                                                                                                                
        </Folder>                                                                                                                                             
    </Document>                                                                                                                                               
</kml>"""

def truncate(x, digits):
    """Returns a string representation of x including a number of places to the right of 
    the decimal equal to digits.
    
    Arguments:
        x - the input float
        digits - the number of places of precision to the right of the decimal
    """
    FORMAT = """.%sf"""
    format_x = FORMAT % digits
    return format(x, format_x)

def createPlacemark(key, polygon):
    """Returns a KML placemark for a polygon as a string.

    Arguments:
        key - the unique identifier for a cell
        polygon - the list of Points defining the boundary of the cell
    """ 
    data = (key, polygon, key)
    placemark = PLACEMARK_1 % data
    for c in polygon:
        point = '                        %s,%s,1\n' % (c[0], c[1])
        placemark = placemark + point
    placemark = placemark + PLACEMARK_2
    return placemark

class Point(object):
    """A degree-based geographic coordinate independent of a coordinate reference system."""

    def __init__(self, lng, lat):
        self._lng = lng
        self._lat = lat

    def get_lng(self):
        return self._lng
    lng = property(get_lng)

    def get_lat(self):
        return self._lat
    lat = property(get_lat)

    def isvalid(self):
        if math.fabs(self.lat) <= 90:
            if math.fabs(self.lng) <= 180:
                return True
        return False

    def __str__(self):
        return str(self.__dict__)

class CellPolygon(object):
    """A polygon for a cell uniquely identifiable by a key."""

    def __init__(self, cellkey, polygon):
        self._cellkey = cellkey
        self._polygon = polygon
        self._hashcode = hash((self._cellkey, self._polygon))

    def getcellkey(self): 
        return self._cellkey
    cellkey = property(getcellkey)

    def getpolygon(self):
        return self._polygon
    polygon = property(getpolygon)

    def __str__(self):
        return str(self.__dict__)

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __eq__(self, other):
        if isinstance(other, CellPolygon):
            return self._hashcode == other._hashcode
        return NotImplemented

    def __hash__(self):
        return self._hashcode

    def __cmp__(self, other):
        if self.cellkey > other.cellkey:
            return 1
        elif self.cellkey < other.cellkey:
            return - 1
        return 0

class RMGCell(object):
    """
    RMGCell is the grid cell of an adaptive Rectangular Mesh Grid (RMG) over an ellipsoid. The cell is 
    identified by the x and y indexes of its location with respect to the north pole and the meridian
    at -180 longitude (index 0,0 is at -180,90). All cells are defined by the coordinate reference system
    and the resolution in cells per degree. The angular y component of the cell is a constant 
    (1/cells per degree), while the angular x component is a function of latitude in order to preserve
    equal area for all cells. Exceptions to the equal area may occur at the poles, depending on the 
    cells per degree. 
    """ 

    @staticmethod
    def distances_per_degree(lat, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """ Returns a tuple containing arc distances in the units of the semi-major axis. 
        The first is the distance along an arc of one degree at a constant latitude (lat).
        The second is the distance along an arc of one degree at a constant lng centered on the given lat.

        Arguments:
            lat - the latitude at which the distances apply
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
        """

        f = 1.0 / inverse_flattening    
        e_squared = f * (2.0 - f) # e^2 = 2f - f^2

        # N - radius of curvature in the prime vertical, (tangent to ellipsoid at latitude)
        # N(lat) = a/(1-e^2*sin^2(lat))^0.5
        N = a / math.sqrt(1.0 - e_squared * (math.pow(math.sin(lat * math.pi / 180.0), 2.0))) 

        # M - radius of curvature in the prime meridian, (tangent to ellipsoid at latitude)
        # M(lat) = a(1-e^2)/(1-e^2*sin^2(lat))^1.5
        M = a * (1.0 - e_squared) / math.pow(1.0 - e_squared * math.pow(math.sin(lat * math.pi / 180.0), 2.0), 1.5)

        # longitude is irrelevant for the calculations to follow so simplify by using longitude = 0, so Y = 0
        # X = Ncos(lat)cos(long). long = 0, so cos(long) = 1.0
        X = N * math.cos(lat * math.pi / 180.0) * 1.0 
        lngdistanceperdegree = math.pi * X / 180.0 
        latdistanceperdegree = math.pi * M / 180.0
        return (lngdistanceperdegree, latdistanceperdegree)

    @staticmethod
    def column_count(lat, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        key = RMGCell.key(-180, lat, cells_per_degree, a, inverse_flattening)
        indexes = key.split('-')
        y_index = int(indexes[1])
        e = RMGCell.east(0, y_index, cells_per_degree, a, inverse_flattening)
        width = e + 180
        columns = math.ceil(360/width)
        return columns        

    @staticmethod
    def key(lng, lat, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the unique identifier for a cell given a longitude, latitude, and parameters 
        of the ellipsoid on which this coordinate occurs.

        Arguments:
            lng - the longitude of a Point for which the enclosing cell is sought
            lat - the latitude of a Point for which the enclosing cell is sought
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
            cells_per_degree - the desired resolution of the grid
        """
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        area = RMGCell.area(cells_per_degree, a)
        y_dist = latdpd / cells_per_degree # y side of cell in units of a
        y_index = math.floor((90.0 - lat) * cells_per_degree)
        if lat <= -90:
            y_index -= 1
        x_dist = area / y_dist
        x_angle = x_dist / lngdpd # x side in degrees
        if x_angle >= 360:
            x_index = 0
        else:
            x_index = math.floor((180.0 + lng) / x_angle)
        return '%s-%s' % (int(x_index), int(y_index))

    @staticmethod
    def north(y_index, cells_per_degree=CELLS_PER_DEGREE):
        """Returns the latitude of the north edge of the cell given a y_index and a cell resolution.

        Arguments:
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
        """
        lat = 90.0 - float(y_index) / cells_per_degree
        return lat 

    @staticmethod
    def south(y_index, cells_per_degree=CELLS_PER_DEGREE):
        """Returns the latitude of the south edge of the cell.

        Arguments:
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
        """
        lat = 90.0 - float(y_index + 1) / cells_per_degree
        return lat 

    @staticmethod
    def west(x_index, y_index, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the longitude of the west edge of the cell.

        Arguments:
            x_index - the zero-based index of the cell measured east from -180
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
        """
        if y_index >= (cells_per_degree * 180) - 1:
            return - 180
        lat = (RMGCell.north(y_index, cells_per_degree) + RMGCell.south(y_index, cells_per_degree)) / 2
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        area = RMGCell.area(cells_per_degree, a)
        y_dist = latdpd / cells_per_degree # y side of cell in units of a
        x_dist = area / y_dist
        x_angle = x_dist / lngdpd # x side in degrees
        lng = -180.0 + (x_index * x_angle)
        return lng

    @staticmethod
    def east(x_index, y_index, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the longitude of the east edge of the cell.

        Arguments:
            x_index - the zero-based index of the cell measured east from -180
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
        """
#        if y_index == 0:
#            return 180
        if y_index >= (cells_per_degree * 180) - 1:
            return 180
        lat = (RMGCell.north(y_index, cells_per_degree) + RMGCell.south(y_index, cells_per_degree)) / 2
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        area = RMGCell.area(cells_per_degree, a)
        y_dist = latdpd / cells_per_degree # x side of cell in units of a
        x_dist = area / y_dist
        x_angle = x_dist / lngdpd # y side in degrees
        lng = -180.0 + ((x_index + 1) * x_angle)
        return lng

    @staticmethod
    def mid_lat(y_index, cells_per_degree=CELLS_PER_DEGREE):
        """Returns the latitude of the center of the cell.

        Arguments:
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
        """
        if y_index >= (cells_per_degree * 180) - 1:
            return 0
        lat = 90.0 - (y_index / cells_per_degree) - (1.0 / (2.0 * cells_per_degree))
        return lat 

    @staticmethod
    def mid_lng(x_index, y_index, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the longitude of the center of the cell.

        Arguments:
            x_index - the zero-based index of the cell measured east from -180
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
        """
        if y_index >= (cells_per_degree * 180) - 1:
            return 0
        lat = (RMGCell.north(y_index, cells_per_degree) + RMGCell.south(y_index, cells_per_degree)) / 2
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        area = RMGCell.area(cells_per_degree, a)
        y_dist = latdpd / cells_per_degree # y side of cell in units of a
        x_dist = area / y_dist
        x_angle = x_dist / lngdpd # x side in degrees
        lng = -180.0 + (x_index * x_angle) + (x_angle / 2.0)
        return lng

    @staticmethod
    def width(y_index, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the width of the cell in degrees.

        Arguments:
            y_index - the zero-based index of the cell measured south from the north pole
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
        """
        lat = (RMGCell.north(y_index, cells_per_degree) + RMGCell.south(y_index, cells_per_degree)) / 2
        lngdpd, latdpd = RMGCell.distances_per_degree(lat, a, inverse_flattening)
        area = RMGCell.area(cells_per_degree, a)
        y_dist = latdpd / cells_per_degree # y side of cell in units of a
        x_dist = area / y_dist
        x_angle = x_dist / lngdpd # x side in degrees
        return x_angle

    @staticmethod
    def area(cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS):
        """Returns the area of the cell in the units of the semi-major axis.

        Arguments:
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
        """
        side = 2 * math.pi * a / 360 / cells_per_degree
        area = side * side
        return area

    @staticmethod
    def center(key, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the Point at the center of the cell having the given key.

        Arguments:
            key - the unique identifier for a cell
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter
        """
        indexes = key.split('-')
        x_index = int(indexes[0])
        y_index = int(indexes[1])
        lat = RMGCell.mid_lat(y_index, cells_per_degree)
        lng = RMGCell.mid_lng(x_index, y_index, cells_per_degree, a, inverse_flattening)
        return Point(lng, lat)

    @staticmethod
    def polygon(key, cells_per_degree=CELLS_PER_DEGREE, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns a polygon (list of Points) of the cell defined by the given key.

        Arguments:
            key - the unique identifier for a cell
            cells_per_degree - the desired resolution of the grid
            digits - the number of digits of precision to retain in the coordinates
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter
        """
        indexes = key.split('-')
        x_index = int(indexes[0])
        y_index = int(indexes[1])
        n = truncate(RMGCell.north(y_index, cells_per_degree), digits)
        s = truncate(RMGCell.south(y_index, cells_per_degree), digits)
        w = truncate(RMGCell.west(x_index, y_index, cells_per_degree, a, inverse_flattening), digits)
        e = truncate(RMGCell.east(x_index, y_index, cells_per_degree, a, inverse_flattening), digits)
        return [(n, w), (s, w), (s, e), (n, e), (n, w)]

    @staticmethod
    def tile(nwcorner, secorner, cells_per_degree=CELLS_PER_DEGREE, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns a set of polygons for cells intersecting a bounding box.

        Arguments:
            nwcorner - the starting Point in the northwest corner of the bounding box
            secorner - the ending Point in the southeast corner of the bounding box
            cells_per_degree - the desired resolution of the grid
            digits - the number of digits of precision to retain in the coordinates
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
        """ 
        cells = set()
        north = nwcorner.lat
        west = nwcorner.lng
        south = secorner.lat
        east = secorner.lng
        lng = west
        lat = north
        # key for the NW corner of the tile
        key = RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
        indexes = key.split('-')
        x_index = int(indexes[0])
        y_index = int(indexes[1])

        while lat >= south:
            while lng <= east:
                polygon = tuple([(float(x[0]), float(x[1])) for x in RMGCell.polygon(key, cells_per_degree, digits, a, inverse_flattening)])
                cells.add(CellPolygon(key, polygon))
                x_index += 1
                lng = RMGCell.east(x_index, y_index, cells_per_degree, a, inverse_flattening)
            lng = west
            y_index += 1
            lat = RMGCell.south(y_index, cells_per_degree)
        return cells

    @staticmethod
    def tile_cellcount(nwcorner, secorner, cells_per_degree=CELLS_PER_DEGREE, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns the number of cells intersecting a bounding box.

        Arguments:
            nwcorner - the starting Point in the northwest corner of the bounding box
            secorner - the ending Point in the southeast corner of the bounding box
            cells_per_degree - the desired resolution of the grid
            a - the semi-major axis of the ellipsoid for the coordinate reference system
            inverse_flattening - the inverse of the ellipsoid's flattening parameter 
                (298.257223563 for WGS84)
        """ 
        key0 = RMGCell.key(nwcorner.lng, nwcorner.lat, cells_per_degree, a, inverse_flattening)
        key1 = RMGCell.key(secorner.lng, secorner.lat, cells_per_degree, a, inverse_flattening)
        indexes = key0.split('-')
        start_y_index = int(indexes[1])
        indexes = key1.split('-')
        end_y_index = int(indexes[1])
        tile_width = secorner.lng - nwcorner.lng
        cellcount = 0
        for y_index in range(start_y_index, end_y_index):
            cell_width = RMGCell.width(y_index, cells_per_degree, a, inverse_flattening)
            cells_this_row = int(math.ceil(tile_width / cell_width))
            cellcount += cells_this_row
        return cellcount

    @staticmethod
    def tileKml(nwcorner, secorner, cells_per_degree=CELLS_PER_DEGREE, digits=DEGREE_DIGITS, a=SEMI_MAJOR_AXIS, inverse_flattening=INVERSE_FLATTENING):
        """Returns a grid of cells overlapping a bounding box, as KML."""
        placemarks = []
        cell_polygon_list = RMGCell.tile(nwcorner, secorner, cells_per_degree, digits, a, inverse_flattening)
        for cell_polygon in cell_polygon_list:
            polygon = RMGCell.polygon(cell_polygon.cellkey, cells_per_degree, digits, a, inverse_flattening)                                        
            p = createPlacemark(cell_polygon.cellkey, polygon)
            placemarks.append(p)
        return KML % ' '.join(placemarks)

    def __init__(self, x_index, y_index):
        self._x_index = x_index
        self._y_index = y_index
        self._hashcode = hash((self._x_index, self._y_index))

    def getx(self):
        return self._x_index
    x_index = property(getx)

    def gety(self):
        return self._y_index
    y_index = property(gety)

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, cell1, cell2):
        return cell1._hashcode == cell2._hashcode

    def __hash__(self):
        return self._hashcode

    def __cmp__(self, other):
        """Compares cell objects by x, then by y."""
        if self.x_index > other.x_index:
            return 1
        elif self.x_index < other.x_index:
            return - 1
        elif self.y_index > other.x_index:
            return 1
        elif self.y_index < other.y_index:
            return - 1
        else: return 0

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)    
    
    f = open('python.out', 'w')
    f.flush()
    f.close()

    parser = OptionParser()
    parser.add_option("-c", "--command", dest="command",
                      help="RMG command",
                      default=None)
    parser.add_option("-f", "--nwcorner", dest="nwcorner",
                      help="NW corner of bounding box",
                      default=None)
    parser.add_option("-t", "--secorner", dest="secorner",
                      help="SW corner of bounding box",
                      default=None)
    parser.add_option("-n", "--cells-per-degree", dest="cells_per_degree",
                      help="Number of cells per degree",
                      default=CELLS_PER_DEGREE)
    parser.add_option("-d", "--digits", dest="digits",
                      help="Number of digits of precision in coordinates",
                      default=DEGREE_DIGITS)
    parser.add_option("-a", "--a", dest="a",
                      help="Semi-major axis of the reference ellipsoid",
                      default=SEMI_MAJOR_AXIS)
    parser.add_option("-i", "--inverse-flattening", dest="inverse_flattening",
                      help="Inverse flattening of the ellipsoid",
                      default=INVERSE_FLATTENING)
    parser.add_option("-k", "--key", dest="key",
                      help="Unique key for a cell",
                      default=None)
    parser.add_option("-x", "--x", dest="x",
                      help="Longitude",
                      default=None)
    parser.add_option("-y", "--y", dest="y",
                      help="Latitude",
                      default=None)

    (options, args) = parser.parse_args()
    command = options.command.lower()
    
    if command == 'distances':
        lat = float(options.y)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.distances_per_degree(lat, a, inverse_flattening)
    if command == 'columns':
        lat = float(options.y)
        cells_per_degree = float(options.cells_per_degree)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.column_count(lat, cells_per_degree, a, inverse_flattening)
    if command == 'key':
        lng = float(options.x)
        lat = float(options.y)
        cells_per_degree = float(options.cells_per_degree)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.key(lng, lat, cells_per_degree, a, inverse_flattening)
    if command == 'area':
        cells_per_degree = float(options.cells_per_degree)
        a = float(options.a)
        print RMGCell.area(cells_per_degree, a)
    if command == 'center':
        key = float(options.key)
        cells_per_degree = float(options.cells_per_degree)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.center(key, cells_per_degree, a, inverse_flattening)
    if command == 'polygon':
        key = float(options.key)
        cells_per_degree = float(options.cells_per_degree)
        digits = float(options.digits)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.polygon(key, cells_per_degree, digits, a, inverse_flattening)        
    if command == 'tile':
        nwcorner = map(float, options.nwcorner.split(','))
        secorner = map(float, options.secorner.split(','))
        cells_per_degree = float(options.cells_per_degree)
        digits = float(options.digits)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.tile(nwcorner, secorner, cells_per_degree, digits, a, inverse_flattening)        
    if command == 'tilecellcount':
        nwcorner = map(float, options.nwcorner.split(','))
        secorner = map(float, options.secorner.split(','))
        cells_per_degree = float(options.cells_per_degree)
        digits = float(options.digits)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.tile_cellcount(nwcorner, secorner, cells_per_degree, digits, a, inverse_flattening)        
    if command == 'tilekml':
        nwcorner = map(float, options.nwcorner.split(','))
        secorner = map(float, options.secorner.split(','))
        cells_per_degree = float(options.cells_per_degree)
        digits = float(options.digits)
        a = float(options.a)
        inverse_flattening = float(options.inverse_flattening)
        print RMGCell.tileKml(nwcorner, secorner, cells_per_degree, digits, a, inverse_flattening)