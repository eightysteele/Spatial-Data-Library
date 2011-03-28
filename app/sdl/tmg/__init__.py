#!/usr/bin/env python
#
# Copyright 2011 Jante LLC and University of Kansas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import logging
import math

'''Set CELL_COUNT to the number of cells on each side of the rhomboid.'''
CELL_COUNT = 3.0

'''
DEGREE_DIGITS is the number of significant digits to the right of the decimal
to use in latitude and longitude equality determination and representation. This 
should be set to 7 to preserve reversible transformations between coordinate systems 
down to a resolution of roughly 1 m.
''' 
DEGREE_DIGITS = 7

'''
SEMI_MAJOR_AXIS is the radius of the sphere at the equator for the WGS84 datum. 
Cell construction and lookup are based on projection of the geographic coordinates 
onto a sphere of this radius. As a result, edge lengths and areas of cells projected
onto the WGS84 ellipsoid will vary slightly by latitude.
'''
SEMI_MAJOR_AXIS = 6378137.0

'''
EDGE_LENGTH is the length of the edge of any face of the icosahedron inscribed on the 
sphere whose radius is SEMI_MAJOR_AXIS.
'''
EDGE_LENGTH = 6706370.116516389

'''
MID_EDGE_RADIUS is the distance from the center of the sphere to the midpoint of any
edge of any face of the icosahedron inscribed on the sphere whose radius is SEMI_MAJOR_AXIS.
'''
MID_EDGE_RADIUS = 5425567.394830056

'''
SURFACE_DISTANCE_BETWEEN_VERTEXES is the great circle distance between any two 
adjacent rhomboid vertexes on the sphere whose radius is SEMI_MAJOR_AXIS.
'''
SURFACE_DISTANCE_BETWEEN_VERTEXES = 7061546.20147
# SURFACE_DISTANCE_BETWEEN_VERTEXES = VERTEX_ANGLE*(2 * math.pi * SEMI_MAJOR_AXIS) / 360.0

'''
CELL_SIDE_LENGTH is the great circle distance between any two 
adjacent cell vertexes on the sphere whose radius is SEMI_MAJOR_AXIS.
'''
CELL_SIDE_LENGTH = SURFACE_DISTANCE_BETWEEN_VERTEXES / CELL_COUNT 

'''
RADIANS is constant representing the number by which to multiply a value in degrees
to get the equivalent value in radians.
'''
RADIANS = 0.017453292519943295

'''
PHI = (1+SQRT(5))/2 - the golden ratio. The vertices of an icosahedron of side length 
2 can be placed at the coordinates: 
(0, +-1, +-PHI)
(+-1, +-PHI, 0)
(+-PHI, 0, +-1)
'''
PHI = 1.618033988749895

'''
VERTEX_LAT is the latitude in degrees north or south of the equator for all
vertexes that are not at one of the poles.
'''
VERTEX_LAT = 26.565051177077997

'''
VERTEX_ANGLE is the angle in degrees from the center of the sphere between any two
adjacent vertexes.
''' 
VERTEX_ANGLE = 63.434948822922

#FACE_CENTER_LAT = 52.622631859350314
#CENTER_ANGLE = 41.810314895778575

# Precalculated useful sines and cosines.
COS_72 = 0.30901699437494734
SIN_72 = -0.9510565162951536
COS_VERTEX_LAT = 0.8944271909999159
SIN_VERTEX_LAT = -0.4472135954999581

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
    
KML = u'''<?xml version="1.0" encoding="UTF-8"?>                                                                                                          
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
</kml>'''

FORMAT = """.%sf"""

def flip(p):
    '''Swap latitude and longitude order in the tuple. KML requires longitude first.'''
    return (p[1], p[0])

def sqr(x):
    '''Square of x.'''
    return x * x

def signum(x):
    '''Sign of x.'''
    if x > 0:
        return 1.0
    if x < 0:
        return -1.0
    return 0

def truncate_lat_lng(p):
    '''Set the precision of a lat long in degrees to DEGREE_DIGITS.'''    
    return (truncate(p[0],DEGREE_DIGITS),truncate(p[1],DEGREE_DIGITS))

def truncate(x, digits):
    '''Set the representational precision of x to digits places to the right of the decimal.'''
    format_x = FORMAT % digits
    return format(x,format_x)

def equal_within_tolerance(a, b, digits):
    '''Determine if two floats are equal within the accuracy given by digits places to the right of the decimal.'''
    if math.fabs(a-b) <= math.pow(10,-digits):
        return True
    return False

def equal_lat_lng(p0, p1):
    '''Determine if two degree lat longs are the same within the tolerance given by digits places to the right of the decimal.'''
    return (equal_within_tolerance(p0[0], p1[0], DEGREE_DIGITS) and equal_within_tolerance(p0[1], p1[1], DEGREE_DIGITS)) 

def lng180(lng):
    '''Given a longitude in degrees, returns a longitude in degrees between {-180, 180].'''
    newlng=lng
    if lng <= -180:
        newlng = lng + 360
    elif lng > 180:
        newlng = lng - 360
    return float(truncate(newlng, DEGREE_DIGITS))

def lng360(lng):
    '''Given a longitude in degrees, returns a longitude in degrees between [0, 360}.'''
    if lng < 0:
        return float(truncate(lng + 360, DEGREE_DIGITS))
    if lng > 360:
        return float(truncate(lng - 360, DEGREE_DIGITS))
    return float(truncate(lng, DEGREE_DIGITS))

def xyz_from_lat_lng((lat, lng)):
    '''Returns the Cartesian coordinates of a lat long on a unit sphere.'''
    x = math.cos(lat * RADIANS) * math.cos(lng * RADIANS)
    y = math.cos(lat * RADIANS) * math.sin(lng * RADIANS)
    z = math.sin(lat * RADIANS)
    return (x, y, z)

def lat_lng_from_xyz((x, y, z)):
    '''Returns the lat long (in degrees) of Cartesian coordinates on a unit sphere.'''
    R=math.sqrt(x*x+y*y+z*z)
    znorm = z/R
    ynorm = y/R
    xnorm = x/R
    lng = math.atan2(ynorm,xnorm) / RADIANS
    lat = math.asin(znorm) / RADIANS
    return (lat, lng)

def great_circle_distance(start_lat_lng, end_lat_lng):
    '''
    Returns the distance along a great circle between two lat longs on the surface of a
    sphere of radius SEMI_MAJOR_AXIS.
    '''
    dLat = (end_lat_lng[0]-start_lat_lng[0]) * RADIANS
    dLon = (end_lat_lng[1]-start_lat_lng[1]) * RADIANS 
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(start_lat_lng[0] * RADIANS) * math.cos(end_lat_lng[0] * RADIANS) * math.sin(dLon/2) * math.sin(dLon/2) 
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return SEMI_MAJOR_AXIS * c 

def great_circle_intersection_xyz(p0, p1, p2, p3):
    '''
    Returns one of the two points on the sphere where two great circles, each defined by two Cartesian points,
    intersect. Result should be checked against expectation for hemisphere, and if not correct, take the 
    antipode.
    '''
    r = SEMI_MAJOR_AXIS
    x0 = p0[0]
    y0 = p0[1]
    z0 = p0[2]
    x1 = p1[0]
    y1 = p1[1]
    z1 = p1[2]
    x2 = p2[0]
    y2 = p2[1]
    z2 = p2[2]
    x3 = p3[0]
    y3 = p3[1]
    z3 = p3[2]
    a = (y0*z1 - y1*z0)
    b = -(x0*z1 - x1*z0)
    c = (x0*y1 - x1*y0)
    d = (y2*z3 - y3*z2)
    e = -(x2*z3 - x3*z2)
    f =  (x2*y3 - x3*y2)
    h = (d*c-f*a)/(e*a-d*b)
    g = ((-b*h - c)/a)
    k = math.sqrt(r*r/(g*g + h*h + 1))
    lat = math.asin(k/r)
    lng = math.atan2(h*k,g*k)
    return (float(truncate(lat/RADIANS, DEGREE_DIGITS)),float(truncate(lng/RADIANS, DEGREE_DIGITS)))

def great_circle_intersection_lat_lngs(p0, p1, p2, p3):
    '''
    Returns the lat long of the point within the same hemisphere where two great circles, defined by 
    two lat longs each, intersect.
    '''
    if p0[0] == -90.0:
        corrected_p0 = p1
        corrected_p1 = p0
    else:
        corrected_p0 = p0
        corrected_p1 = p1
    p = great_circle_intersection_xyz(xyz_from_lat_lng(corrected_p0), xyz_from_lat_lng(corrected_p1), xyz_from_lat_lng(p2), xyz_from_lat_lng(p3))
    corrected_lat = p[0]
    corrected_lng = p[1]
    # If the calculated intersection is in the opposite hemisphere, return the antipode instead 
    if not same_hemisphere(p,p0):
        corrected_lng = lng180(p[1] + 180)
        if p[0] != 90.0:
            corrected_lat = -1.0 * p[0] 
    return (corrected_lat,corrected_lng)

def same_hemisphere(p0, p1):    
    if great_circle_distance(p0,p1) > math.pi*SEMI_MAJOR_AXIS/2:
        return False
    return True
    
class Cell(object):
    '''
    Cell is the rhomboidal grid cell of a Triangular Mesh Grid over a sphere. The cell is 
    identified by the rhomboid within which it lies and the x and y offset indexes from the 
    southern vertex.
    ''' 
    
    @staticmethod
    def rotate(lat_lng, axis_lat_lng, rotation_angle):
        '''
        Return the lat long in degrees of the point determined by right-handed rotating an input
        point lat_lng through an angle rotation_angle about the axis axis_lat_lng.
        '''
        x = math.cos(lat_lng[0] * RADIANS) * math.cos(lat_lng[1] * RADIANS)
        y = math.cos(lat_lng[0] * RADIANS) * math.sin(lat_lng[1] * RADIANS)
        z = math.sin(lat_lng[0] * RADIANS)

        c1 = math.cos(axis_lat_lng[0] * RADIANS) * math.cos(axis_lat_lng[1] * RADIANS)
        c2 = math.cos(axis_lat_lng[0] * RADIANS) * math.sin(axis_lat_lng[1] * RADIANS)
        c3 = math.sin(axis_lat_lng[0] * RADIANS)
        cosa = math.cos(rotation_angle * RADIANS)
        sina = math.sin(rotation_angle * RADIANS)

        x1 = x * cosa
        y1 = y * cosa
        z1 = z * cosa

        x1 = x1 + ((1 - cosa) * (c1 * c1 * x + c1 * c2 * y + c1 * c3 * z))
        y1 = y1 + ((1 - cosa) * (c2 * c1 * x + c2 * c2 * y + c2 * c3 * z))
        z1 = z1 + ((1 - cosa) * (c3 * c1 * x + c3 * c2 * y + c3 * c3 * z))

        x1 = x1 + (c2 * z - c3 * y) * sina
        y1 = y1 + (c3 * x - c1 * z) * sina
        z1 = z1 + (c1 * y - c2 * x) * sina

        newlat = math.asin(z1) / RADIANS
        newlng = math.atan2(y1, x1) / RADIANS
        return (newlat, newlng)
           
    @staticmethod
    def polygon(rhomboid_num, x_index, y_index, cell_count=None):
        '''
        Return a list of points (lat longs in degrees) for the vertexes of a cell
        given by the rhomboid index number and x, y offset indexes from the southern vertex of the rhomboid.
        '''
        if cell_count == None:
            cell_count=CELL_COUNT
        s = Rhomboid.south_lat_lng(rhomboid_num)
        e = Rhomboid.east_lat_lng(rhomboid_num)
        w = Rhomboid.west_lat_lng(rhomboid_num)
        n = Rhomboid.north_lat_lng(rhomboid_num)
        out = 's=%s e=%s w=%s n=%s' % (s, e, w, n)
        logging.debug(out)

        bearing_we0 = Rhomboid.get_bearing(w,e)
        distance_we0 = CELL_SIDE_LENGTH*x_index
        distance_we1 = CELL_SIDE_LENGTH*(x_index+1)
        we0 = Rhomboid.get_point_from_distance_at_bearing(w, distance_we0, bearing_we0)
        we1 = Rhomboid.get_point_from_distance_at_bearing(w, distance_we1, bearing_we0)
        out = 'dist_we0=%s bearing_we0=%s we0=%s' % (distance_we0, bearing_we0, we0)
        logging.debug(out)
        out = 'dist_we1=%s bearing_we0=%s we1=%s' % (distance_we1, bearing_we0, we1)
        logging.debug(out)

        bearing_ew0 = Rhomboid.get_bearing(e,w)
        distance_ew0 = CELL_SIDE_LENGTH*y_index
        distance_ew1 = CELL_SIDE_LENGTH*(y_index+1)
        ew0 = Rhomboid.get_point_from_distance_at_bearing(e, distance_ew0, bearing_ew0)
        ew1 = Rhomboid.get_point_from_distance_at_bearing(e, distance_ew1, bearing_ew0)
        out = 'dist_ew0=%s bearing_ew0=%s ew0=%s' % (distance_ew0, bearing_ew0, ew0)
        logging.debug(out)
        out = 'dist_ew1=%s bearing_ew0=%s ew1=%s' % (distance_ew1, bearing_ew0, ew1)
        logging.debug(out)

        bearing_se0 = Rhomboid.get_bearing(s,e)
        distance_se0 = CELL_SIDE_LENGTH*x_index
        se0 = Rhomboid.get_point_from_distance_at_bearing(s, distance_se0, bearing_se0)
        out = 'dist_se0=%s bearing_se0=%s se0=%s' % (distance_se0, bearing_se0, se0)
        logging.debug(out)

        bearing_sw0 = Rhomboid.get_bearing(s,w)
        distance_sw0 = CELL_SIDE_LENGTH*y_index
        if s[0] == -90.0:
            sw0 = Rhomboid.get_point_from_distance_at_bearing((-90.0,w[1]), distance_sw0, bearing_sw0)
        else:
            sw0 = Rhomboid.get_point_from_distance_at_bearing(s, distance_sw0, bearing_sw0)
        out = 'dist_sw0=%s bearing_sw0=%s sw0=%s' % (distance_sw0, bearing_sw0, sw0)
        logging.debug(out)

        distance_se1 = CELL_SIDE_LENGTH*(x_index+1)
        se1 = Rhomboid.get_point_from_distance_at_bearing(s, distance_se1, bearing_se0)
        out = 'dist_se1=%s bearing_se0=%s se1=%s' % (distance_se1, bearing_se0, se1)
        logging.debug(out)
        
        distance_sw1 = CELL_SIDE_LENGTH*(y_index+1)
        if s[0] == -90.0:
            sw1 = Rhomboid.get_point_from_distance_at_bearing((-90.0,w[1]), distance_sw1, bearing_sw0)
        else:
            sw1 = Rhomboid.get_point_from_distance_at_bearing(s, distance_sw1, bearing_sw0)
        out = 'dist_sw1=%s bearing_sw0=%s sw1=%s' % (distance_sw1, bearing_sw0, sw1)
        logging.debug(out)
        
        bearing_en0 = Rhomboid.get_bearing(e,n)
        distance_en0 = CELL_SIDE_LENGTH*y_index
        en0 = Rhomboid.get_point_from_distance_at_bearing(e, distance_en0, bearing_en0)
        out = 'dist_en0=%s bearing_en0=%s en0=%s' % (distance_en0, bearing_en0, en0)
        logging.debug(out)

        distance_en1 = CELL_SIDE_LENGTH*(y_index+1)
        en1 = Rhomboid.get_point_from_distance_at_bearing(e, distance_en1, bearing_en0)
        out = 'dist_en1=%s bearing_en0=%s en1=%s' % (distance_en1, bearing_en0, en1)
        logging.debug(out)
        
        bearing_wn0 = Rhomboid.get_bearing(w,n)
        distance_wn0 = CELL_SIDE_LENGTH*x_index
        wn0 = Rhomboid.get_point_from_distance_at_bearing(w, distance_wn0, bearing_wn0)
        out = 'dist_wn0=%s bearing_wn0=%s wn0=%s' % (distance_wn0, bearing_wn0, wn0)
        logging.debug(out)

        distance_wn1 = CELL_SIDE_LENGTH*(x_index+1)
        pre_wn1 = Rhomboid.get_point_from_distance_at_bearing(w, distance_wn1, bearing_wn0)
        if pre_wn1[0] == 90.0:
            wn1 = (90.0, wn0[1])
        else:
            wn1 = pre_wn1
        out = 'dist_wn1=%s bearing_wn0=%s wn1=%s' % (distance_wn1, bearing_wn0, wn1)
        logging.debug(out)
                
        if x_index + y_index == 0 and rhomboid_num >= 5:
            # S vertex on S pole
            s_vertex = (-90.0,e[1])
            e_vertex = se1
            w_vertex = sw1
            out = '0a) S cell vertex on S pole e_vertex=%s w_vertex=%s' % (e_vertex, w_vertex)
            logging.debug(out)
            if x_index + y_index < CELL_COUNT-1:
                # N cell vertex in or on S triangle also
                n_vertex = great_circle_intersection_lat_lngs(se1,we1, sw1,ew1)
                out = '0b) N cell vertex in or on S triangle n_vertex=%s' % (str(n_vertex))
                logging.debug(out)
            else:
                # N cell vertex in N triangle
                if equal_lat_lng(en1,wn1):
                    n_vertex = en1
                else:
                    n_vertex = great_circle_intersection_lat_lngs(we1,wn1, ew1,en1)
                out = '0c) N cell vertex in N triangle n_vertex=%s' % (str(n_vertex))
                logging.debug(out)
                
        elif x_index + y_index <= CELL_COUNT:
            # S cell vertex in or on S triangle
            if equal_lat_lng(se0,sw0):
                s_vertex = se0
                out = '1a) S cell vertex in or on S triangle s_vertex=%s' % (str(s_vertex))
                logging.debug(out)
            else:
                # S cell vertex in N triangle
                s_vertex = great_circle_intersection_lat_lngs(se0,we0, sw0,ew0)
                out = '1b) S cell vertex in S triangle s_vertex=%s' % (str(s_vertex))
                logging.debug(out)
            if x_index + y_index < CELL_COUNT:
                # E, W cell vertexes in or on S triangle   
                if equal_lat_lng(se1,we1):
                    e_vertex = se1
                else:
                    e_vertex = great_circle_intersection_lat_lngs(se1,we1, sw0,ew0)
                if equal_lat_lng(sw1,we0):
                    w_vertex = sw1
                else:
                    w_vertex = great_circle_intersection_lat_lngs(se0,we0, sw1,ew1)
                out = '2) E,W cell vertexes in or on S triangle e_vertex=%s w_vertex=%s' % (e_vertex, w_vertex)
                logging.debug(out)
            else:
                e_vertex = great_circle_intersection_lat_lngs(ew0,en0, we1,wn1)
                w_vertex = great_circle_intersection_lat_lngs(we0,wn0, ew1,en1)
                out = '3) E,W cell vertexes in N triangle e_vertex=%s w_vertex=%s' % (e_vertex, w_vertex)
                logging.debug(out)
            if x_index + y_index < CELL_COUNT-1:
                # N cell vertex in or on S triangle also
                n_vertex = great_circle_intersection_lat_lngs(se1,we1, sw1,ew1)
                out = '4) N cell vertex in or on S triangle also n_vertex=%s' % (str(n_vertex))
                logging.debug(out)
            else:
                # N cell vertex in N triangle
                if equal_lat_lng(en1,wn1):
                    n_vertex = en1
                    out = '5a)'
                    logging.debug(out)
                else:
                    n_vertex = great_circle_intersection_lat_lngs(we1,wn1, ew1,en1)
                    out = '5d)'
                    logging.debug(out)
                out = '5) N cell vertex in N triangle n_vertex=%s' % (str(n_vertex))
                logging.debug(out)
        else:
            # Whole cell in N triangle
            e_vertex = great_circle_intersection_lat_lngs(ew0,en0, we1,wn1)
            w_vertex = great_circle_intersection_lat_lngs(we0,wn0, ew1,en1)
            n_vertex = great_circle_intersection_lat_lngs(we1,wn1, ew1,en1)
            s_vertex = great_circle_intersection_lat_lngs(ew0,en0, we0,wn0)
            out = '6) Whole cell in N triangle n_vertex=%s' % (str(n_vertex))
            logging.debug(out)

        out = 's_vertex=%s n_vertex=%s e_vertex=%s w_vertex=%s' % (s_vertex, n_vertex, e_vertex, w_vertex)
        logging.debug(out)
                 
        if s_vertex[0] == -90:
            p0 = (s_vertex[0],e_vertex[1])
            p1 = (s_vertex[0],w_vertex[1])
#            out = 'S half sv=-90 p0=%s p1=%s' % (p0, p1)
#            logging.debug(out)
            '''S-out, E, N, W, S-in'''
            return [flip(truncate_lat_lng(p0)), flip(truncate_lat_lng(e_vertex)), flip(truncate_lat_lng(n_vertex)), flip(truncate_lat_lng(w_vertex)), flip(truncate_lat_lng(p1))]
        elif n_vertex[0] == 90:
            p2 = (90.0,e_vertex[1])
            p3 = (90.0,w_vertex[1])
#            out = 'S half nv=+90: p2=%s p3=%s' % (p2, p3)
#            logging.debug(out)
            '''S, E, N-in, N-out, W, S'''
            return [flip(truncate_lat_lng(s_vertex)), flip(truncate_lat_lng(e_vertex)), flip(truncate_lat_lng(p2)), flip(truncate_lat_lng(p3)), flip(truncate_lat_lng(w_vertex)), flip(truncate_lat_lng(s_vertex))]
#        out = 'S half sv<>-90 nv<>90 s_vertex=%s' % (str(s_vertex))
#        logging.debug(out)
        '''S, E, N, W, S'''
        return [flip(truncate_lat_lng(s_vertex)), flip(truncate_lat_lng(e_vertex)), flip(truncate_lat_lng(n_vertex)), flip(truncate_lat_lng(w_vertex)), flip(truncate_lat_lng(s_vertex))]
        
    @staticmethod
    def createPlacemark(key, polygon):
        '''Render a KML placemark for a polygon.''' 
        data = (key, polygon, key)
        placemark = PLACEMARK_1 % data
        for c in polygon:
            point = '                        %s,%s,1\n' % (c[0], c[1])
            placemark = placemark + point
        placemark = placemark + PLACEMARK_2
        return placemark
    
    @staticmethod
    def createKmlMesh(cell_count):
        '''Render a triangular mesh of cells in KML.''' 
        placemarks = []
        for n in range(1):
            for x in range(cell_count):
                for y in range(cell_count):
#            for x in range(1):
#                for y in range(1):
                    polygon = Cell.polygon(n, x, y, cell_count)                                        
                    key = '%s-%s-%s' % (n, x, y)
                    p=Cell.createPlacemark(key, polygon)
                    placemarks.append(p)
                    out = 'n=%s x=%s y=%s' % (n, x, y)
                    logging.debug(out)
        return KML % ' '.join(placemarks)

    @staticmethod
    def get_canonical_south_point(x_index, y_index, cell_count=None):
        # This isn't correct. Rotations depend on the index.
        if not cell_count:
            cell_count = CELL_COUNT
        # Canonical rhomboid 0 has origin at -VERTEX_LAT, 0
        lat_lng0 = (-VERTEX_LAT, 0)
        axis_lat_lng = (VERTEX_LAT, -108)

        # Rotate NE from the origin by x_index cell widths
        lat_lng = Cell.rotate(lat_lng0, axis_lat_lng, 72 * x_index / cell_count)
            
        # Rotate NW from there by y_index cell widths
        lat_lng = Cell.rotate(lat_lng, (-VERTEX_LAT, -108), 72 * y_index / cell_count)
        return (lat_lng)

    @staticmethod
    def get_canonical_east_point(x_index, y_index, cell_count=None):
        if not cell_count:
            cell_count = CELL_COUNT
        # Start at canonical south point
        lat_lng0 = Cell.get_canonical_south_point(x_index, y_index)

        newlat = lat_lng0[0] + (90 + VERTEX_LAT) / (2 * cell_count)
        newlng = lat_lng0[1] + 36 / cell_count
        lat_lng = (newlat, newlng)
        return (lat_lng)

    @staticmethod
    def get_canonical_west_point(x_index, y_index, cell_count=None):
        if not cell_count:
            cell_count = CELL_COUNT
#        # Start at canonical south point
        lat_lng0 = Cell.get_canonical_south_point(x_index, y_index)
#        # Rotate NW from there by one cell width
        newlat = lat_lng0[0] + (90 + VERTEX_LAT) / (2 * cell_count)
        newlng = lat_lng0[1] - 36 / cell_count
        lat_lng = (newlat, newlng)
        return (lat_lng)

    @staticmethod
    def get_canonical_north_point(x_index, y_index, cell_count=None):
        if not cell_count:
            cell_count = CELL_COUNT

        # Start at canonical south point
        lat_lng0 = Cell.get_canonical_south_point(x_index, y_index)

        # Rotate N from there by one cell diagonal width
        newlat = lat_lng0[0] + (90 + VERTEX_LAT) / cell_count
        lat_lng = (newlat, lat_lng0[1])
        return (lat_lng)

    def __init__(self, rhomboid_num, x_index, y_index):
        self.rhomboid_num = rhomboid_num
        self.x_index = x_index
        self.y_index = y_index
    
class Face(object):
    @staticmethod
    def get_northern(lng):
        '''One of the get_northern five faces: 0=[-36,36}, 1=[36,108}, 2=[108,180},
            3=[180,-108}, 4=[-108,-36}.
        '''
        clng = lng180(lng)
        face = int(math.floor((clng + 360 + 36) / 72))
        if clng >= -36:
            face = face - 5
        return face

    @staticmethod
    def get_southern(lng):
        '''One of the get_southern five faces: 15=[0,72}, 16=[72,144}, 17=[144,-144},
        18=[-144,-72}, 19=[-72,0}.
        '''
        clng = lng180(lng)
        face = int(math.floor((clng + 360) / 72) + 15)
        if clng >= 0:
            face = face - 5
        return face

    @staticmethod
    def get_equatorial(lat, lng):
        ''' lat and lng are somewhere in the get_equatorial zone. Find which of ten
        equal sections around the equator lng falls into by rotating lng to
        equivalent longitude at lat=0.
        '''
        clng = lng180(lng)
        section = 0
        if clng >= -18:
            section = int(math.floor((18 + clng + (18 * (lat / VERTEX_LAT))) / 36))
        else:
            section = int(math.floor((18 + 360 + clng + (18 * (lat / VERTEX_LAT))) / 36))
        # Even numbered sections are northern facets [5,9]:
        if section % 2 == 0:
            face = 5 + (section / 2)
        # Odd numbered sections are get_southern facets [10,14]:
        else:
            face = 10 + ((section - 1) / 2)
        return face

    @staticmethod
    def get_face_number(lat, lng):
        '''The face is one of the sides of an icosahedron get_face_number consecutively
        from 0 beginning with the get_northern face having its center at lng=0.
        Return the face of the icosahedron onto which the lat lng projects.
        '''
        if lat >= VERTEX_LAT:
            face = Face.get_northern(lng)
        elif lat < -(VERTEX_LAT):
            face = Face.get_southern(lng)
        else:
            face = Face.get_equatorial(lat, lng)
        return face

    def __init__(self, lat, lng):
        self.lat = lat
        self.lng = lng
        self.num = self.get_face_number(lat, lng)

class Rhomboid(object):
    @staticmethod
    def canonical_lat(face):
        if face.num < 10:
            return face.lat
        clng = Rhomboid.canonical_lng(face)
        x = math.cos(face.lat * RADIANS) * math.cos(clng * RADIANS)
        y = math.cos(face.lat * RADIANS) * math.sin(clng * RADIANS)
        z = math.sin(face.lat * RADIANS)
        c1 = COS_VERTEX_LAT * COS_72
        c2 = COS_VERTEX_LAT * SIN_72
        c3 = SIN_VERTEX_LAT
        cosa = COS_72
        sina = -(SIN_72)
        z1 = z * cosa
        z1 = z1 + ((1 - cosa) * (c3 * c1 * x + c3 * c2 * y + c3 * c3 * z))
        z1 = z1 + (((c1 * y) - (c2 * x)) * sina)
        newlat = math.asin(z1) / RADIANS
#        out = 'canonical_lat(face.num=%s) x=%s, y=%s, z=%s z1=%s face.lat=%s clng=%s' % (face.num, x, y, z, z1, face.lat, clng)
#        logging.debug(out)
        return newlat
    
    @staticmethod
    def get_bearing(start_lat_lng, end_lat_lng):
        if start_lat_lng[0] == -90.0:
            return 0
        y = math.sin((end_lat_lng[1] - start_lat_lng[1]) * RADIANS) * math.cos(end_lat_lng[0] * RADIANS)
        x = math.cos(start_lat_lng[0] * RADIANS) * math.sin(end_lat_lng[0] * RADIANS) - math.sin(start_lat_lng[0] * RADIANS) * math.cos(end_lat_lng[0] * RADIANS) * math.cos((end_lat_lng[1] - start_lat_lng[1]) * RADIANS)
        b = math.atan2(y, x) / RADIANS
        return float(truncate(b, DEGREE_DIGITS))
    
    @staticmethod
    def get_point_from_distance_at_bearing(start_lat_lng, distance, bearing):
        ad = distance/SEMI_MAJOR_AXIS
        lat1 = start_lat_lng[0] * RADIANS
        lng1 = start_lat_lng[1] * RADIANS
        b = bearing * RADIANS
        lat2 = math.asin( math.sin(lat1) * math.cos(ad) + \
                          math.cos(lat1)*math.sin(ad)*math.cos(b) )
        lng2 = lng1 + math.atan2( math.sin(b)*math.sin(ad)*math.cos(lat1), \
                                  math.cos(ad)-math.sin(lat1)*math.sin(lat2))
        return (float(truncate(lat2/RADIANS,DEGREE_DIGITS)), float(truncate(lng2/RADIANS,DEGREE_DIGITS)))
                                 
    @staticmethod
    def rotation_axis_xyz(start_lat_lng, end_lat_lng):
        
#        In 2- or 3-dimensional Euclidean space, two vectors are orthogonal if their dot product is zero.
        
        start_xyz = xyz_from_lat_lng(start_lat_lng)
        end_xyz = xyz_from_lat_lng(end_lat_lng)

#        out = 'ROTATION_AXIS_XYZ: start_xyz=%s end_xyz=%s' % (start_xyz, end_xyz)
#        logging.debug(out)

        a1 = start_xyz[0]
        a2 = start_xyz[1]
        a3 = start_xyz[2]
        b1 = end_xyz[0]
        b2 = end_xyz[1]
        b3 = end_xyz[2]
    
        c1 = a2*b3 - a3*b2
        c2 = a3*b1 - a1*b3
        c3 = a1*b2 - a2*b1
        
        # test orthogonality by dot products
#        R = math.sqrt(c1*c1 + c2*c2 + c3*c3)
#        a.c = a1c1+a2c2+a3c3 = 0
#        b.c = b1c1+b2c2+b3c3 = 0
#        c1(a1-b1) + c2(a2-b2) + c3(a3-b3) = 0

#        c1 = c1/R
#        c2 = c2/R
#        c3 = c3/R
#        adotc = a1*c1 + a2*c2 + a3*c3
#        bdotc = b1*c1 + b2*c2 + b3*c3
#
#        out = ' adotc=%s bdotc=%s' % (adotc, bdotc)
#        logging.debug(out)

        return (c1, c2, c3)

    @staticmethod
    def get_bearing_ne(rhomboid_num):
        start_lat_lng = Rhomboid.south_lat_lng(rhomboid_num)
        end_lat_lng = Rhomboid.east_lat_lng(rhomboid_num)
        return Rhomboid.get_bearing(start_lat_lng, end_lat_lng, rhomboid_num)

    @staticmethod
    def get_bearing_nw(rhomboid_num):
        start_lat_lng = Rhomboid.south_lat_lng(rhomboid_num)
        end_lat_lng = Rhomboid.west_lat_lng(rhomboid_num)
        return Rhomboid.get_bearing(start_lat_lng, end_lat_lng, rhomboid_num)

    @staticmethod
    def south_lat_lng(rhomboid_num):
        lat = -VERTEX_LAT
        lng = 0
        if rhomboid_num < 5:
            # Latitude is correct already.
            lng = lng180(72 * (rhomboid_num % 5))
        else:
            lat = -90
            lng = lng180(72 + 72 * (rhomboid_num % 5))
        # Longitude might as well be zero for the pole.
        # But note that polygon construction from vertices at the poles may 
        # be sensitive to the longitude when connecting to an adjacent vertex.
        return (lat, lng)

    @staticmethod
    def east_lat_lng(rhomboid_num):
        lat = VERTEX_LAT
        lng = 0
        if rhomboid_num < 5:
            # Latitude is correct already.
            lng = lng180(36 + 72 * (rhomboid_num % 5))
        else:
            lat = -VERTEX_LAT
            lng = lng180(72 + 72 * (rhomboid_num % 5))
        return (lat, lng)

    @staticmethod
    def west_lat_lng(rhomboid_num):
        lat = VERTEX_LAT
        lng = 0
        if rhomboid_num < 5:
            # Latitude is correct already
            lng = lng180(-36 + 72 * (rhomboid_num % 5))
        else:
            lat = -VERTEX_LAT
            lng = lng180(72 * (rhomboid_num % 5))
        return (lat, lng)

    @staticmethod
    def north_lat_lng(rhomboid_num):
        lat = 90
        lng = 0
        # For rhomboid_num < 5:
        # Latitude is correct already.
        # Longitude might as well be zero for the pole.
        # But note that polygon construction from vertices at the poles may 
        # be sensitive to the longitude when connecting to an adjacent vertex.
        if rhomboid_num >= 5:
            lat = VERTEX_LAT
            lng = lng180(36 + 72 * (rhomboid_num % 5))
        return (lat, lng)

    @staticmethod
    def south_xyz(rhomboid_num):
        return xyz_from_lat_lng(Rhomboid.south_lat_lng(rhomboid_num))

    @staticmethod
    def east_xyz(rhomboid_num):
        return xyz_from_lat_lng(Rhomboid.east_lat_lng(rhomboid_num))

    @staticmethod
    def west_xyz(rhomboid_num):
        return xyz_from_lat_lng(Rhomboid.west_lat_lng(rhomboid_num))

    @staticmethod
    def north_xyz(rhomboid_num):
        return xyz_from_lat_lng(Rhomboid.north_lat_lng(rhomboid_num))

    @staticmethod
    def canonical_lng(face):
        clng = lng360(face.lng)
        if face.num >= 10 and face.num < 15:
            clng = clng - 36 - (face.num % 5) * 72;
        else:
            clng = clng - (face.num % 5) * 72;
#        out = 'canonical_lng(face.num=%s), clng=%s' % (face.num, clng)
#        logging.debug(out)
        return lng360(clng)
    
    @staticmethod
    def cartesian_dist(fromlat, fromlng, tolat, tolng, radius):
        from_x = radius * math.cos(fromlat * RADIANS) * math.cos(fromlng * RADIANS)
        from_y = radius * math.cos(fromlat * RADIANS) * math.sin(fromlng * RADIANS)
        from_z = radius * math.sin(fromlat * RADIANS)
        to_x = radius * math.cos(tolat * RADIANS) * math.cos(tolng * RADIANS)
        to_y = radius * math.cos(tolat * RADIANS) * math.sin(tolng * RADIANS)
        to_z = radius * math.sin(tolat * RADIANS)
        dist = sqr(from_x - to_x)
        dist = dist + sqr(from_y - to_y)
        dist = dist + sqr(from_z - to_z)
        dist = math.sqrt(dist)
#        out = 'fromlat = %s fromlng = %s \ntolat = %s tolng = %s' % (fromlat, fromlng, tolat, tolng)
#        logging.debug(out)
#        out = 'dist = %s \n  from_x = %s, from_y = %s from_z = %s \n  to_x = %s to_y = %s to_z = %s' % (dist, from_x, from_y, from_z, to_x, to_y, to_z)
#        logging.debug(out)
        return dist
    
    @staticmethod
    def calc_num(face):
        if face.num < 10:
            facet = face.num % 5
        else:
            facet = (face.num % 5) + 5
        return facet

    @staticmethod
    def get_x_dist(clat, clng):
        cface = Face.get_face_number(clat, clng)
#        out = 'get_x_dist(clat = %s clng = %s) VERTEX_LAT = %s diff = %s' % (clat, clng, VERTEX_LAT, VERTEX_LAT-clat)
#        logging.debug(out)
        if cface == 5:
            d3 = Rhomboid.cartesian_dist(clat, clng, VERTEX_LAT, 36, SEMI_MAJOR_AXIS)
            d0 = Rhomboid.cartesian_dist(clat, clng, -(VERTEX_LAT), 0, SEMI_MAJOR_AXIS)
        else:
            d0 = Rhomboid.cartesian_dist(clat, clng, VERTEX_LAT, -36, SEMI_MAJOR_AXIS)
            d3 = Rhomboid.cartesian_dist(clat, clng, 90, 0, SEMI_MAJOR_AXIS)
        return (d0, d3)

    @staticmethod
    def get_y_dist(clat, clng):
        cface = Face.get_face_number(clat, clng)
        if cface == 5:
            d0 = Rhomboid.cartesian_dist(clat, clng, -(VERTEX_LAT), 0, SEMI_MAJOR_AXIS)
            d0 = math.floor(d0)
            d1 = Rhomboid.cartesian_dist(clat, clng, VERTEX_LAT, -36, SEMI_MAJOR_AXIS)
            d1 = math.floor(d1)
        else:
            d0 = Rhomboid.cartesian_dist(clat, clng, VERTEX_LAT, 36, SEMI_MAJOR_AXIS)
            d0 = math.floor(d0)
            d1 = Rhomboid.cartesian_dist(clat, clng, 90, 0, SEMI_MAJOR_AXIS)
            d1 = math.floor(d1)
        return (d0, d1)

    @staticmethod
    def get_x_edge_fraction(d0, d3):
        dist_from_d0 = (EDGE_LENGTH * EDGE_LENGTH - d3 * d3 + d0 * d0) / (2 * EDGE_LENGTH)
        h = math.sqrt(d0 * d0 - dist_from_d0 * dist_from_d0)
        dp = h * math.tan(30 * RADIANS)        
        dist_from_d0 = dist_from_d0 + dp
        if dist_from_d0 > (d0 / 2):
            edge_angle = VERTEX_ANGLE / 2
            edge_angle = edge_angle + math.atan((dist_from_d0 - (d0 / 2)) / MID_EDGE_RADIUS) / RADIANS
        else:
            edge_angle = VERTEX_ANGLE / 2
            edge_angle = edge_angle - math.atan(((d0 / 2) - dist_from_d0) / MID_EDGE_RADIUS) / RADIANS
        edge_fraction = edge_angle / VERTEX_ANGLE
        return edge_fraction

    @staticmethod
    def get_y_edge_fraction(d0, d1):
        dist_from_d0 = (EDGE_LENGTH * EDGE_LENGTH - d1 * d1 + d0 * d0) / (EDGE_LENGTH * 2)
        h = math.sqrt(d0 * d0 - dist_from_d0 * dist_from_d0)        
        dp = h * math.tan(30 * RADIANS) 
        dist_from_d0 = dist_from_d0 + dp       
        dist_from_d0 = math.floor(dist_from_d0)
        if dist_from_d0 > (d0 / 2):
            edge_angle = VERTEX_ANGLE / 2
            edge_angle = edge_angle + math.atan((dist_from_d0 - (d0 / 2)) / MID_EDGE_RADIUS) / RADIANS
        else:
            edge_angle = VERTEX_ANGLE / 2
            edge_angle = edge_angle - math.atan(((d0 / 2) - dist_from_d0) / MID_EDGE_RADIUS) / RADIANS
        edge_fraction = edge_angle / VERTEX_ANGLE
        return edge_fraction
        
    @staticmethod
    def calc_x_index(clat, clng, cell_count = None):
        if cell_count == None:
            cell_count = CELL_COUNT
            
        d0, d3 = Rhomboid.get_x_dist(clat, clng)
#        out = 'd0=%s, d3=%s' % (d0, d3)
#        logging.debug(out)
        if d3 < 1:
            return cell_count - 1;
        if d0 < 1:
            return 0;        
        edge_fraction = Rhomboid.get_x_edge_fraction(d0, d3)
        cell_index = int(cell_count * edge_fraction)
        return cell_index - 1
        
    @staticmethod
    def calc_y_index(clat, clng, cell_count = None):
        if cell_count == None:
            cell_count = CELL_COUNT

        d0, d1 = Rhomboid.get_y_dist(clat, clng)
#        out = 'd0=%s, d1=%s' % (d0, d1)
#        logging.debug(out)
        if d0 < 1:
            return 0
        if d1 < 1:
            return cell_count - 1
        edge_fraction = Rhomboid.get_y_edge_fraction(d0, d1)
        cell_index = int(cell_count * edge_fraction)
        return cell_index
        
    def __init__(self, lat, lng, cell_count = None):
        self.face = Face(lat, lng)
        self.clat = self.canonical_lat(self.face)
        self.clng = self.canonical_lng(self.face)
        self.num = self.calc_num(self.face)
        self.x = self.calc_x_index(self.clat, self.clng, cell_count)
        self.y = self.calc_y_index(self.clat, self.clng, cell_count)
        self.key = '%s-%s-%s' % (self.num, self.x, self.y)

def get_cell_polygon(cell_key, cell_count = None):
    rhomboid_num, x, y = cell_key.split('-')
    return Cell.polygon(rhomboid_num, x, y, cell_count)

def get_cell_key(lat, lng, cell_count = None):
    return Rhomboid(lat,lng, cell_count).key

def get_tile(from_ll, to_ll, cell_count = None):
    '''
    Return a list of TMG keys given a direction-sensitive geographic coordinate bounding box. 
    Orientation of the bounding box matters, as this function is meant to support bounding 
    boxes that span more than 180 degrees of longitude. The tile will be defined by the 
    direction proceeding from the corner given by from_ll to the corner given by to_ll. Then
    calculations will be done from NW to SE corner of bounding box.
    '''
    if from_ll[0] >= to_ll[0]:
        bb_n = from_ll[0]
        bb_s = to_ll[0]
    else:
        bb_n = to_ll[0]
        bb_s = from_ll[0]
    from_lng = lng180(from_ll[1])
    to_lng = lng180(to_ll[1])
    from_sign = signum(from_lng)
    to_sign = signum(to_lng)
    if from_sign*to_sign >= 0:
        '''
        from and to points are in the same E or W hemisphere.
        positive lng_dist means go E from from_ll to to_ll.
        '''
        lng_dist = from_lng - to_lng
        if lng_dist >= 0:
            bb_w = from_lng
            bb_e = to_lng
        else:
            bb_w = to_lng
            bb_e = from_lng
    else:
        '''
        from and to points are in opposite E and W hemispheres.
        positive lng_dist means go E from from_ll to to_ll.
        '''
        if to_lng < 0:
            lng_dist = -180 - to_lng + from_lng - 180 
        else:
            lng_dist = -180 - from_lng + to_lng - 180 

        if lng_dist >= 0:
            bb_w = from_lng
            bb_e = to_lng
        else:
            bb_w = to_lng
            bb_e = from_lng
            
    '''Bounding box defined by vertexes (bb_n, bb_w) and (bb_s, bb_e).'''
    bb = ((bb_n, bb_w), (bb_s, bb_e))

    if cell_count == None:
        cell_count = CELL_COUNT

    start_cell_key = get_cell_key(bb_n, bb_w, cell_count)
    '''Tile will be stored as a list of keys.'''
    key_list = []
    while start_cell_key != None and cell_in_bb(start_cell_key, bb, cell_count):
        cell_key_sw = next_cell_sw(start_cell_key, cell_count)
        if cell_key_sw != None and cell_in_bb(cell_key_sw, bb, cell_count):
            start_cell_key = cell_key_sw
        cursor = start_cell_key
        while cell_in_bb(cursor, bb, cell_count):
            key_list.append(cursor)
            cursor = next_cell_ne(cursor, cell_count)
        last_cell_key = start_cell_key
        start_cell_key = next_cell_south(bb_w, start_cell_key, cell_count)
    next_key_ne = next_cell_ne(last_cell_key, cell_count)
    if next_cell_east(last_cell_key, cell_count) == next_key_ne:
        start_cell_key = next_cell_east(next_key_ne, cell_count)
    else:
        start_cell_key = next_cell_east(last_cell_key, cell_count)
    cursor = start_cell_key
    while cell_in_bb(cursor, bb, cell_count):
        key_list.append(cursor)
        cursor = next_cell_ne(cursor, cell_count)
    start_cell_key = next_cell_east(bb_s, start_cell_key,cell_count)
    return key_list

def cell_in_bb(cell_key, bb, cell_count = None):
    '''
    Return True if any part of the cell referred to by cell_key lies within the NW to SE-oriented 
    bounding box bb, False otherwise.
    '''

    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT

    bb_n = bb[0][0]
    bb_w = bb[0][1]
    bb_s = bb[1][0]
    bb_e = bb[1][1]

    p = get_cell_polygon(cell_key, cell_count)
    '''
    For South vertex at -90 latitude:
        S-out, E, N, W, S-in
    For North vertex at 90 latitude:
        S, E, N-in, N-out, W, S
    For all other cells:
        S, E, N, W, S
    '''
    s = p[0][0]
    e = p[1][1]
    n = p[2][0]
    w = p[len(p)-2][1]
    
    if s > bb_n:
        '''Southern vertex of cell is north of northern limit of bounding box'''
        return False
    if n < bb_s:
        '''Northern vertex of cell is south of southern limit of bounding box'''
        return False

    '''
    e_w_dist is the distance (in degrees) between the east and
    west sides of the bounding box.
    '''
    e_w_dist = bb_e - bb_w
    if e_w_dist < 0:
        '''
        Bounding box crosses longitude 180. Translate all longitudes east by the 
        offset value so that the simple test for in bounding box can be done.
        '''
        # e_w_dist = e_w_dist + 360
        offset = bb_w + 180
        e = e + offset
        w = w + offset
        bb_w = bb_w + offset
        bb_e = bb_e + offset
    if e > bb_w and w < bb_e:
        '''
        East vertex of the cell is east of the west side of the bounding box and
        west vertex of the cell is west of the east side of the bounding box.
        ''' 
        return True
    return False
    
def next_cell_ne(cell_key, cell_count = None):        
    '''
    Return the cell (if any) adjacent to the cell given by rhomboid_num, x, y, cell_count
    in the NE direction.
    ''' 

    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT

    rhomboid_num, x, y = cell_key.split('-')

    if x < cell_count - 1:
        '''Next cell NE is in the same rhomboid.'''
        return (rhomboid_num, x + 1, y)

    '''Next cell NE is in an adjacent rhomboid.'''
    if rhomboid_num < 5:
        '''The rhomboid is one of the 5 northern ones.'''
        if x == cell_count - 1 and y == cell_count - 1:
            '''The given cell is the northernmost cell in the given rhomboid.'''
            return None
        return ((rhomboid_num + 1) % 5, 0, y)
    
    '''The rhomboid is one of the 5 southern ones.'''
    return ((rhomboid_num + 1) % 5, 0, y)

def next_cell_se(cell_key, cell_count = None):
    '''
    Return the cell (if any) adjacent to the cell given by rhomboid_num, x, y, cell_count
    in the SE direction.
    ''' 
    
    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT
    
    rhomboid_num, x, y = cell_key.split('-')

    if y > 0:
        '''Next cell SE is in the same rhomboid.'''
        return (rhomboid_num, x, y - 1)

    '''Next cell SE is in an adjacent rhomboid.'''
    if rhomboid_num < 5:
        '''The rhomboid is one of the 5 northern ones.'''
        return (5 + (rhomboid_num % 5), x, cell_count - 1)

    '''The rhomboid is one of the 5 southern ones.'''
    if x == 0 and y == 0:
        '''The given cell is the southernmost cell in the given rhomboid.'''
        return None
            
    return (5 + (rhomboid_num + 1) % 5, x, cell_count - 1)

def next_cell_sw(cell_key, cell_count = None):
    '''
    Return the cell (if any) adjacent to the cell given by rhomboid_num, x, y, cell_count
    in the SW direction.
    ''' 
    
    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT
    
    rhomboid_num, x, y = cell_key.split('-')

    if x > 0:
        '''Next cell SW is in the same rhomboid.'''
        return (rhomboid_num, x - 1, y)

    '''Next cell SW is in an adjacent rhomboid.'''
    if rhomboid_num < 5:
        '''The rhomboid is one of the 5 northern ones.'''
        return (5 + (rhomboid_num + 4) % 5, cell_count - 1, y)

    '''The rhomboid is one of the 5 southern ones.'''
    if x == 0 and y == 0:
        '''The given cell is the southernmost cell in the given rhomboid.'''
        return None
            
    return (5 + (rhomboid_num - 1) % 5, cell_count - 1, y)

def next_cell_south(lng, cell_key, cell_count = None):
    '''
    Return the cell (if any) due south of the given longitude and
    cell given by cell_key, cell_count.
    ''' 
    
    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT

    next_se = next_cell_se(cell_key, cell_count)
    if next_se == None:
        '''The given cell is the southernmost cell in a southern rhomboid.'''
        return None
    if lng_in_cell(lng, next_se, cell_count):
        '''lng crosses the cell to the SE of the given cell.'''
        return next_se
    '''lng crosses the cell to the SW of the given cell.'''
    return next_cell_sw(cell_key, cell_count)

def lng_in_cell(lng, cell_key, cell_count = None):

    p = get_cell_polygon(cell_key, cell_count)
    '''
    For South vertex at -90 latitude:
        S-out, E, N, W, S-in
    For North vertex at 90 latitude:
        S, E, N-in, N-out, W, S
    For all other cells:
        S, E, N, W, S
    So, the east vertex is always the second point in the cell polygon and
    the west vertex is always the penultimate point in the cell polygon. 
    '''
    e = p[1][1]
    w = p[len(p)-2][1]
    
    '''
    e_w_dist is the distance (in degrees) between the east and
    west vertexes of the given cell.
    '''
    e_w_dist = e - w
    if e_w_dist < 0:
        '''Cell crosses longitude 180.'''
        e_w_dist = e_w_dist + 360
    '''
    lng_w_dist is the distance (in degrees) between the given lng and the
    west vertex of the given cell.
    '''
    lng_w_dist = lng - w
    if lng_w_dist < 0:
        '''West vertex and lng in opposite hemispheres.'''
        lng_w_dist = lng_w_dist + 360
    if lng_w_dist > e_w_dist:
        '''lng is not between east and west vertexes'''
        return False
    
    '''lng is between east and west vertexes'''
    return True

def next_cell_east(lat, cell_key, cell_count = None):
    '''
    Return the cell (if any) due east of the given longitude and
    cell given by cell_key, cell_count.
    ''' 
    
    '''If no cell_count is provided, use the default constant CELL_COUNT.'''
    if cell_count == None:
        cell_count = CELL_COUNT

    next_se = next_cell_se(cell_key, cell_count)
    if next_se == None:
        '''
        The given cell is the southernmost cell in a southern rhomboid.
        next_se is east of cell_key.
        '''
        rhomboid_num, x, y = cell_key.split('-')
        next_se = '%s-%s-%s' % (5 + (rhomboid_num + 1) % 5, x, y)
        return next_se
    if lat_in_cell(lat, next_se, cell_count):
        '''lat crosses the cell to the SE of the given cell.'''
        return next_se
    '''lat crosses the cell to the NE of the given cell.'''
    return next_cell_ne(cell_key, cell_count)

def lat_in_cell(lat, cell_key, cell_count = None):
    p = get_cell_polygon(cell_key, cell_count)
    '''
    For South vertex at -90 latitude:
        S-out, E, N, W, S-in
    For North vertex at 90 latitude:
        S, E, N-in, N-out, W, S
    For all other cells:
        S, E, N, W, S
    '''
    s = p[0][0]
    n = p[2][0]
    
    if lat < s:
        '''lat is south of southern vertex of cell'''
        return False
    if lat > n:
        '''lat is north of northern vertex of cell'''
        return False
    return True

#def intersect_east(p, rhomboid_num, x, y, cell_count=None):
#    if cell_count == None:
#        cell_count = CELL_COUNT
#    lat = p[0]
#    lng = p[1]    
#    polygon = Cell.polygon(rhomboid_num, x, y)
#    s_vertex = polygon[0]
#    e_vertex = polygon[1]
#    n_vertex = polygon[2]
#    if lat == e_vertex[0]:
#        return (lat, e_vertex[1])
#    if lat > e_vertex[0]:
#        '''
#        lat is north of e_vertex, find lat intersect with great circle between
#        e_vertex and n_vertex.
#        '''
#        lng_intersection = nearer_parallel_great_circle_intersection(p, e_vertex, n_vertex)
#    else:
#        '''
#        lat is south of e_vertex, find lat intersect with great circle betweem
#        e_vertex and s_vertex.
#        '''
#        lng_intersection = nearer_parallel_great_circle_intersection(p, s_vertex, e_vertex)
#    return lng_intersection
    
#def nearer_parallel_great_circle_intersection(p0, p1, p2):
#    '''
#    Return the nearer (if any) of two intersections where a parallel given by point p0[0] intersect
#    a great circle defined by p1 and p2.
#    '''
#    lat1 = p1[0] * RADIANS
#    lat2 = p2[0] * RADIANS
#    lat3 = p0[0] * RADIANS
#    lon1 = p1[1] * RADIANS
#    lon2 = p2[1] * RADIANS
#    l12 = lon1 - lon2
#    A = math.sin(lat1) * math.cos(lat2) * math.cos(lat3) * math.sin(l12)
#    B = math.sin(lat1) * math.cos(lat2) * math.cos(lat3) * math.cos(l12) - math.cos(lat1) * math.sin(lat2) * math.cos(lat3)
#    C = math.cos(lat1) * math.cos(lat2) * math.sin(lat3) * math.sin(l12)
#    lon = math.atan2(B, A) # atan2(y,x) convention
#    if (math.abs(C) > math.sqrt(A * A + B * B)):
#        "no crossing"
#        return
#    else:
#        dlon = math.acos(C / math.sqrt(A * A + B * B))
#        lon3_1 = ((lon1 + dlon + lon + math.pi) % 2 * math.pi) - math.pi / RADIANS
#        lon3_2 = ((lon1 - dlon + lon + math.pi) % 2 * math.pi) - math.pi / RADIANS
##        lon3_1 = math.mod(lon1 + dlon+lon + math.pi, 2*math.pi) - math.pi
##        lon3_2 = math.mod(lon1 - dlon+lon + math.pi, 2*math.pi) - math.pi
#    if great_circle_distance(p0, (p0[0], lon3_1)) < great_circle_distance(p0, (p0[0], lon3_2)):
#        return (p0[0], lon3_1)
#    return (p0[0], lon3_2) 

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)    
    
    f = open('python.out', 'w')

    out = '%s' % (Cell.createKmlMesh(CELL_COUNT))
    logging.debug(out)

    f.flush()
    f.close()