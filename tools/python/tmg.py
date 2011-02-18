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

VERTEX_LAT = 26.565051177077997
CELL_COUNT = 4
PHI = 1.618033988749895
RADIANS = 0.017453292519943295
SEMI_MAJOR_AXIS = 6378137.0
VERTEX_ANGLE = 63.434948822922
EDGE_LENGTH = 6706370.116516389
MID_EDGE_RADIUS = 5425567.394830056
FACE_CENTER_LAT = 52.622631859350314
CENTER_ANGLE = 41.810314895778575
VERTEX_LAT = 26.565051177077997
COS_72 = 0.30901699437494734
SIN_72 = -0.9510565162951536
COS_VERTEX_LAT = 0.8944271909999159
SIN_VERTEX_LAT = -0.4472135954999581

def lng180(lng):
    '''Returns a longitude in {-180, 180].'''
    if lng <= -180:
        return lng + 360
    if lng > 180:
        return lng - 360
    return lng
   
def lng360(lng):
    '''Returns a longitude in [0, 360}.'''
    if lng < 0:
        return lng + 360
    if lng > 360:
        return lng - 360
    return lng

def sqr(x):
    return x * x

class Cell(object):
    @staticmethod
    def rotate(lat_lng, axis_lat_lng, rotation_angle):
        x = math.cos(lat_lng[0] * RADIANS) * math.cos(lat_lng[1] * RADIANS)
        y = math.cos(lat_lng[0] * RADIANS) * math.sin(lat_lng[1] * RADIANS)
        z = math.sin(lat_lng[0] * RADIANS)
        c1 = math.cos(axis_lat_lng[0] * RADIANS) * math.cos(axis_lat_lng[1] * RADIANS)
        c2 = math.cos(axis_lat_lng[0] * RADIANS) * math.sin(axis_lat_lng[1] * RADIANS)
        c3 = math.sin(axis_lat_lng[0] * RADIANS)
        cosa = math.cos(rotation_angle * RADIANS)
        sina = math.sin(rotation_angle * RADIANS)
        x1 = x * cosa
        x1 = x1 + ((1 - cosa) * (c1 * c1 * x + c1 * c2 * y + c1 * c3 * z))
        x1 = x1 + (c2 * z - c3 * y) * sina
        y1 = y * cosa
        y1 = y1 + ((1 - cosa) * (c2 * c1 * x + c2 * c2 * y + c2 * c3 * z))
        y1 = y1 + (c3 * x - c1 * z) * sina
        z1 = z * cosa
        z1 = z1 + ((1 - cosa) * (c3 * c1 * x + c3 * c2 * y + c3 * c3 * z))
        z1 = z1 + (c1 * y - c2 * x) * sina
        newlat = math.asin(z1) / RADIANS
        newlng = math.atan2(y1, x1)
        return (newlat, newlng)

    @staticmethod
    def polygon(rhomboid_num, x_index, y_index):
        return [Cell.get_point_from_canonical(rhomboid_num, Cell.get_canonical_south_point(x_index, y_index)),
                Cell.get_point_from_canonical(rhomboid_num, Cell.get_canonical_east_point(x_index, y_index)),
                Cell.get_point_from_canonical(rhomboid_num, Cell.get_canonical_north_point(x_index, y_index)),
                Cell.get_point_from_canonical(rhomboid_num, Cell.get_canonical_west_point(x_index, y_index)),
                Cell.get_point_from_canonical(rhomboid_num, Cell.get_canonical_south_point(x_index, y_index))
                ]

    @staticmethod
    def get_canonical_south_point(x_index, y_index):
        # Canonical rhomboid 0 has origin at -VERTEX_LAT, 0
        lat_lng0 = (-VERTEX_LAT, 0)
        axis_lat_lng = (VERTEX_ANGLE, -108)
        lat0 = -VERTEX_LAT
        lng0 = 0
        # Rotate NE from the origin by x_index cell widths
        lat_lng = Cell.rotate(lat_lng0, axis_lat_lng, VERTEX_ANGLE * x_index / CELL_COUNT)
            
        # Rotate NW from there by y_index cell widths
        lat_lng = Cell.rotate(lat_lng, (-VERTEX_ANGLE, -108), VERTEX_ANGLE * y_index / CELL_COUNT)

        return (lat_lng)

    @staticmethod
    def get_canonical_east_point(x_index, y_index):
        # Start at canonical south point
        lat_lng0 = Cell.get_canonical_south_point(x_index, y_index)
        # Rotate NE from there by one cell width
        lat_lng = Cell.rotate(lat_lng0, (VERTEX_ANGLE, -108), VERTEX_ANGLE / CELL_COUNT)

        return (lat_lng)

    @staticmethod
    def get_canonical_west_point(x_index, y_index):
        # Start at canonical south point
        lat_lng0 = Cell.get_canonical_south_point(x_index, y_index)
        # Rotate NW from there by one cell width
        lat_lng = Cell.rotate(lat_lng0, (-VERTEX_ANGLE, -108), VERTEX_ANGLE / CELL_COUNT)

        return (lat_lng)

    @staticmethod
    def get_canonical_north_point(x_index, y_index):
        # Start at canonical east point
        lat_lng0 = Cell.get_canonical_east_point(x_index, y_index)
        # Rotate NW from there by one cell width
        lat_lng = Cell.rotate(lat_lng0, (-VERTEX_ANGLE, -108), VERTEX_ANGLE / CELL_COUNT)

        return (lat_lng)

    @staticmethod
    def get_point_from_canonical(rhomboid, lat_lng):
        if rhomboid < 5:
            # Latitude is correct already
            return (lat_lng[0], lat_lng[1] + rhomboid % 5 * 72)
        else:
            # Start by rotating SE to rhomboid 5
            se_lat_lng = Cell.rotate(lat_lng, (VERTEX_ANGLE, 108), 72)
            # Then do Longitude rotation
            return (se_lat_lng[0], se_lat_lng[1] + 36 + rhomboid % 5 * 72)

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
    def calc_x(clat, clng):
        d0, d3 = Rhomboid.get_x_dist(clat, clng)
#        out = 'd0=%s, d3=%s' % (d0, d3)
#        logging.debug(out)
        if d3 < 1:
            return CELL_COUNT - 1;
        if d0 < 1:
            return 0;        
        edge_fraction = Rhomboid.get_x_edge_fraction(d0, d3)
        cell_index = int(CELL_COUNT * edge_fraction)
        return cell_index - 1
        
        
    @staticmethod
    def calc_y(clat, clng):
        d0, d1 = Rhomboid.get_y_dist(clat, clng)
#        out = 'd0=%s, d1=%s' % (d0, d1)
#        logging.debug(out)
        if d0 < 1:
            return 0
        if d1 < 1:
            return CELL_COUNT - 1
        edge_fraction = Rhomboid.get_y_edge_fraction(d0, d1)
        cell_index = int(CELL_COUNT * edge_fraction)
        return cell_index
        

    def __init__(self, lat, lng):
        self.face = Face(lat, lng)
        self.clat = self.canonical_lat(self.face)
        self.clng = self.canonical_lng(self.face)
        self.num = self.calc_num(self.face)
        self.x = self.calc_x(self.clat, self.clng)
        self.y = self.calc_y(self.clat, self.clng)
        self.key = '%s-%s-%s' % (self.num, self.x, self.y)
