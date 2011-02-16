/**
 * Copyright 2007 University of California at Berkeley.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */

public abstract class TMG {
  // Digits of precision are important to maintain for these static variables.

  // cellcount is the number of cells along each axis in the TMG
  static final int cellcount = 4;

  // phi is (1+Math.sqrt(5))/2;
  static final double phi = 1.618033988749895;

  // Multiply value in degrees by radians to convert degrees to radians.
  // Divide value in radians by radians to convert radians to degrees.
  // = Math.PI/180;
  // = 0.017453292519943295 radians/degree
  static final double radians = 0.017453292519943295;

  // semimajoraxis is the equatorial radius of WGS84 in meters to the nearest
  // 0.1 meters.
  // = 6378137.0 meters
  static final double semimajoraxis = 6378137.0;

  // vertexangle is the angle in degrees from the center of the sphere
  // between two adjacent vertexes.
  // = 2*Math.asin(1/(Math.sqrt(phi*phi+1)))/radians;
  // = 63.434948822922 degrees
  static final double vertexangle = 63.434948822922;

  // edgelength is the direct distance between adjacent vertexes
  // = 2*semimajoraxis/Math.sqrt(phi*phi+1);
  // = 6706370.116516389 m for WGS84.
  static final double edgelength = 6706370.116516389;

  // midedgeradius is the distance from the centroid to the middle of an edge.
  // = edgelength*phi/2;
  // = 5425567.394830056 for WGS84
  static final double midedgeradius = 5425567.394830056;

  // ==========
  // facecenterlat is the latitude of the center of each of the 5 northern
  // faces in degrees. This puts the latitudes of the centers of the
  // equatorial faces at +- (90-facecenterlat) = +- 10.8123171
  // The distance to each vertex from the center of the face on the surface of
  // the sphere is 4160829.589 m.
  // = Math.asin((phi*phi)/(Math.sqrt(3)*Math.sqrt(phi*phi+1)))/radians;
  // = 52.622631859350314;
  static final double facecenterlat = 52.622631859350314;

  // centerangle is the angle in degrees from the center of the sphere
  // between the centers of two adjacent faces.
  // = 2*Math.acos(phi/Math.sqrt(3))/radians;
  // = 41.810314895778575
  static final double centerangle = 41.810314895778575;

  // vertexlat is the latitude in degrees N or S for all non-polar vertexes.
  // = 90-vertexangle;
  // = 26.565051177077997;
  static final double vertexlat = 26.565051177077997;

  // cos72 = Math.cos(-72 * radians);
  // = 0.30901699437494734;
  static final double cos72 = 0.30901699437494734;
  // sin72 = Math.sin(-72 * radians);
  // = -0.9510565162951536;
  static final double sin72 = -0.9510565162951536;
  // cosvertexlat=Math.cos(-vertexlat * radians);
  // = 0.8944271909999159;
  static final double cosvertexlat = 0.8944271909999159;
  // sinvertexlat=Math.sin(-vertexlat * radians);
  // = -0.4472135954999581;
  static final double sinvertexlat = -0.4472135954999581;

  public static double canonicalLatitude(double lat, double lng) {
    int face = getFace(lat, lng);
    // For faces [0-9], canonical latitude is the same as latitude
    if (face < 10) {
      return lat;
    }
    // Get the canonical longitude.
    double clng = canonicalLongitude(face, lng);
    double newlat = 0;
    // double newlng = 0;
    // Original Cartesian coordinates of lat, lng.
    double x = 0, y = 0, z = 0;

    // For the southern faces [10-19], canonical latitude can be determined by
    // right-handed rotating face 10 or 15 by 72 degrees around the axis passing
    // through -vertexlat, -72.
    // The rotation axis is as follows where lat and lng are in radians:
    // c1=cos(lat)cos(long)
    // c2=cos(lat)sin(long)
    // c3=sin(lat)
    double c1 = cosvertexlat * cos72;
    double c2 = cosvertexlat * sin72;
    double c3 = sinvertexlat;
    double cosa = cos72;
    double sina = -sin72;

    x = Math.cos(lat * radians) * Math.cos(clng * radians);
    y = Math.cos(lat * radians) * Math.sin(clng * radians);
    z = Math.sin(lat * radians);

    // Cartesian coordinates after rotation.
    // double x1 = x * cosa + (1 - cosa)
    // * (c1 * c1 * x + c1 * c2 * y + c1 * c3 * z) + (c2 * z - c3 * y) * sina;
    // double y1 = y1 = y * cosa + (1 - cosa)
    // * (c2 * c1 * x + c2 * c2 * y + c2 * c3 * z) + (c3 * x - c1 * z) * sina;
    double z1 = z * cosa + (1 - cosa)
        * (c3 * c1 * x + c3 * c2 * y + c3 * c3 * z) + (c1 * y - c2 * x) * sina;
    // newlat and newlng are the equivalent to lat and lng in face 0
    newlat = Math.asin(z1) / radians;
    // newlng = Math.atan2(y1, x1) / radians;
    return newlat;
  }

  public static double canonicalLongitude(double lat, double lng) { // checked
    int face = getFace(lat, lng);
    double clng = lng360(lng);
    if (face < 5) {
      clng -= face * 72;
    } else if (face < 10) {
      clng = clng - (face - 5) * 72;
    } else if (face < 15) {
      clng = clng - 36 - (face - 10) * 72;
    } else {
      clng = clng - 36 - (face - 15) * 72;
    }
    return lng360(clng);
  }

  // Direct distance between to points on a sphere.
  public static double cartesianDistance(double fromlat, double fromlng,
      double tolat, double tolng, double radius) { // checked
    double from_x = radius * Math.cos(fromlat * radians)
        * Math.cos(fromlng * radians);
    double from_y = radius * Math.cos(fromlat * radians)
        * Math.sin(fromlng * radians);
    double from_z = radius * Math.sin(fromlat * radians);

    double to_x = radius * Math.cos(tolat * radians)
        * Math.cos(tolng * radians);
    double to_y = radius * Math.cos(tolat * radians)
        * Math.sin(tolng * radians);
    double to_z = radius * Math.sin(tolat * radians);

    return Math.sqrt((from_x - to_x) * (from_x - to_x) + (from_y - to_y)
        * (from_y - to_y) + (from_z - to_z) * (from_z - to_z));
  }

  public static double centerlat(int face) { // checked
    if (face < 5) {
      return facecenterlat;
    }
    if (face < 10) {
      return facecenterlat - centerangle;
    }
    if (face < 15) {
      return -(facecenterlat - centerangle);
    }
    return -facecenterlat;
  }

  public static double centerlng(int face) { // checked
    if (face < 10) {
      return (face % 5) * 72;
    }
    return lng360(36 + (face % 5) * 72);
  }

  // The face is one of the sides of an icosahedron number consecutively from 0
  // beginning with the northern face having its center at lng=0.
  // Return the face of the icosahedron onto which the lat lng projects.
  public static int getFace(double lat, double lng) { // checked
    int face = 0;
    double clng = lng180(lng);
    if (lat >= vertexlat) {
      // One of the northern five faces.
      // 0=[-36,36} 1=[36,108}, 2=[108,180} 3=[180,-108}, 4=[-108,-36}
      face = (int) Math.floor((clng + 360 + 36) / 72);
      if (clng >= -36) {
        face -= 5;
      }
    } else if (lat < -vertexlat) {
      // One of the southern five faces.
      // 15=[0,72} 16=[72,144}, 17=[144,-144} 18=[-144,-72}, 19=[-72,0}
      face = (int) (Math.floor((clng + 360) / 72)) + 15;
      if (clng >= 0) {
        face -= 5;
      }
    } else {

      // lat, lng is somewhere in the equatorial zone.
      // Find which of ten equal sections around the equator lng falls into
      // by rotating lng to equivalent longitude at lat=0.
      int section = 0;
      if (clng >= -18) {
        section = (int) Math.floor((18 + clng + 18 * (lat) / vertexlat) / 36);
      } else {
        section = (int) Math.floor((18 + 360 + clng + 18 * (lat) / vertexlat) / 36);
      }
      // Even numbered sections are northern facets [5,9].
      if (section % 2 == 0) {
        face = 5 + section / 2;
      } else {
        // Odd numbered sections are southern facets [10,14].
        face = 10 + (section - 1) / 2;
      }
    }
    return face;
  }

  // The facet is a spherical rhombus consisting of two adjacent faces, one of
  // which has a vertex at a pole. The sphere consists of ten such facets.
  public static int getRhomboid(double lat, double lng) {
    int face = getFace(lat, lng);
    if (face < 10) {
      return face % 5;
    }
    return face % 5 + 5;
  }

  // RhomboidX is the index of the cell along the SE or NW edge of the
  // equivalent
  // canonical facet.
  // Use only canonicalLatitude and canonicalLongitude in calls to getRhomboidX.
  public static int getRhomboidX(double clat, double clng) {
    // double clat = canonicalLatitude(lat, lng);
    // double clng = canonicalLongitude(getFace(lat, lng), lng);
    double d0 = 0;
    double d3 = 0;
    if (clat < vertexlat) { // face 5, use SE edge for i index
      // d3 is the distance from the canonical lat lng to the E vertex
      d3 = cartesianDistance(clat, clng, vertexlat, 36, semimajoraxis);
      if (d3 < 1) {
        return 0;
      }
      // d0 is the distance from the canonical lat lng to the S vertex
      d0 = cartesianDistance(clat, clng, -vertexlat, 0, semimajoraxis);
      if (d0 < 1) {
        return cellcount - 1;
      }
    } else { // face 0, use NW edge for i index
      // d0 is the distance from the canonical lat lng to the W vertex
      d0 = cartesianDistance(clat, clng, vertexlat, -36, semimajoraxis);
      if (d0 < 1) {
        return 0;
      }
      // d3 is the distance from the canonical lat lng to the N vertex
      d3 = cartesianDistance(clat, clng, 90, 0, semimajoraxis);
      if (d3 < 1) {
        return cellcount - 1;
      }
    }

    double distfromd0 = (edgelength * edgelength - d3 * d3 + d0 * d0)
        / (2 * edgelength);
    double h = Math.sqrt(d0 * d0 - distfromd0 * distfromd0);
    double dp = h * Math.tan(30 * radians);
    double edgeangle = 0;
    distfromd0 += dp;
    if (distfromd0 > d0 / 2) {
      edgeangle = vertexangle / 2;
      edgeangle += Math.atan((distfromd0 - (d0 / 2)) / midedgeradius) / radians;
    } else {
      edgeangle = vertexangle / 2;
      edgeangle -= Math.atan(((d0 / 2) - distfromd0) / midedgeradius) / radians;
    }
    double edgefraction = edgeangle / vertexangle;
    double cell = cellcount * edgefraction;
    int index = (int) cell;
    return index - 1;
  }

  // RhomboidY is the index of the cell along the NE or SW edge of the
  // equivalent
  // canonical facet.
  // Use only canonicalLatitude and canonicalLongitude in calls to getRhomboidY.
  public static int getRhomboidY(double clat, double clng) {
    // double clat = canonicalLatitude(lat, lng);
    // double clng = canonicalLongitude(getFace(lat, lng), lng);
    double d0;
    double d1;
    if (clat < vertexlat) { // face 5, use SW edge for j index
      // d0 is the distance from the canonical lat lng to the S vertex
      d0 = Math.floor(cartesianDistance(clat, clng, -vertexlat, 0,
          semimajoraxis));
      if (d0 < 1) {
        return 0;
      }
      // d1 is the distance from the canonical lat lng to the W vertex
      d1 = Math.floor(cartesianDistance(clat, clng, vertexlat, -36,
          semimajoraxis));
      if (d1 < 1) {
        return cellcount - 1;
      }
    } else { // face 0, use NE edge for j index
      d0 = Math.floor(cartesianDistance(clat, clng, vertexlat, 36,
          semimajoraxis));
      if (d0 < 1) {
        return 0;
      }
      // d1 is the distance from the canonical lat lng to the N vertex
      d1 = Math.floor(cartesianDistance(clat, clng, 90, 0, semimajoraxis));
      if (d1 < 1) {
        return cellcount - 1;
      }
    }
    double distfromd0 = (edgelength * edgelength - d1 * d1 + d0 * d0)
        / (2 * edgelength);
    double h = Math.sqrt(d0 * d0 - distfromd0 * distfromd0);
    double dp = h * Math.tan(30 * radians);
    double edgeangle = 0;
    distfromd0 += dp;
    distfromd0 = Math.floor(distfromd0);
    if (distfromd0 > d0 / 2) {
      edgeangle = vertexangle / 2;
      edgeangle += Math.atan((distfromd0 - (d0 / 2)) / midedgeradius) / radians;
    } else {
      edgeangle = vertexangle / 2;
      edgeangle -= Math.atan(((d0 / 2) - distfromd0) / midedgeradius) / radians;
    }
    double edgefraction = edgeangle / vertexangle;
    double cell = cellcount * edgefraction;
    int index = (int) cell;
    return index;
  }

  public static String getTMGKey(double lat, double lng) {
    // TMG is a Triangular Mesh Grid system for a sphere based on the
    // faces of an inscribed icosahedron.
    // A facet is a spherical rhomboid defined by two adjacent faces of
    // an icosahedron circumscribed by the sphere. Here the orientation
    // of the icosahedron is defined to have two vertexes at the
    // poles and the remaining ten vertexes arranged every 36 degrees of
    // longitude with the vertex of facet 0 (face 5) at longitude=0,
    // latitude=-26.565051177077997.

    int rhomboid = getRhomboid(lat, lng);

    // Get equivalent lat, long for a canonical face (0 or 5).
    double clat = canonicalLatitude(lat, lng);
    double clng = canonicalLongitude(lat, lng);

    // Get the x and y offsets from vertex 0 of the rhombus.
    int RhomboidX = getRhomboidX(clat, clng);
    int RhomboidY = getRhomboidY(clat, clng);

    return rhomboid + "-" + RhomboidX + "-" + RhomboidY;
  }

  // Distance along great circle on sphere between two points using the law of
  // cosines.
  // public static double greatCircleDistance(double lat1, double lng1,
  // double lat2, double lng2, double radius) {
  // // latitudes and longitudes assumed to be degrees
  // return Math.acos(Math.sin(lat1 * radians) * Math.sin(lat2 * radians)
  // + Math.cos(lat1 * radians) * Math.cos(lat2 * radians)
  // * Math.cos((lng2 - lng1) * radians))
  // * radius;
  // }

  // public static int index(double percent) {
  // double cellvalue = cellcount * percent;
  // int index = (int) cellvalue;
  // double indexfloor = Math.floor(cellvalue);
  // int indexfloorint = (int) indexfloor;
  // BigDecimal bd = BigDecimal.valueOf(cellvalue);
  // int indexbd = bd.intValue();
  // System.out.println("\npercent: " + percent + " cellcount: " + cellcount);
  // System.out.println("cellvalue: " + cellvalue);
  // System.out.println("index: " + index);
  // System.out.println("indexfloor: " + indexfloor);
  // System.out.println("indexfloorint: " + indexfloorint);
  // System.out.println("bd: " + bd);
  // System.out.println("indexbd: " + indexbd);
  // return 0;
  // }

  // Return a longitude in {-180, 180]
  public static double lng180(double lng) {
    if (lng <= -180) {
      return lng + 360;
    }
    if (lng > 180) {
      return lng - 360;
    }
    return lng;
  }

  // Return a longitude in [0, 360}
  public static double lng360(double lng) {
    if (lng < 0) {
      return lng + 360;
    }
    if (lng > 360) {
      return lng - 360;
    }
    return lng;
  }

  public static void main(String args[]) {
    // System.out.println("phi: " + phi + "\nvertexangle: " + vertexangle
    // + "\nvertexlat: " + vertexlat + "\ncenterangle: " + centerangle
    // + "\nfacecenterlat: " + facecenterlat + "\nedgelength: " + edgelength
    // + "\nmidedgeradius: " + midedgeradius + "\nfacecenterlat: "
    // + facecenterlat + "\ncosvertexlat: " + cosvertexlat
    // + "\nsinvertexlat: " + sinvertexlat);
    long starttime = System.currentTimeMillis();

    double testlat = 0;
    double testlng = 0;
    // int xcellcount = 360 * 120;
    int xcellcount = 1;
    int ycellcount = 180 * 120;
    for (int j = 0; j < xcellcount; j++) {
      testlng = -180 + j * 360 / xcellcount;
      for (int i = 0; i < ycellcount; i++) {
        testlat = -90 + i * 180 / ycellcount;
        System.out.println("key: " + getTMGKey(testlat, testlng) + " lat: "
            + testlat + " lng: " + testlng);
      }
    }
    long endtime = System.currentTimeMillis();
    System.out.println("Run time: " + (endtime - starttime) + " ms");
  }
}