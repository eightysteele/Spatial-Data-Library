from gdalconst import GA_ReadOnly
import osgeo.gdal as gdal

def experiment_with_gdal():
    '''Just experimenting with reading files with GDAL.'''
    
    # Opens a tile:
    path = '../../data/bioclim/tiles/bio_37/bio1_37.bil'
    tile = gdal.Open(path, GA_ReadOnly)
    
    # TODO: How to read metadata contained in .hdr via tile?
    
    # Gets bounding box
    # TODO: This is inconsistent with bio1_37.hdr
    gt = tile.GetGeoTransform()
    geotrans = gt
    box = {'maxLat': max(list(geotrans)[1::2]),
           'minLat': min(list(geotrans)[1::2]),
           'maxLon': max(list(geotrans)[0::2]),
           'minLon': min(list(geotrans)[0::2])
           }
    print 'Bounding box: ' + str(box)
    
    # Reads all data:
    # data = tile.ReadRaster(0, 0, tile.RasterXSize, tile.RasterYSize)
    
    # Reads 10x10 of data:
    data = tile.ReadAsArray(0, 0, 10, 10)
    print data


def my_helper(lat, lng):
    print 'lat=%f, lng=%f' % (lat, lng)

def point2key(lat, lng):
    '''Converts the lat/lng into a TGM key and returns it.
    
    Arguments:
        lat - the decimal latitude
        lng - the decimal longitude
    
    Returns:
        The TGM key as a string
    '''
    
    # You can define and call helper methods if it helps:
    my_helper(lat, lng)

    
if __name__ == '__main__':    
    #experiment_with_gdal()
    point2key(37.323345, -4.0334334)
