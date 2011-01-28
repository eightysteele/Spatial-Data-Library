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

if __name__ == '__main__':    
    experiment_with_gdal()
