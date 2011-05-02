import fixpath
fixpath.fix_sys_path()

from sdl import tmg
import shapefile

if __name__ == '__main__':
    w = shapefile.Writer(shapefile.POLYGON)
    for key in tmg.get_tile([-120, 40], [-110, 37]):
        pass
    
