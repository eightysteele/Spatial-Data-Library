import fixpath
fixpath.fix_sys_path()

from sdl.tmg import gettile
import shapefile

if __name__ == '__main__':
    w = shapefile.Writer(shapefile.POLYGON)
    w.field('CellKey','C','255')    
    for x in gettile((30, 0), (60, -30), 0.008333):
        key = x.cellkey
        parts = [list(x) for x in x.polygon]
        w.poly(parts=[parts])
        w.record(CellKey=key)
    w.save('tmg.shp')
