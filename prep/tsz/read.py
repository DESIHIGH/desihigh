import glob
import pandas               as pd
import pylab                as pl
import astropy.io.fits      as fits
import numpy                as np 

from   desimodel.io         import load_tiles
from   scipy.spatial        import cKDTree as KDTree
from   desimodel.focalplane import get_tile_radius_deg
from   astropy.table        import Table
from   desimodel.footprint  import is_point_in_desi


dat    = pd.read_csv('pszmmf1.dat', sep='\s+')
ra     = dat.iloc[:,5]
dec    = dat.iloc[:,6] 

tsz    = Table(np.c_[ra, dec], names=['RA', 'DEC'])

## 
files  = glob.glob('/global/cfs/cdirs/desi/spectro/redux/andes/tiles/*/*/cframe-b0-*.fits')
plates = []

for p in files:
    p  = fits.open(p)[0]
    
    plates.append([p.header['TILEID'], p.header['TILERA'], p.header['TILEDEC']])
    
plates = np.array(plates)

plates = Table(plates, names=['TILEID', 'RA', 'DEC'])

# tiles    = load_tiles(extra=True, onlydesi=False)
isin       = is_point_in_desi(plates, tsz['RA'], tsz['DEC'], return_tile_index=False)

# indx     = np.unique(indx).astype(np.int)
# plates   = plates[indx]

desitsz    = tsz[isin]

print(plates)

desitsz.write('andes_tsz.txt', format='ascii', overwrite=True)
