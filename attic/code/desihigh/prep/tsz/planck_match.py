import glob
import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt
import matplotlib.image  as img 

from   legacy            import cutout


# stampmap_ra213.52_dec71.3.png
paths     = glob.glob('/Users/MJWilson/DESIATHIGHSCHOOL/planck_tsz_cutouts/*.png')
files     = [x.split('/')[-1] for x in paths]

ras       = np.array([np.float(x.split('_')[1].replace('ra', '').replace('.png', '')) for x in files])
decs      = np.array([np.float(x.split('_')[2].replace('dec', '').replace('.png', '')) for x in files])

isin      = (decs > -10.) & (decs < 20.)

ras       = ras[isin]
paths     = np.array(paths)[isin]
decs      = decs[isin]

ngal      = np.count_nonzero(isin)
ngal      = 50

fig, axes = plt.subplots(ngal, 2, figsize=(15., 3.5 * ngal))

for i, (mra, mdec, mpath) in enumerate(zip(ras, decs, paths)): 
  url     = cutout(axes[i,0], mra, mdec)

  im      = img.imread(mpath) 
  
  axes[i,1].imshow(im)
  axes[i,1].set_title('{:d}: {:.3f} {:.3f}'.format(i, mra, mdec))
  
  if i == (ngal - 1):
      break
  
pl.savefig('legacy_cutout.png')
