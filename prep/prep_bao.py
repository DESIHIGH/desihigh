import os
import glob
import astropy.io.fits as fits
import numpy as np

from   astropy.table import Table, join

# https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/sweep/
#
# os.system('wget https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/sweep/8.0/sweep-000m005-010p000.fits')
# os.system('wget https://portal.nersc.gov/project/cosmo/data/legacysurvey/dr8/south/sweep/8.0-photo-z/sweep-000m005-010p000-pz.fits')
#
# Moved both to ../../../legacy_sweeps

sweep = fits.open('../../../legacy_sweeps/sweep-000m005-010p000.fits')
pzs   = fits.open('../../../legacy_sweeps/sweep-000m005-010p000-pz.fits')

# Line matched.
assert  (sweep[1].data.shape == pzs[1].data.shape)

# Cut to those with valid pzs (NOBSR > 1 in GRZ).
# np.unique(pzs[1].data['z_phot_mean'])
# array([-9.9000000e+01,  5.5704173e-03,  6.6326247e-03, ...,
#        2.3876209e+00,  2.3987792e+00,  2.4036493e+00], dtype=float32)

dsweep = sweep[1].data
dpzs   = pzs[1].data

isin   = (dpzs['z_phot_mean'] > 0.0)

dpzs   = dpzs[isin]
dsweep = dsweep[isin]

# Detected in Z for removing mag. errors.
isin   = dsweep['FLUX_Z'] > 0.0

dpzs   = dpzs[isin]
dsweep = dsweep[isin]

# Photometric Redshift files @ 
# https://www.legacysurvey.org/dr8/files/#photometric-redshifts

# z-band magnitude brighter than 21 are mostly reliable.
zmag   = 22.5 - 2.5 * np.log10(dsweep['FLUX_Z'] / dsweep['MW_TRANSMISSION_Z'])

isin   = (zmag < 21.1)

dpzs   = dpzs[isin]
dsweep = dsweep[isin]

# Morphology cut for star-galaxy separation. 
isin   = dsweep['TYPE'] != 'PSF'

dpzs   = dpzs[isin]
dsweep = dsweep[isin]

dpzs   = Table(dpzs)
dsweep = Table(dsweep)

dpzs['ID']   = np.arange(len(dpzs))
dsweep['ID'] = np.arange(len(dpzs))

assert  (len(dsweep) == len(dpzs))

cat          = join(dpzs, dsweep, join_type='left', keys='ID') 
cat          = cat['RA', 'DEC', 'z_phot_mean']

cat['REDSHIFT'] = cat['z_phot_mean']

del cat['z_phot_mean']

# Redshift limits.
isin         = (cat['REDSHIFT'] > 0.6) & (cat['REDSHIFT'] < 1.05) 

print(cat)

cat.write('../dat/lrg_cat.fits', format='fits', overwrite=True)

indices      = np.random.choice(len(cat), 50, replace=False)

cat          = cat[indices]

cat['Dec']   = cat['DEC']

del cat['DEC']

cat['COLOR'] = 'FUCHSIA'
# cat['ALPHA'] =    0.2

print(cat)

cat.write('../dat/lrg_cat_viewer.fits', format='fits', overwrite=True)

#
# np.savetxt('../dat/lrg_cat_viewer.txt', cat['RA', 'Dec'])
