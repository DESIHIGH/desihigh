import astropy.io.fits as     fits
import numpy           as     np
import pylab           as     pl

from   astropy.table         import Table
from   scipy.ndimage.filters import gaussian_filter


tiles   = {'bgs': 66003, 'elg': 67230, 'lrgqso': 68002} 

for tracer, band in zip(['bgs', 'elg', 'lrgqso'], ['B', 'Z', 'Z']): 
  zbest = Table.read('../../../andes/zbest-0-{}-20200315.fits'.format(tiles[tracer]))
  coadd = fits.open('../../../andes/coadd-0-{}-20200315.fits'.format(tiles[tracer]))

  assert  np.all(coadd['FIBERMAP'].data['TARGETID'] == zbest['TARGETID'])

  SNRs  = coadd['SCORES'].data['MEDIAN_COADD_SNR_{:s}'.format(band)]
  dChs  = zbest['DELTACHI2']
  
  rank  = np.argsort(dChs)

  lcut  = np.percentile(dChs[rank], 80)
  hcut  = np.percentile(dChs[rank], 85)

  cut   = (dChs[rank] > lcut) & (dChs[rank] < hcut)

  zs    = zbest[rank][cut]
  zs.sort('TARGETID')
  
  tids  = zs['TARGETID']

  print(tracer, len(tids))

  isin   = np.isin(coadd['FIBERMAP'].data['TARGETID'], tids)

  wave    = []
  flux    = []
  
  for arm in ['B', 'R', 'Z']:
    wave += coadd['{:s}_WAVELENGTH'.format(arm)].data.tolist()
    flux += coadd['{:s}_FLUX'.format(arm)].data[isin].T.tolist()
    
  wave    = np.array(wave)
  rank    = np.argsort(wave)
  wave    = wave[rank]

  flux    = np.array(flux)
  flux    = flux[rank, :]

  for i, x in enumerate(flux.T):
    flux[:,i] = gaussian_filter(x, 3)

  # pl.clf()
  # pl.plot(wave, flux, lw=1.0)
  # pl.show()

  print(zbest)
  
print('\n\nDone.\n\n')
