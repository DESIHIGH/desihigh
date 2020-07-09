import astropy.io.fits as     fits
import numpy           as     np
import pylab           as     pl

from   astropy.table         import Table, join
from   scipy.ndimage.filters import gaussian_filter


tiles   = {'bgs': 66003, 'elg': 67230, 'lrgqso': 68002} 

for tracer, band in zip(['bgs', 'elg', 'lrgqso'], ['B', 'Z', 'Z']): 
  zbest = Table.read('../../../andes/zbest-0-{}-20200315.fits'.format(tiles[tracer]))
  coadd = fits.open('../../../andes/coadd-0-{}-20200315.fits'.format(tiles[tracer]))

  assert  np.all(coadd['FIBERMAP'].data['TARGETID'] == zbest['TARGETID'])

  tinfo = Table(coadd['FIBERMAP'].data)['TARGETID', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET']
  zbest = join(zbest, tinfo, join_type='left', keys='TARGETID')

  # No stars to start. 
  zbest = zbest[(zbest['SPECTYPE'] == 'GALAXY')]
  
  dChs  = zbest['DELTACHI2']
  
  rank  = np.argsort(dChs)

  lcut  = np.percentile(dChs[rank], 80)
  hcut  = np.percentile(dChs[rank], 85)

  cut   = (dChs[rank] > lcut) & (dChs[rank] < hcut)

  zs    = zbest[rank][cut]
  zs.sort('TARGETID')

  for x in ['NPIXELS', 'NUMEXP', 'NUMTILE', 'NCOEFF']:
    del zs[x]

  assert  np.all(zs['ZWARN'] == 0)

  print('\n\n')
  print(zs)
  
  zs.write('../student_andes/zbest-{}-{}-20200315.fits'.format(tracer, tiles[tracer]), format='fits', overwrite=True)
  
  tids    = zs['TARGETID']

  isin    = np.isin(coadd['FIBERMAP'].data['TARGETID'], tids)

  assert  np.all(coadd['FIBERMAP'].data['TARGETID'][isin] == tids)
  
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

  result  = Table()
  result['WAVE'] = wave

  for i, x in enumerate(flux.T):
    result['TARGET{:d}'.format(tids[i])] = flux[:,i]
    
  # print(result)

  result.write('../student_andes/coadd-{}-{}-20200315.fits'.format(tracer, tiles[tracer]), format='fits', overwrite=True)  
  
print('\n\nDone.\n\n')
