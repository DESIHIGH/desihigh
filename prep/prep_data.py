import astropy.io.fits as     fits
import numpy           as     np
import pylab           as     pl

from   astropy.table         import Table, join, vstack
from   scipy.ndimage.filters import gaussian_filter
from   desitarget.cmx        import cmx_targetmask

# Coma on petal 0 of 70510.
tiles   = {'mws': 66003, 'bgs': 66003, 'elg': 67230, 'lrg': 68002, 'qso': 68002} 

for tracer, band in zip(['mws', 'bgs', 'elg', 'lrg', 'qso'], ['B', 'B', 'Z', 'Z', 'Z']): 
  zbest = Table.read('../../../andes/zbest-0-{}-20200315.fits'.format(tiles[tracer]))
  coadd = fits.open('../../../andes/coadd-0-{}-20200315.fits'.format(tiles[tracer]))

  assert  np.all(coadd['FIBERMAP'].data['TARGETID'] == zbest['TARGETID'])
  
  tinfo = Table(coadd['FIBERMAP'].data)['TARGETID', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'CMX_TARGET', 'TARGET_RA', 'TARGET_DEC']
  zbest = join(zbest, tinfo, join_type='left', keys='TARGETID')

  if tracer == 'mws':  
    zbest = zbest[(zbest['CMX_TARGET'] & cmx_targetmask.cmx_mask.mask('SV0_MWS')) != 0]
    zbest = zbest[(zbest['SPECTYPE'] == 'STAR')]

  elif tracer == 'bgs':
    zbest = zbest[(zbest['CMX_TARGET'] & cmx_targetmask.cmx_mask.mask('SV0_BGS')) != 0]
    zbest = zbest[(zbest['SPECTYPE'] == 'GALAXY')]

  elif tracer in ['elg', 'lrg']:
    zbest = zbest[(zbest['CMX_TARGET'] & cmx_targetmask.cmx_mask.mask('SV0_{}'.format(tracer.upper()))) != 0]
    zbest = zbest[(zbest['SPECTYPE'] == 'GALAXY')]

  else:
    zbest = zbest[(zbest['CMX_TARGET'] & cmx_targetmask.cmx_mask.mask('SV0_QSO')) != 0]
    zbest = zbest[(zbest['SPECTYPE'] == 'QSO')]
      
  dChs  = zbest['DELTACHI2']
  
  rank  = np.argsort(dChs)

  lcut  = np.percentile(dChs[rank], 60)
  hcut  = np.percentile(dChs[rank], 95)

  cut   = (dChs[rank] > lcut) & (dChs[rank] < hcut)

  zs    = zbest[rank][cut]
  zs.sort('TARGETID')

  for x in ['NPIXELS', 'NUMEXP', 'NUMTILE', 'NCOEFF']:
    del zs[x]

  assert  np.all(zs['ZWARN'] == 0)

  # Limit number of rows.
  zs = zs[:5]
  
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

  result               = Table()
  result['WAVELENGTH'] = wave

  for i, x in enumerate(flux.T):
    result['TARGET{:d}'.format(tids[i])] = flux[:,i]
    
  # print(result)

  result.write('../student_andes/coadd-{}-{}-20200315.fits'.format(tracer, tiles[tracer]), format='fits', overwrite=True)  
  
print('\n\nDone.\n\n')
