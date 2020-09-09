import glob
import astropy.io.fits as fits
import pylab           as pl 

from   astropy.table   import Table


specs  = glob.glob('boss_spec/0p96/*.fits')

fluxs  = Table()
models = Table() 

zbest  = Table()

for tid, ss in enumerate(specs):
    ss     = fits.open(ss)

    # ss.info()

    coadd  = ss['COADD']
    spall  = ss['SPALL']

    # spall.header
    # 'PLUG_RA', 'PLUG_DEC', 'Z', 'flux', 'loglam', 'model'

    pra    = spall.data['PLUG_RA']
    pdec   = spall.data['PLUG_DEC']
    pz     = spall.data['Z']

    loglam = coadd.data['loglam']

    flux   = coadd.data['flux']
    model  = coadd.data['model']

    if tid == 0:
      fluxs['WAVE']  = 10. ** loglam
      models['WAVE'] = 10. ** loglam

      zbest['Z']     = pz
      zbest['RA']    = pra
      zbest['Dec']   = pdec

    else:
        # 
        # fluxs['FLUX_{}'.format(tid)]   = flux

        #
        # models['MODEL_{}'.format(tid)] = model

        zbest.add_row((pz, pra, pdec))
    
    print('Solved for {} of {}'.format(tid, len(specs)))

##                                                                                                                                                                                                                   
pl.hist(zbest['Z'], bins=50)
pl.show()

zbest = zbest[(zbest['Z'] > 0.3) & (zbest['Z'] < 0.7)]    
zbest.write('tsz_boss_cluster_zs.fits', format='fits', overwrite=True)
