import numpy as np
import fitsio
import numpy.lib.recfunctions as rfn
from astropy.table import Table

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def generate_tile_data():

    # right ascension coordinates of tiles for the 2021-09-22 DESI observing plan
    ra = [
       266.0, 261.4, 312.4, 317.0, 318.7, 321.7, 326.3, 327.6, 329.6, 333.0, 336.3, 337.7, 340.1, 
       337.4, 338.6, 336.3, 335.8, 336.1, 335.2, 335.3, 336.6, 96.6, 94.3, 98.5, 102.0,
       103.2, 103.8, 105.3, 106.6, 109.1, 109.3
    ]
    
    # transform the right ascention values to fall between 110 and -100 degrees 
    ra = np.array(ra)
    ra = (ra - 150)%360 + 150 - 360
    
    # declination coordinates of tiles for the 2021-09-22 DESI observing plan
    declination = [
       24.8, 12.8, 0.5, -2.7, 2.3, 4.1, 0.6, -2.5, 0.2, -0.2, -0.6, -6.4, -12.6, 26.1, 31.9, 19.4, 
        15.3, 29.3, 23.4, 9.3, 6.2, 62.6, 65.9, 69.0, 64.9, 61.3, 52.7, 49.5, 43.7 , 35.3, 40.0
    ]
    
    ra.tofile('../data/20210922_tiles_ra.BIN')
    np.array(declination).tofile('../data/20210922_tiles_dec.BIN')
    

def generate_DR1_BGS_sample():

        
    gals = '/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron-main-bright.fits'
    
    # Redshift limits
    zmin = 0.
    zmax = 0.24
    
    
    # Save information about the mask and masked galaxies
    
    #Open the FastSpecFit VAC
    
    with fitsio.FITS(gals) as full_catalog:
        
        print("reading data")
        metadata = full_catalog[2]['TARGETID','Z','ZWARN','DELTACHI2','SPECTYPE','RA','DEC','BGS_TARGET', 'SURVEY', 'PROGRAM'
                           ][:]
        specphot = full_catalog[1][
                           'ABSMAG01_SDSS_R', 
                          ][:]
    
    
        catalog = rfn.merge_arrays([metadata, specphot], flatten=True, usemask=False)
    
        del metadata, specphot
    
        print(len(catalog), "target observations read in")
        
        # Select BGS Bright galaxies
        select = np.where(
                   (catalog['SPECTYPE']=='GALAXY') 
                   & (catalog['SURVEY']=='main')
                   & (catalog['PROGRAM']=='bright')
                 )
        
        catalog=catalog[select]
            
        #check for duplicate targets
        _, select = np.unique(catalog['TARGETID'], return_index=True)
        if len(catalog) != len(select):
            raise ValueError(f'Duplicate galaxies detected. {len(select)} out of {len(catalog)} are unique')
                
        print(len(catalog), "bright time galaxies")
    
        # Impose redshift limits
        print("Imposing redshift limits")
        select = np.where((catalog['Z']>zmin)  # > zmin and not >=zmin to avoid galaxies at origin
                        & (catalog['Z']<=zmax) 
                     )
    
        catalog=catalog[select]
        
        print(len(catalog), "galaxies in redshift limits")
            
        #Quality cuts 
        #made to match Ross 2024, The Dark Energy Spectroscopic Instrument: Construction of Large-scale Structure Catalogs
        select = np.where(
                   (catalog['ZWARN']==0) 
                   & (catalog['DELTACHI2']>40) 
        )   
        
        catalog=catalog[select]
        
        print(len(catalog), "galaxies in final catalog")
            
        #save catalog
        out=Table([catalog['TARGETID'],
                   catalog['RA'],
                   catalog['DEC'],
                   catalog['Z'],
                 
                  ],
                   
                   names=['TARGETID','RA','DEC','Z',]
                 )
    select_region = (out['RA']>190) * (out['RA']<200) * (out['DEC']>-5) * (out['DEC']<5) * (out['Z']>.15) * (out['Z']<.24)
    region = np.array([out['TARGETID'][select], out['RA'][select], out['DEC'][select], out['Z'][select]])
    region.tofile('../data/DR1_BGS_sample_galaxies.BIN')