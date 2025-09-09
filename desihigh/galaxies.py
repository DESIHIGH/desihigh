import numpy as np
import numpy.lib.recfunctions as rfn
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import fitsio
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def get_ra_dec_z_region():

    galaxies = np.fromfile('../data/DR1_BGS_sample_galaxies.BIN')
    num_gals = 28284

    ra = galaxies[num_gals:2*num_gals]

    dec = galaxies[2*num_gals:3*num_gals]

    redshift = galaxies[3*num_gals:]

    return ra, dec, redshift


def get_x_y_z_region(Om0=.315, H0=100):
    

    ra, dec, redshift = get_ra_dec_z_region()

    cosmology_model = FlatLambdaCDM(Om0=Om0, H0=H0)
    distance = cosmology_model.comoving_distance(redshift).value

    x, y, z = ra_dec_dist_to_xyz(ra, dec, distance)

    return x, y, z

def ra_dec_dist_to_xyz(ra_degrees, dec_degrees, distance):

    ra = ra_degrees*np.pi/180.
    
    dec = dec_degrees*np.pi/180.
    
    x = distance*np.cos(ra)*np.cos(dec)
    
    y = distance*np.sin(ra)*np.cos(dec)
    
    z = distance*np.sin(dec)
        
    return x, y, z


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