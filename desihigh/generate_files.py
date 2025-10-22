import pickle
import numpy as np
import numpy.lib.recfunctions as rfn
import fitsio
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM


import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def generate_tile_data(ra_output_file='../data/20210922_tiles_ra.BIN', 
                       dec_output_file='../data/20210922_tiles_dec.BIN'):

    """

    Writes the locations of a night of DESI tiles to file. The tiles are 
    hard-coded to the 2021-09-22 DESI observing plan. This data is used
    in the Mapping The Universe notebook.


    parameters:
    ---------------------------------------------------------------------

    ra_output_file: string
        The path to a .BIN file that will save the R.A. coordinates.
        
    dec_output_file: string
        The path to a .BIN file that will save the declination 
        coordinates.

    """

    # right ascension coordinates of tiles for the 
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
    
    ra.tofile(ra_output_file)
    np.array(declination).tofile(dec_output_file)
    

def generate_DR1_BGS_sample(gals = '/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron-main-bright.fits',
                            output_path = '../data/DR1_BGS_sample_galaxies.BIN',
                            z_min = 0.15,
                            z_max = 0.24,
                            ra_min = 190.,
                            ra_max = 200.,
                            dec_min = -5.,
                            dec_max = 5.,
                            
                           ):

    """

    Writes the locations of a region of DESI BGS galaxies of file. This 
    data is used in the Mapping The Universe and Python Packages 
    notebooks.


    parameters:
    ---------------------------------------------------------------------

    gals: string
        The path to the fastspec catalog that contains the galaxy 
        information. Defaults to the v2.1 fastspec-iron-main-bright.fits
        file on NERSC
        
    output_path: string
        The path to a .BIN file that will save the galaxy coordinates

    z_min: float
        The minimum redshift for the saved region of galaxies. Defaults 
        to 0.15

    z_max: float
        The maximum redshift for the saved region of galaxies. Defaults 
        to 0.24

    ra_min: float
        The minimum R.A. for the saved region of galaxies. Defaults 
        to 190.

    ra_max: float
        The maximum R.A. for the saved region of galaxies. Defaults 
        to 200.

    dec_min: float
        The minimum dec. for the saved region of galaxies. Defaults 
        to -5.

    dec_max: float
        The maximum dec. for the saved region of galaxies. Defaults 
        to 5.

    """
        
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
    
        # Impose survey region limits
        print("Imposing survey region limits")
        select = np.where((catalog['Z']>z_min) 
                        & (catalog['Z']<z_max) 
                        & (catalog['RA']>ra_min)
                        & (catalog['RA']<ra_max)
                        & (catalog['DEC']>dec_min)
                        & (catalog['DEC']<dec_max)
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
    region = np.array([out['TARGETID'], out['RA'], out['DEC'], out['Z']])
    region.tofile(output_path)

def generate_lss_nz(output_path = '../data/lss_catalogs_nz.pickle',
                   bgs_path='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/BGS_ANY_clustering.dat.fits',
                   elg_path='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/ELG_LOPnotqso_clustering.dat.fits',
                   lrg_path='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/LRG_clustering.dat.fits',
                   qso_path='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/QSO_clustering.dat.fits',
                   ):

    """

    Calculates the n(z) of DESI LSS catalogs and saves them to file.


    parameters:
    ---------------------------------------------------------------------

    output_path: string
        The path to a .pickle file that will save the n(z) profile

    bgs_path: string
        The path to a the BGS LSS catalog. Defaults to the v1.5 catalog
        on NERSC

    elg_path: string
        The path to a the ELG LSS catalog. Defaults to the v1.5 catalog
        on NERSC

    lrg_path: string
        The path to a the LRG LSS catalog. Defaults to the v1.5 catalog
        on NERSC

    qso_path: string
        The path to a the QSO LSS catalog. Defaults to the v1.5 catalog
        on NERSC
        
    """

    H0 = 67.4
    csm0 = FlatLambdaCDM(Om0=.315, H0=H0)
    
    bgs = Table.read(bgs_path)
    elg = Table.read(elg_path)
    lrg = Table.read(lrg_path)
    qso = Table.read(qso_path)
    
    bgs_Mlyr = csm0.comoving_distance(bgs['Z']).value * 3.26 # convert to Mlyr
    elg_Mlyr = csm0.comoving_distance(elg['Z']).value * 3.26 # convert to Mlyr
    lrg_Mlyr = csm0.comoving_distance(lrg['Z']).value * 3.26 # convert to Mlyr
    qso_Mlyr = csm0.comoving_distance(qso['Z']).value * 3.26 # convert to Mlyr
    
    bins = np.linspace(0, 15280/(H0/100), 100)
    
    bin_centers = bins[:-1] + np.diff(bins)/2
    
    r_min = bins[:-1]
    r_max = bins[1:]
    
    bgs_Mlyr_hist, _ = np.histogram(bgs_Mlyr, bins=bins)
    lrg_Mlyr_hist, _ = np.histogram(lrg_Mlyr, bins=bins)
    elg_Mlyr_hist, _ = np.histogram(elg_Mlyr, bins=bins)
    qso_Mlyr_hist, _ = np.histogram(qso_Mlyr, bins=bins)
    
    with open(output_path, 'wb') as file:
        pickle.dump((bin_centers, r_min, r_max, bgs_Mlyr_hist, lrg_Mlyr_hist, elg_Mlyr_hist, qso_Mlyr_hist), file)
