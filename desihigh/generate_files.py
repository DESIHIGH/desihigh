# This file will contain all the functions created to generate the DESIHIGH files
# Some of those functions might have hardcoded paths or parameters.

import os
import pickle
from glob import glob
from pathlib import Path

import fitsio
import getdist
import numpy as np
import numpy.lib.recfunctions as rfn
from getdist import plots, MCSamples
from scipy.ndimage import gaussian_filter1d
from astropy.table import Table, vstack
from astropy.cosmology import FlatLambdaCDM
from PIL import Image
from matplotlib.image import pil_to_array
import healpy as hp

from desispec.io import read_spectra
from desispec.coaddition import coadd_cameras
from desispec.interpolation import resample_flux
from desispec.resolution import Resolution
import redrock.templates

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def generate_tile_data(
    ra_output_file: str = '../data/20210922_tiles_ra.BIN', 
    dec_output_file: str = '../data/20210922_tiles_dec.BIN',
):
    """
    Writes the locations of a night of DESI tiles to file. 
    
    The tiles are  hard-coded to the 2021-09-22 DESI observing plan. This data 
    is used in the Mapping The Universe notebook.

    Parameters
    ----------
    ra_output_file : str
        The path to a .BIN file that will save the R.A. coordinates.
    dec_output_file : str
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
    ra = (ra - 290)%360 + 290 - 360
    
    # declination coordinates of tiles for the 2021-09-22 DESI observing plan
    declination = [
       24.8, 12.8, 0.5, -2.7, 2.3, 4.1, 0.6, -2.5, 0.2, -0.2, -0.6, -6.4, -12.6, 26.1, 31.9, 19.4, 
        15.3, 29.3, 23.4, 9.3, 6.2, 62.6, 65.9, 69.0, 64.9, 61.3, 52.7, 49.5, 43.7 , 35.3, 40.0
    ]
    
    ra.tofile(ra_output_file)
    np.array(declination).tofile(dec_output_file)
    
def generate_DR1_BGS_sample(
    gals: str = '/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron-main-bright.fits',
    output_path: str = '../data/DR1_BGS_sample_galaxies.BIN',
    z_min: float = 0.15,
    z_max: float = 0.24,
    ra_min: float = 190.,
    ra_max: float = 200.,
    dec_min: float = -5.,
    dec_max: float = 5.,
):
    """
    Writes the locations of a region of DESI BGS galaxies of file. 
    
    This data is used in the Mapping The Universe and Python Packages 
    notebooks.

    Parameters
    ----------
    gals : str
        The path to the fastspec catalog that contains the galaxy information. 
        Defaults to the v2.1 fastspec-iron-main-bright.fits file on NERSC     
    output_path : str
        The path to a .BIN file that will save the galaxy coordinates
    z_min : float
        The minimum redshift for the saved region of galaxies. Defaults 
        to 0.15
    z_max : float
        The maximum redshift for the saved region of galaxies. Defaults 
        to 0.24
    ra_min : float
        The minimum R.A. for the saved region of galaxies. Defaults 
        to 190.
    ra_max : float
        The maximum R.A. for the saved region of galaxies. Defaults 
        to 200.
    dec_min : float
        The minimum dec. for the saved region of galaxies. Defaults 
        to -5.
    dec_max : float
        The maximum dec. for the saved region of galaxies. Defaults 
        to 5.
    """ 
    #Open the FastSpecFit VAC
    with fitsio.FITS(gals) as full_catalog:
        
        print("reading data")
        metadata = full_catalog[2][
            'TARGETID',
            'Z',
            'ZWARN',
            'DELTACHI2',
            'SPECTYPE',
            'RA',
            'DEC',
            'BGS_TARGET', 
            'SURVEY', 
            'PROGRAM',
        ][:]
        specphot = full_catalog[1][
            'ABSMAG01_SDSS_R', 
        ][:]
    
    
        catalog = rfn.merge_arrays([metadata, specphot], flatten=True, usemask=False)
        del metadata, specphot # free memory
    
        print(len(catalog), "target observations read in")
        
        # Select BGS Bright galaxies
        select = np.where(
            (catalog['SPECTYPE']=='GALAXY') &
            (catalog['SURVEY']=='main') &
            (catalog['PROGRAM']=='bright')
        )
        
        catalog=catalog[select]
            
        #check for duplicate targets
        _, select = np.unique(catalog['TARGETID'], return_index=True)
        if len(catalog) != len(select):
            raise ValueError(f'Duplicate galaxies detected. {len(select)} out of {len(catalog)} are unique')
                
        print(len(catalog), "bright time galaxies")
    
        # Impose survey region limits
        print("Imposing survey region limits")
        select = np.where(
            (catalog['Z']>z_min) &
            (catalog['Z']<z_max) &
            (catalog['RA']>ra_min) &
            (catalog['RA']<ra_max) &
            (catalog['DEC']>dec_min) &
            (catalog['DEC']<dec_max)
        )
    
        catalog=catalog[select]
        
        print(len(catalog), "galaxies in redshift limits")
            
        #Quality cuts 
        #made to match Ross 2024, The Dark Energy Spectroscopic Instrument: Construction of Large-scale Structure Catalogs
        select = np.where(
            (catalog['ZWARN']==0)  &
            (catalog['DELTACHI2']>40) 
        )   
        
        catalog=catalog[select]
        
        print(len(catalog), "galaxies in final catalog")
            
        #save catalog
        out=Table(
            [
                catalog['TARGETID'],
                catalog['RA'],
                catalog['DEC'],
                catalog['Z'],
            ],
            names=['TARGETID','RA','DEC','Z',]
        )
    region = np.array([out['TARGETID'], out['RA'], out['DEC'], out['Z']])
    region.tofile(output_path)

def generate_lss_nz(
    output_path: str = '../data/lss_catalogs_nz.pickle',
    bgs_path: str = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/BGS_ANY_clustering.dat.fits',
    elg_path: str = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/ELG_LOPnotqso_clustering.dat.fits',
    lrg_path: str = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/LRG_clustering.dat.fits',
    qso_path: str = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/QSO_clustering.dat.fits',
):
    """
    Calculates the n(z) of DESI LSS catalogs and saves them to file.

    Parameters
    ----------
    output_path : str
        The path to a .pickle file that will save the n(z) profile
    bgs_path : str
        The path to a the BGS LSS catalog. Defaults to the v1.5 catalog on NERSC
    elg_path : str
        The path to a the ELG LSS catalog. Defaults to the v1.5 catalog on NERSC
    lrg_path : str
        The path to a the LRG LSS catalog. Defaults to the v1.5 catalog on NERSC
    qso_path : str
        The path to a the QSO LSS catalog. Defaults to the v1.5 catalog on NERSC
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

def generate_sample_fibers(
    tile_path_petal_0: str = '/global/cfs/cdirs/desi//public/dr1/spectro/redux/iron/tiles/pernight/153/20210504/coadd-0-153-20210504.fits',
    output_path: str = '../data/fibers-153-20210504.fits',
):
    """
    Writes fiber info for a DESI tile to file. 
    
    This data is used in the DESI Survey Overview notebook.

    Parameters
    ----------
    tile_path_petal_0 : str
        The path to the file that cotains the tile/fiber information for petal 0. 
        The petal number in the file path should be denoted with the string '-0-'. 
        Subsequent petals are read in automatically.
    output_path : str
        The path to a .FITS file that will save the combined tile/fiber info.

    """
    # read the first petal
    fiber_table =Table.read(tile_path_petal_0, hdu=1)
    fiber_table=fiber_table['TARGETID', 'FIBER', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PETAL_LOC']
    # read subsequent petals
    for i in range(1,10):
        current_petal = tile_path_petal_0.replace('-0-',f'-{i}-')
        temp_fiber_table =Table.read(current_petal, hdu=1)
        temp_fiber_table=temp_fiber_table['TARGETID', 'FIBER', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PETAL_LOC']
        fiber_table = vstack([fiber_table, temp_fiber_table])

    fiber_table.write(output_path)
    
def gen_black_body(filename: str, save_to: str, shift: float = 0.3) -> None:
    """
    Generate a csv file containing the black body spectrum redshifted by a given value.
    Input file found at https://www.pveducation.org/pvcdrom/appendices/standard-solar-spectra

    Parameters
    ----------
    filename : str
        Path to the input file containing the black body spectrum.
    save_to : str
        Path to the output file where the redshifted spectrum will be saved.
    shift : float, optional
        Redshift value to apply to the spectrum, by default 0.3
    """
    black_body = np.genfromtxt(filename, delimiter=";", names=True)
    flux = black_body['Wm2nm1']
    wavelength = black_body['Wavelength_nm']
    shifted_wavelength = wavelength * np.sqrt(1+shift/(1-shift))
    
    mask = (wavelength < 4000) # Keep only wavelengths below 4000 nm
    wavelength = wavelength[mask]
    flux = flux[mask]
    shifted_wavelength = shifted_wavelength[mask]

    stack = np.column_stack((wavelength, shifted_wavelength, flux))
    headers = "wavelength_nm, shifted_wavelength_nm, flux_Wm2nm1"
    np.savetxt(save_to, stack, delimiter=",", header=headers)

def get_desi_spectrum(
    targetid: int, 
    specprod: str = 'iron', 
    dir_from_prod: str = 'healpix/main/dark/21/2196', # DR1 example directory structure
    smoothing_sigma: int = 1.0,
    save_dir: str | None = None,
):
    """
    Get the observed spectrum and best-fit redrock model for a given DESI target ID.
    Requires to be run at NERSC in the desimodules environment with access to the DESI spectroscopic reduction products.
    Inspired from: https://github.com/desihub/timedomain/blob/master/desitrip/docs/nb/RedrockResiduals.ipynb
    
    Parameters:
    -----------
    targetid: int
        The DESI target ID for which to retrieve the spectrum and model.
    specprod: str, optional
        The DESI spectroscopic reduction product (e.g., 'iron', 'fuji').
    dir_from_prod: str, optional
        The directory path from the specprod root to the coadd and zbest files (e.g., 'healpix/main/dark/21/2196').
    smoothing_sigma: float, optional
        The sigma for Gaussian smoothing applied to the spectra for better visualization (in pixels).
    save_dir: str or None, optional
        If provided, the directory where the raw and model spectra will be saved as text files. If None, the spectra will not be saved to files.
    
    Returns:
    --------
    wave: 1D numpy array
        The observed wavelength array (in Angstroms).
    flux: 1D numpy array
        The observed flux array (in erg s^-1 cm^-2 Angstrom^-1),
        smoothed with a Gaussian kernel for better visualization.
    txflux: 1D numpy array
        The best-fit redrock model flux array (in erg s^-1 cm^-2 Angstrom^-1),
        convolved with the instrument resolution and smoothed for visualization.
    z: float
        The redshift of the best-fit model.
    
    Example:
    --------
    Using the target ID 39628462612284326 from DESI, you can call the function as follows:
    >>> wave, flux, txflux, z = get_desi_spectrum(39628462612284326, specprod='iron', dir_from_prod='healpix/main/dark/21/2196')
    
    This will return the spectrum and model for: https://www.legacysurvey.org/viewer/desi-spectrum/dr1/targetid39628462612284326
    """
    # Load Redrock Templates
    templates = dict()
    for filename in redrock.templates.find_templates():
        t = redrock.templates.Template(filename)
        templates[(t.template_type, t.sub_type)] = t
    
    # Get the coadd and zbest files for the specified product and directory.
    coadd_pattern = '/'.join([os.environ['DESI_SPECTRO_REDUX'], specprod, dir_from_prod, 'coadd*.fits'])
    zbest_pattern = '/'.join([os.environ['DESI_SPECTRO_REDUX'], specprod, dir_from_prod, 'redrock*.fits'])
    
    # Loop over coadd and zbest files to find the spectrum for the given targetid.
    coadd_files = sorted(glob(coadd_pattern))
    zbest_files = sorted(glob(zbest_pattern))
    for cafile, zbfile in zip(coadd_files, zbest_files):
        # Access data per petal.
        zbest = Table.read(zbfile, hdu=1).filled('')
        pspectra = read_spectra(cafile)      # coadded exposures, separate cameras
        cspectra = coadd_cameras(pspectra)   # coadd B, R, Z
        fibermap = cspectra.fibermap

        # Check if the targetid is in this petal's fibermap.
        
        if targetid in fibermap['TARGETID']:
            break # Found the petal containing the targetid, exit loop.
        
    idx = np.where(fibermap['TARGETID'] == targetid)[0][0]
    
    z = zbest['Z'][idx]
    wave = cspectra.wave['brz']
    flux = cspectra.flux['brz'][idx]
    res = cspectra.resolution_data['brz'][idx]
    
    spectype = zbest['SPECTYPE'][idx].strip() if zbest['SPECTYPE'][idx] else ''
    subtype = zbest['SUBTYPE'][idx].strip() if zbest['SUBTYPE'][idx] else ''
    
    # Evaluate the best-fit template at the observed wavelengths, applying redshift and resolution.
    fulltype = (spectype, subtype)
    ncoeff = templates[fulltype].flux.shape[0]
    coeff = zbest['COEFF'][idx][0:ncoeff] 
    tflux = templates[fulltype].flux.T.dot(coeff)
    twave = templates[fulltype].wave * (1 + z)
    
    R = Resolution(res)
    txflux = R.dot(resample_flux(wave, twave, tflux))
    
    # Apply small gaussian smoothing to the fluxes for better visualization.
    flux = gaussian_filter1d(flux, sigma=smoothing_sigma) 
    txflux = gaussian_filter1d(txflux, sigma=smoothing_sigma)
    
    if save_dir is not None:
        save_dir = Path(save_dir)
        
        # Save raw data to a file with a header containing the target information.
        header = f'TARGETID={targetid}, Z={z:.3f}, SPECTYPE={spectype}, SUBTYPE={subtype}, smoothing={smoothing_sigma}\n'
        header += 'wave_A flux\n'
        np.savetxt(save_dir / 'desi_spectra_data.txt', np.column_stack([wave, flux]), header=header)

        # Save model data to a file with a header containing the model information.
        header = f'Model for TARGETID={targetid} at z=0, SPECTYPE={spectype}, SUBTYPE={subtype}, smoothing={smoothing_sigma}\n'
        header += 'wave_A flux\n'
        wave_model_z0 = wave / (1 + z) # Shift observed wavelengths to z=0 frame
        np.savetxt(save_dir / 'redrock_model_data.txt', np.column_stack([wave_model_z0, txflux]), header=header)
    
    return wave, flux, txflux, z

CHAIN_MAPPING = {
    'omegam': 'Omega_m',
    'w': 'w0_fld',
    'wa': 'wa_fld',
    'H0rdrag': 'H0_rd',
    'rdrag': 'rd',
}

def get_desi_chain(
    chains_dir: str = '/global/cfs/cdirs/desi/public/papers/y3/bao-cosmo-params/cobaya', 
    model: str = 'base_w_wa', 
    dataset: str = 'desi-bao-all_CMB-compressed-theta-ombh2-ombch2_desy5sn',
    parameters = ['age', 'rdrag', 'H0rdrag', 'omegam', 'w', 'wa'],
    mapping = CHAIN_MAPPING,
    save_dir: str | None = None,
) -> MCSamples:
    """
    Gets a DESI cosmological parameter chain from the specified directory, extracts the specified parameters, and optionally saves the new chain to a file.
    Data is from: https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html
    
    Parameters
    ----------
    chains_dir : str, optional
        The directory containing the cosmological parameter chains, by default '/global/cfs/cdirs/desi/public/papers/y3/bao-cosmo-params/cobaya'
    model : str, optional
        The model name, by default 'base_w_wa'
    dataset : str, optional
        The dataset name, by default 'desi-bao-all_CMB-compressed-theta-ombh2-ombch2_desy5sn'
    parameters : list, optional
        The parameters to extract from the chain, by default ['age', 'rdrag', 'H0rdrag', 'omegam', 'w', 'wa']
    mapping : _type_, optional
        A mapping from the original parameter names to the desired names, by default CHAIN_MAPPING
    save_dir : str | None, optional
        The directory to save the new chain to, by default None

    Returns
    -------
    MCSamples
        A new MCSamples object containing only the specified parameters, with names mapped according to the provided mapping. 
        If save_dir is not None, the new chain is also saved to a .npy file in the specified directory.
    """
    chains_dir = Path(chains_dir) # Convert to Path object
    
    chain = getdist.loadMCSamples(str(chains_dir / model / dataset / 'chain'))

    # Print the names of the parameters in the chain
    print("Parameters in the chain:")
    print(chain.getParamNames().list())

    # Create a new McSamples object with the selected parameters
    idx = [chain.getParamNames().list().index(param) for param in parameters]
    names = [mapping.get(param, param) for param in parameters]
    labels = [chain.parLabel(param) for param in parameters]
    new_chain = MCSamples(
        samples = chain.samples[:, idx],
        names = names,
        labels = labels,
    )

    if save_dir is not None:
        save_dir = Path(save_dir) # Convert to Path object
        np.save(save_dir / 'desi_dr2_chain.npy', new_chain) # Save the new chain to a file
    
    return new_chain

def generate_earth_healpix(image_path='../images/MTU-Earth.jpg', output_path='../data/earth_map_healpix.BIN', nside = 256):
    """
    Converts a image of Earth to a healpix grayscale representation. 
    
    Based on tutorial at https://www.zonca.dev/posts/2013-08-08-healpix-map-of-earth-using-healpy

    Parameters
    ----------
    image_path : str
        The path to an image file to convert to healpix 
        representation.
    output_path : str
        The path to a .BIN file that will save the healpix pixel 
        representation.
    nside : int
        The healpix pixelation NSIDE parameter. Defaults to 256
    """

    earth_array = pil_to_array(Image.open(image_path).convert("L"))

    # grid of coorinates to sample from image
    theta = np.linspace(0, np.pi, earth_array.shape[0])[:, None]
    phi = np.linspace(-np.pi, np.pi, earth_array.shape[1])

    # convert coordinates to pixels
    pix = hp.ang2pix(nside, theta, phi)

    # convert to healpix map
    earth_healpix = np.zeros(hp.nside2npix(nside), dtype=np.float32)
    earth_healpix[pix] = earth_array

    # adjust grayscale to fit well with the gist_earth color map
    earth_healpix = earth_healpix**(.75)

    earth_healpix.tofile(output_path)
