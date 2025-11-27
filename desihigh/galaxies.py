import numpy as np
import numpy.typing as npt
from astropy.cosmology import FlatLambdaCDM

FloatArray = npt.NDArray[np.float32]

def get_ra_dec_z_region(galaxy_file: str):
    """
    Loads sky coordinate data for DESI galaxies in a user-specified path.

    Parameters
    ----------
    galaxy_file: string
        The path to a .BIN file containing the galaxies to be read in. The 
        galaxies should be stored as a one dimensional numpy array, with target 
        IDs followed by R.A. vaues followed by Dec values followed by redshift.

    Returns
    -------
    ra : numpy array of shape (N,)
        The R.A. coordinates in degrees  
    dec : numpy array of shape (N,)
        The declination coordinates in degrees   
    redshift : numpy array of shape (N,)
        The redshift coordinates
        
    """

    galaxies = np.fromfile(galaxy_file)
    num_gals = galaxies.shape[0]//4

    ra = galaxies[num_gals:2*num_gals]

    dec = galaxies[2*num_gals:3*num_gals]

    redshift = galaxies[3*num_gals:]

    return ra, dec, redshift


def get_x_y_z_region(galaxy_file: str, Om0: float = .315, H0: float = 100):
    
    """
    Loads cartiesian coordinate data for DESI galaxies in a user-specified path. 
    
    The coordinates are stored as RA-dec-redshift and are transformed to 
    carteisan coordiantes with a user-specified cosmology.

    Parameters
    ----------
    galaxy_file : string
        The path to a .BIN file containing the galaxies to be read in.  The 
        galaxies should be stored as a one dimensional numpy array, with target 
        IDs followed by R.A. vaues followed by Dec values followed by redshift.
    Om0 : float
        The omega matter value for the specified cosmology. Defaults to 0.315
    H0 : float
        The Hubble constant for the specified cosmology. Defaults to 100 
        km/s/Mpc

    Returns
    -------
    x : numpy array of shape (N,)
        The x coordinates in Mpc (or Mpc/h if H0=100) 
    y : numpy array of shape (N,)
        The y coordinates in Mpc (or Mpc/h if H0=100) 
    z : numpy array of shape (N,)
        The z coordinates in Mpc (or Mpc/h if H0=100)
        
    """
    
    ra, dec, redshift = get_ra_dec_z_region(galaxy_file)

    cosmology_model = FlatLambdaCDM(Om0=Om0, H0=H0)
    distance = cosmology_model.comoving_distance(redshift).value

    x, y, z = ra_dec_dist_to_xyz(ra, dec, distance)

    return x, y, z

def ra_dec_dist_to_xyz(ra_degrees: FloatArray, dec_degrees: FloatArray, distance: FloatArray):

    """
    Transforms sky coordinates to cartesian coordinates

    Parameters
    ----------
    ra_degrees : numpy array of shape (N,)
        The R.A. coordinates in degrees
    dec_degrees : numpy array of shape (N,)
        The declination coordinates in degrees
    distance : numpy array of shape (N,)
        The distance coordinates in user-specified units

    Returns
    -------
    x : numpy array of shape (N,)
        The x coordinates in user-specified units
    y : numpy array of shape (N,)
        The y coordinates in user-specified units 
    z : numpy array of shape (N,)
        The z coordinates in user-specified units
        
    """

    ra = ra_degrees*np.pi/180.
    
    dec = dec_degrees*np.pi/180.
    
    x = distance*np.cos(ra)*np.cos(dec)
    
    y = distance*np.sin(ra)*np.cos(dec)
    
    z = distance*np.sin(dec)
        
    return x, y, z

