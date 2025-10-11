import numpy as np
from astropy.cosmology import FlatLambdaCDM

def get_ra_dec_z_region(galaxy_file):

    """

    galaxy_file: string
        The path to a .BIN file containing the galaxies to be read in. 
        The galaxies should be stored as a one dimesnional numpy array,
        with target IDs followed by RA vaues followed by Dec values
        followed by redshift.

    """

    galaxies = np.fromfile(galaxy_file)
    num_gals = galaxies.shape[0]//4

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

