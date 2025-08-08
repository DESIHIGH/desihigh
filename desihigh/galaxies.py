import numpy as np

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def get_ra_dec_z_region():

    galaxies = np.fromfile('./galaxies.BIN')
    num_gals = 28284

    ra = galaxies[num_gals:2*num_gals]

    dec = galaxies[2*num_gals:3*num_gals]

    redshift = galaxies[3*num_gals:]

    return ra, dec, redshift


def get_x_y_z_region(Om0=.315, H0=100):
    
    from astropy.cosmology import FlatLambdaCDM

    ra, dec, redshift = get_ra_dec_z_region()

    cosmology_model = FlatLambdaCDM(Om0=Om0, H0=H0)
    distance = cosmology_model.comoving_distance(redshifts).value

    x, y, z = ra_dec_dist_to_xyz(ra, dec, redshift)

    return x, y, z

def ra_dec_dist_to_xyz(ra_degrees, dec_degrees, distance):

    ra = ra_degrees*np.pi/180.
    
    dec = dec_degrees*np.pi/180.
    
    x = distance*np.cos(ra)*np.cos(dec)
    
    y = distance*np.sin(ra)*np.cos(dec)
    
    z = distance*np.sin(dec)
        
    return x, y, z