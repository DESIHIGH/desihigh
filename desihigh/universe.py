import numpy as np
import numpy.typing as npt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
from astropy.cosmology import FlatLambdaCDM

def get_age_and_size(model: FlatLambdaCDM):
    """
    Reports the cosmic time vs scale factor relation for a cosmological model.

    Parameters:
    -----------
    model : astropy.cosmology.FlatLambdaCDM
        The astropy cosmology model

    Returns
    -------
    age_of_universe : numpy array of shape (1000,)
        The cosmic time ranging from a redshift range of 0 to 100
    size_of_universe : numpy array of shape (1000,)
        The scale factor ranging from a redshift range of 0 to 100
        
    """

    redshift_range = np.linspace(0, 100, 1000)
    age_of_universe = model.age(redshift_range).value
    size_of_universe = model.scale_factor(redshift_range)
    
    return age_of_universe, size_of_universe

def time_to_redshift(model: FlatLambdaCDM, interp_redshift: npt.ArrayLike, time: npt.ArrayLike):
    """
    Converts cosmic time to redshift, assuming a cosmological model

    Parameters:
    -----------
    model : astropy.cosmology.FlatLambdaCDM
        The astropy cosmology model
    interp_redshift: array-like of shape (N,)
        The sampled redshifts used for interpolating time to redshift, ordered 
        by increasing redshift
    time: array-like of shape (M,)
        The cosmic time that will be interpolated to redshift values

    Returns
    -------
    numpy array of shape (M,)
        The interpolated redshifts
        
    """

    
    t_to_z = InterpolatedUnivariateSpline(model.age(interp_redshift)[::-1], interp_redshift[::-1])

    return t_to_z(time)

def set_y_scale():
    """
    Replaces negative y axis values with positive values for a mirrored graph
        
    """
    
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{abs(x):.3g}'))

def x_position_to_y_position(model: FlatLambdaCDM, redshift_range: npt.ArrayLike, x_position: npt.ArrayLike):
    """
    Generates y-axis locations given a-axis locations for an expansion history 
    plot.
    
    Generates y-axis positions for galaxies in a decorative plot of the 
    universe's expansion history, given the galaxies' x-axis positions
    representing cosmic time.

    Parameters:
    -----------
    model : astropy.cosmology.FlatLambdaCDM
        The astropy cosmology model
    redshift_range : array-like of shape (N,)
        The sampled redshifts used for interpolating time to redshift, ordered 
        by increasing redshift
    x_position : array-like of shape (M,)
        The x-axis locations on the galaxies on the graph, corresponding to 
        cosmic time

    Returns
    -------
    y_position: numpy array of shape (M,)
        The y-axis locations on the galaxies on the graph, corresponding 
        to a random spatial distribution
        
    """
    
    galaxy_redshifts = time_to_redshift(model, redshift_range, x_position)
    
    y_position_bounds = model.scale_factor(galaxy_redshifts)
    
    y_position = np.random.uniform(-y_position_bounds, y_position_bounds)

    return y_position
