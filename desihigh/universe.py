import numpy as np
import numpy.typing as npt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
from astropy.cosmology import FlatLambdaCDM

def get_age_and_size(
    model: FlatLambdaCDM,
    redshift_range: np.ndarray = np.linspace(0, 100, 1000),
) -> tuple:
    """
    Reports the cosmic time vs scale factor relation for a cosmological model.

    Parameters:
    -----------
    model : astropy.cosmology.FlatLambdaCDM
        The astropy cosmology model
    redshift_range: np.ndarray
        Redshift range to compute the age and size of the universe for.

    Returns
    -------
    age_of_universe : np.ndarray
        The cosmic time ranging from a redshift range of 0 to 100, 
        with shape matching redshift_range.
    size_of_universe : np.ndarray
        The scale factor ranging from a redshift range of 0 to 100, 
        with shape matching redshift_range.
    """
    age_of_universe = model.age(redshift_range).value
    size_of_universe = model.scale_factor(redshift_range)
    return age_of_universe, size_of_universe

def time_to_redshift(
    model: FlatLambdaCDM, 
    interp_redshift: npt.ArrayLike, 
    time: npt.ArrayLike
) -> np.ndarray:
    """
    Converts cosmic time to redshift, assuming a cosmological model

    Parameters:
    -----------
    model : astropy.cosmology.FlatLambdaCDM
        The astropy cosmology model
    interp_redshift: array-like
        The sampled redshifts used for interpolating time to redshift, 
        ordered by increasing redshift
    time: array-like
        The cosmic time that will be interpolated to redshift values

    Returns
    -------
    np.ndarray
        The interpolated redshifts, with shape matching the time array.
    """
    t_to_z = InterpolatedUnivariateSpline(model.age(interp_redshift)[::-1], interp_redshift[::-1])
    return t_to_z(time)

def set_y_scale():
    """
    Replaces negative y axis values with positive values for a mirrored graph
    """
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{abs(x):.3g}'))

def x_position_to_y_position(
    model: FlatLambdaCDM, 
    redshift_range: npt.ArrayLike, 
    x_position: npt.ArrayLike
) -> np.ndarray:
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
    redshift_range : array-like
        The sampled redshifts used for interpolating time to redshift, ordered 
        by increasing redshift
    x_position : array-like
        The x-axis locations on the galaxies on the graph, corresponding to 
        cosmic time

    Returns
    -------
    y_position: np.ndarray
        The y-axis locations on the galaxies on the graph, 
        corresponding to a random spatial distribution.
        Shape matches the x_positions array.
    """
    
    galaxy_redshifts = time_to_redshift(model, redshift_range, x_position)
    
    y_position_bounds = model.scale_factor(galaxy_redshifts)
    y_position = np.random.uniform(-y_position_bounds, y_position_bounds)

    return y_position
