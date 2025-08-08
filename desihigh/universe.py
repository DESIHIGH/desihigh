import numpy as np

def get_age_and_size(model):

    redshift_range = np.linspace(0, 100, 1000)
    age_of_universe = model.age(redshift_range).value
    size_of_universe = model.scale_factor(redshift_range)
    
    return age_of_universe, size_of_universe

def time_to_redshift(model, interp_redshift, time):

    from scipy.interpolate import InterpolatedUnivariateSpline
    t_to_z = InterpolatedUnivariateSpline(model.age(interp_redshift)[::-1], interp_redshift[::-1])

    return t_to_z(time)

def set_y_scale():
    from matplotlib.ticker import FuncFormatter
    from matplotlib import pyplot as plt
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{abs(x):.3g}'))

def x_position_to_y_position(model, redshift_range, x_position):
    
    galaxy_redshifts = time_to_redshift(model, redshift_range, x_position)
    
    y_position_bounds = model.scale_factor(galaxy_redshifts)
    
    y_position = np.random.uniform(-y_position_bounds, y_position_bounds)

    return y_position