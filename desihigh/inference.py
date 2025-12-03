import numpy as np
import matplotlib.pyplot as plt

from getdist import plots, MCSamples

# Functions to handle under-the-hood operations with getdist 

def plot_chain(sample: MCSamples) -> tuple:
    """
    Plot the chain of values for each parameter.

    Parameters
    ----------
    values : np.ndarray
        Array of values of shape (nparams, nsteps)
    names : list, optional
        List of names of the parameters to display as y-axis labels, by default None

    Returns
    -------
    fig, ax : tuple
        Figure and axes objects for the plot
    """
    values = sample.samples.T # Transpose the values to get shape (nparams, nsteps)
    names = sample.paramNames.list() # Get the parameter names from the sample object
    nparams = len(names)
    
    fig, ax = plt.subplots(nparams, 1, figsize=(5, nparams), sharex=True)
    for i in range(nparams):
        ax[i].plot(values[i], color='black')
        
        ytick_labels = [f'{tick:.2f}' for tick in ax[i].get_yticks()] # yticks precision to 2 decimal places
        ax[i].set_yticks(ax[i].get_yticks())
        ax[i].set_yticklabels(ytick_labels, fontsize=8)
        ax[i].set_ylabel(rf'${sample.parLabel(names[i])}$', fontsize=8)
        
        ax[i].tick_params(axis='x', bottom=False) # Remove x-tick labels for all but the last plot

    ax[-1].set_xlabel('Step')
    ax[-1].tick_params(axis='x', bottom=True, labelsize=8)
    fig.subplots_adjust(hspace=0.2)
    
    return fig, ax

def plot_histogram(sample: MCSamples, name: str, **kwargs):
    """
    Plot the histogram of a single parameter from the MCSamples object.
    """
    g = plots.get_single_plotter()
    g.plot_1d(sample, name, **kwargs)
    return g

def plot_contour(sample: MCSamples, names: list, **kwargs):
    """
    Plot the contour of two parameters from the MCSamples object.
    """
    width_inch = kwargs.pop('width_inch', None)
    ratio = kwargs.pop('ratio', None)
    g = plots.get_single_plotter(width_inch=width_inch, ratio=ratio)
    g.plot_2d(sample, names[0], names[1], **kwargs)
    return g

def plot_triangle(sample: MCSamples, **kwargs):
    """
    Plot the triangle plot of the MCSamples object.
    """
    colors = kwargs.get('colors', ['k'])
    g = plots.get_subplot_plotter()
    g.triangle_plot(
        sample, 
        line_args=[{'color': colors[0]}],
        contour_colors=colors, 
        **kwargs,
    )
    return g