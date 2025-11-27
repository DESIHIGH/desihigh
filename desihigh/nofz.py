import matplotlib.pyplot as plt
import numpy.typing as npt

def get_sky_coverage_fraction (tracer: str):
    """
    Reports the sky coverage fraction for a DESI tracer in DR1

    Parameters
    ----------
    tracer : string
        The DESI tracer. One of 'BGS', 'ELG', 'LRG', or 'QSO'

    Returns
    -------
    float
        The DR1 sky coverage for the tracer
        
    """
    
    # taken from n(z) files at /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/
    # NGC area + SGC area
    if tracer.lower() == 'bgs':
        # bgs
        return (5299.5428 + 2173.1756) / 41252.95
        
    elif tracer.lower() == 'elg':
        # elg
        return (3887.9288 + 2036.108) / 41252.95
        
    elif tracer.lower() == 'lrg':
        #lrg
        return (3755.8352 + 1983.8884) / 41252.95
        
    elif tracer.lower() == 'qso':
        # qso
        return (4644.2888 + 2604.7908) / 41252.95

def plot_galaxy_distribution(bin_centers: npt.ArrayLike, histogram: npt.ArrayLike, label: str=None):
    """
    Plots the n(r) distribution for a galaxy catalog as a histogram

    Parameters
    ----------
    bin_centers : array-like of shape (N,)
        The location of the bin centers
    histogram : array-like of shape (N,)
        The histogram values
    label : string
        A label for the plot. Defaults to None.
        
    """
    plt.step(bin_centers, histogram, where='mid', linewidth = 1, label = label);
