'''
Generate publication ready plots
'''
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from   pkg_resources   import resource_filename

def set_size(width=240, fraction=2):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5 ** 0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt

    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    return  (fig_width_in, fig_height_in)


plt.style.use(resource_filename('tools', 'style.mplstyle'))

nice_fonts = {
    "text.usetex": False,
    "font.family": "serif",
    "font.serif" : "Times New Roman",
}

plt.rcParams["figure.figsize"] = set_size()

matplotlib.rcParams.update(nice_fonts)

