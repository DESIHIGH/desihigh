

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np

def angle_degrees_formatter(tick_value, tick_position):
    """
    Custom matplotlib formatter for depicting graph axes with
    modular units that wrap to 0 at 360 degrees

    Parameters
    ----------
    tick_value : float
        The tick label as an angle.

    tick_position: float
        The tick position (unused)

    Returns
    -------
    float
        Reformatted tick as an angle between 0 and 360 degrees
    """
    
    if tick_value < 0:
        return f"{tick_value + 360:.0f}"
    else:
        return f"{tick_value:.0f}"

def add_spectral_lines (min_wavelength, max_wavelength, linecolor='k', fontsize=12, text_zorder=10, text_offset = 2.5, y_height = 0.8):
    
    line_transform = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)

    bbox=dict(
        facecolor="white",  
        alpha=0.7,          
        edgecolor="none",   
        boxstyle="round,pad=0.0"  
    )

    if min_wavelength < 656.46 and max_wavelength > 656.46:
        
        plt.axvline(656.46, color=linecolor, linestyle=':')
        plt.gca().text(
            656.46+text_offset,
            y_height,
            "H$\\alpha$",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            #ha="center",
            #va="bottom",
            bbox=bbox,
        )

    if min_wavelength < 486.27 and max_wavelength > 486.27:
            
        plt.axvline(486.27, color=linecolor, linestyle=':')
        plt.gca().text(
            486.27-text_offset,
            y_height,
            "H$\\beta$",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            ha="right",
            bbox=bbox,
        )

    if min_wavelength < 372.7 and max_wavelength > 372.9:
    
        plt.axvline(372.7, color=linecolor, linestyle=':')
        plt.gca().text(
            372.7+text_offset,
            y_height,
            "[OII]",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )
        plt.axvline(372.9, color=linecolor, linestyle=':')
        plt.gca().text(
            372.9+text_offset,
            y_height -.06,
            "[OII]",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )

    if min_wavelength < 495.9 and max_wavelength > 500.7:
        
        plt.axvline(495.9, color=linecolor, linestyle=':')
        plt.gca().text(
            495.9+text_offset,
            y_height,
            "[OIII]",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )
        plt.axvline(500.7, color=linecolor, linestyle=':')
        plt.gca().text(
            500.7+text_offset,
            y_height - 0.06,
            "[OIII]",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )
        
    if min_wavelength < 393.4 and max_wavelength > 396.9:
        
        plt.axvline(393.4, color=linecolor, linestyle=':')
        plt.gca().text(
            393.4+text_offset,
            y_height -.06*2,
            "Ca II K",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )

        plt.axvline(396.9, color=linecolor, linestyle=':')
        plt.gca().text(
            396.9+text_offset,
            y_height - 0.06*3,
            "Ca II H",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )

    if min_wavelength < 517.5 and max_wavelength > 517.5:
        
        plt.axvline(517.5, color=linecolor, linestyle=':')
        plt.gca().text(
            517.5+text_offset,
            y_height -.06*2,
            "Mg b",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            bbox=bbox,
        )

    if min_wavelength < 121.6 and max_wavelength > 121.6:
        
        plt.axvline(121.6, color=linecolor, linestyle=':')
        plt.gca().text(
            121.6+text_offset,
            y_height,
            "Ly$\\alpha$",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            #ha="center",
            #va="bottom",
            bbox=bbox,
        )

    if min_wavelength < 190.9 and max_wavelength > 190.9:
        
        plt.axvline(190.9, color=linecolor, linestyle=':')
        plt.gca().text(
            190.9+text_offset,
            y_height,
            "C III]",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            #ha="center",
            #va="bottom",
            bbox=bbox,
        )

    if min_wavelength < 154.9 and max_wavelength > 154.9:
        
        plt.axvline(154.9, color=linecolor, linestyle=':')
        plt.gca().text(
            154.9+text_offset,
            y_height,
            "C IV",
            transform=line_transform,
            fontsize = fontsize,
            zorder = text_zorder,
            #ha="center",
            #va="bottom",
            bbox=bbox,
        )

        



def gaussian(x, mu, sigma, amp):
    """
    Calculate an unnormalized Gaussian function

    Parameters
    ----------
    x : float or array of flaots
        The values input to the Gaussian function

    mu: float
        The mean of the Gaussian function

    signma: float
        The standar deviation of the Gaussian function

    amp: float
        The amplitude of the Gaussian function

    Returns
    -------
    float
        The evaulated Gaussian function
    """
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
