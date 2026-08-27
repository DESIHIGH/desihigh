

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

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

def add_spectral_lines (linecolor='k', fontsize=12, text_zorder=10, text_offset = 5, y_height = 0.8):
    
    line_transform = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
    
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
    )

    plt.axvline(486.27, color=linecolor, linestyle=':')
    plt.gca().text(
        486.27-text_offset,
        y_height,
        "H$\\beta$",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
        ha="right",
    )

    plt.axvline(372.7, color=linecolor, linestyle=':')
    plt.gca().text(
        372.7+text_offset,
        y_height,
        "[OII]",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
    )
    plt.axvline(372.9, color=linecolor, linestyle=':')
    plt.gca().text(
        372.9+text_offset,
        y_height -.06,
        "[OII]",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
    )

    plt.axvline(495.9, color=linecolor, linestyle=':')
    plt.gca().text(
        495.9+text_offset,
        y_height,
        "[OIII]",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
    )
    plt.axvline(500.7, color=linecolor, linestyle=':')
    plt.gca().text(
        500.7+text_offset,
        y_height - 0.06,
        "[OIII]",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
    )

    """plt.axvline(422.6727, color=linecolor, linestyle=':')
    plt.gca().text(
        422.6727+text_offset,
        y_height,
        "Ca",
        transform=line_transform,
        fontsize = fontsize,
        zorder = text_zorder,
    )"""


    