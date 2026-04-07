import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects


def add_lat(latitudes, rotate=None, color='gray', coord='earth') -> 'PathCollection' :
    """
    Adds lines of lattitude to a healpy graph.

    Parameters:
    -----------
    latitudes : arraylike of floats
        The lattitudes to plot
    rotate : healpy Rotator
        The rotation applied to the of the map. Defaults to None,
        which indicates a rotation of [0, 0, 0]
    color : string
        The color of the latitide line, If the color is gray, the 
        alpha value is set to 0.2 and the equator is bolded. If set 
        to any other color, the alpha value is set to 1.
    coord : string
        Either 'earth' or 'sky'. Indicates the direction of east and
        west on the map.

    Returns
    -------
    plot : PathCollection
        The scatter plot that constitutes the lattitude lines.
    """

    if coord == 'earth':
        flip = 1
    elif coord == 'sky':
        flip = -1
    else:
        err_txt = f'invalid cood value: {coord}'
        raise ValueError(err_txt)

    for latitude in latitudes:
    
        lons = np.linspace(-180, 180, 1000)
        lats = latitude + np.zeros(1000)

        if rotate is None:
            lons_transform, lats_transform = lons, lats
        else:
            lons_transform, lats_transform = rotate(lons, lats, lonlat=True)
    
        
        if color == 'gray':
            linewidth = 3 if np.isclose(latitude%360, 0) else 1
            plot = plt.scatter(flip*np.deg2rad(lons_transform), np.deg2rad(lats_transform), color='gray', alpha=.2, s=linewidth)
        else:
            plot = plt.scatter(flip*np.deg2rad(lons_transform), np.deg2rad(lats_transform), color=color, s=2)
    return plot

def add_lon(longitudes, rotate=None, color='gray', coord='earth') -> 'PathCollection' :
    """
    Adds lines of longitude to a healpy graph.

    Parameters:
    -----------
    longitudes : arraylike of floats
        The longitudes to plot
    rotate : healpy Rotator
        The rotation applied to the of the map. Defaults to None,
        which indicates a rotation of [0, 0, 0]
    color : string
        The color of the longitude line, If the color is gray, the 
        alpha value is set to 0.2 and the prime meridian is bolded. 
        If set to any other color, the alpha value is set to 1.
    coord : string
        Either 'earth' or 'sky'. Indicates the direction of east and
        west on the map.

    Returns
    -------
    plot : PathCollection
        The scatter plot that constitutes the longitude lines.
    """
    
    if coord == 'earth':
        flip = 1
    elif coord == 'sky':
        flip = -1
    else:
        err_txt = f'invalid cood value: {coord}'
        raise ValueError(err_txt)

    for longitude in longitudes:

        if coord == 'sky':
            longitude=(longitude+180)%360 - 180
            
        lons = longitude + np.zeros(1000) 
        lats = np.linspace(-90, 90, 1000) 
        
        if rotate is None:
            lons_transform, lats_transform = lons, lats
        else:
            lons_transform, lats_transform = rotate(lons, lats, lonlat=True)

        if color == 'gray':
            linewidth = 3 if np.isclose(longitude%360, 0) else 1
            plot = plt.scatter(flip*np.deg2rad(lons_transform), np.deg2rad(lats_transform) , color='gray', alpha=.2, s=linewidth)
        else:
            plot = plt.scatter(flip*np.deg2rad(lons_transform), np.deg2rad(lats_transform), color=color, s=2)
    return plot

def add_coordinate(latitude, longitude, rotate, zorder, coord='earth'):
    """
    Highlights a coordinate on a healpy graph with a star.

    Parameters:
    -----------
    latitude : float
        The latitude to plot
    longitude : float
        The longitude to plot
    rotate : healpy Rotator
        The rotation applied to the of the map. Defaults to None,
        which indicates a rotation of [0, 0, 0]
    coord : string
        Either 'earth' or 'sky'. Indicates the direction of east and
        west on the map.
    """
    if coord == 'earth':
        text_longitude = f'long = {longitude}°'
        text_latitude = f'lat = {latitude}°'
        flip=1
    elif coord == 'sky':
        text_longitude = f'RA = {longitude}°'
        text_latitude = f'dec = {latitude}°'
        longitude=(longitude+180)%360 - 180
        flip=-1
    else:
        err_txt = f'invalid cood value: {coord}'
        raise ValueError(err_txt)

    if rotate is None:
        longitude_transform, latitude_transform = longitude, latitude
    else:
        longitude_transform, latitude_transform = rotate(longitude, latitude, lonlat=True)
    plt.scatter(flip*np.deg2rad([longitude_transform]), np.deg2rad([latitude_transform]), color='red', marker='*', s=500, edgecolors='black', zorder=zorder+1);

    if rotate is None:
        text_longitude_transform, text_latitude_transform = longitude, latitude-10
    else:
        text_longitude_transform, text_latitude_transform = rotate(longitude, latitude-10, lonlat=True)
    long_label = plt.text(flip*np.deg2rad([text_longitude_transform]), np.deg2rad([text_latitude_transform]), text_longitude,color='cyan', weight='bold');
    
    if rotate is None:
        text_longitude_transform, text_latitude_transform = longitude+10*flip, latitude
    else:
        text_longitude_transform, text_latitude_transform = rotate(longitude+10*flip, latitude, lonlat=True)
        
    lat_label = plt.text(flip*np.deg2rad([text_longitude_transform]), np.deg2rad([text_latitude_transform]), text_latitude,color='magenta', weight='bold');
    
    long_label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])
    lat_label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])

def plot_earth_map(earth_map, longitude, latitude):
    """
    Plots a map of the Earth with labeled graticules and a 
    highlighted coordinate.

    Parameters:
    -----------
    earth_map : numpy array
        The healpix array representing the earth
    longitude : float
        The longitude to highlight
    latitude : float
        The latitude to highlight
    """

    if longitude < -180 or longitude > 180:
        print(f'Make sure your longitude value is between -180 and 180 before making your map of Earth! The longitude value you provided is {longitude}')
        return
    if latitude < -90 or latitude > 90:
        print(f'Make sure your latitude value is between -90 and 90 before making your map of Earth! The latitude value you provided is {latitude}')
        return

    hp.newvisufunc.projview(earth_map, cmap="gist_earth", flip="geo", cbar=None,
                        coord='C',
                        graticule=True, #turn on for graticule labels
                        graticule_labels=True,
                        #graticule_coord='C',
                        phi_convention="counterclockwise",
                        fontsize={'xlabel':15,'ylabel':15,'xtick_label':15,'ytick_label':15},
                        custom_xtick_labels=["-120°", "-60°", "0°", "60°", "120°"],
                       )

    
    # user-chosen coordinates
    add_lat([latitude], None, 'magenta')
    await_z_order = add_lon([longitude], None, 'cyan')
    
    zorder = await_z_order.get_zorder()
    
    add_coordinate(latitude, longitude, None, zorder)
    
    # plot settings
    axis = plt.gca()
    
    axis.set_xlabel('longitude (long.) [degrees]', color='cyan', weight='bold')
    axis.set_ylabel('latitude (lat.) [degrees]', color='magenta', weight='bold')
    for label in axis.get_xticklabels():
        label.set_fontweight('bold')
        label.set_color('cyan')
        label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])
    for label in axis.get_yticklabels():
        label.set_fontweight('bold')
        label.set_color('magenta')


def plot_sky_map(sky_map, longitude, latitude):
    """
    Plots a map of the night sky with labeled graticules and a 
    highlighted coordinate.

    Parameters:
    -----------
    earth_map : numpy array
        The healpix array representing the night sky
    longitude : float
        The longitude to highlight
    latitude : float
        The latitude to highlight
    """


    if longitude < 0 or longitude > 360:
        print(f'Make sure your right ascension value is between 0 and 360 before making your map of the sky! The right ascension value you provided is {longitude}')
        return
    if latitude < -90 or latitude > 90:
        print(f'Make sure your declination value is between -90 and 90 before making your map of the sky! The declination value you provided is {latitude}')
        return

    hp.newvisufunc.projview(sky_map, cbar=None,
                        graticule=True, #turn on for graticule labels
                        graticule_labels=True,
                        coord=('G','C'),
                        fontsize={'xlabel':15,'ylabel':15,'xtick_label':15,'ytick_label':15},
                        custom_xtick_labels=["120°", "60°", "0°", "300°", "240°"],
                        min=-.15,
                        max=2,
                        cmap = mpl.colormaps['bone']
                       )
    
    # user-chosen coordinates
    add_lat([latitude], None, 'magenta', coord='sky')
    await_z_order = add_lon([longitude], None, 'cyan', coord='sky')
    
    zorder = await_z_order.get_zorder()
    
    add_coordinate(latitude, longitude, None, zorder, coord='sky')
    
    # plot settings
    axis = plt.gca()
    
    axis.set_xlabel('right ascension (RA) [degrees]', color='cyan', weight='bold')
    axis.set_ylabel('declination (Dec.) [degrees]', color='magenta', weight='bold')
    for label in axis.get_xticklabels():
        label.set_fontweight('bold')
        label.set_color('cyan')
        label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])
    for label in axis.get_yticklabels():
        label.set_fontweight('bold')
        label.set_color('magenta')

def plot_transformed_earth_map(earth_map, longitude, latitude, long_offset, lat_offset, psi_offset):
    """
    Plots a rotated map of the Earth with unlabeled graticules and a 
    highlighted coordinate.

    Parameters:
    -----------
    earth_map : numpy array
        The healpix array representing the earth
    longitude : float
        The longitude to highlight
    latitude : float
        The latitude to highlight
    long_offset : float
        The longitude to rotate the map by
    lat_offset : float
        The latitude to rotate the map by
    psi_offset : float
        The psi coordinate normal to the Earth's surface to rotate 
        the map by
    """

    if longitude < -180 or longitude > 180:
        print(f'Make sure your longitude value is between -180 and 180 before making your map of Earth! The longitude value you provided is {longitude}')
        return
    if latitude < -90 or latitude > 90:
        print(f'Make sure your latitude value is between -90 and 90 before making your map of Earth! The latitude value you provided is {latitude}')
        return
    
    rotate = hp.Rotator(rot=(long_offset, lat_offset, psi_offset), deg=True)
    
    longitude_transform, latitude_transform = rotate(longitude, latitude, lonlat=True)
    
    hp.newvisufunc.projview(earth_map, cmap="gist_earth", flip="geo", cbar=None,
                            coord='C',
                            #graticule=True, #turn on for graticule labels
                            #graticule_labels=True,
                            rot=[long_offset, lat_offset, psi_offset],
                            #rot_graticule=True, 
                            #graticule_coord='C',
                            phi_convention="counterclockwise",
                            fontsize={'xlabel':15,'ylabel':15,'xtick_label':15,'ytick_label':15},
                            custom_xtick_labels=["-120°", "-60°", "0°", "60°", "120°"],
                           )
    
    #graticules
    add_lat([-60,-30,0,30,60], rotate)
    add_lon([-120,-60,0,60,120, 180], rotate)
    
    # user-chosen coordinates
    add_lat([latitude], rotate, 'magenta')
    await_z_order = add_lon([longitude], rotate, 'cyan')
    
    zorder = await_z_order.get_zorder()
    
    add_coordinate(latitude, longitude, rotate, zorder)
    
    # plot settings
    axis = plt.gca()
    
    axis.set_xlabel('longitude (long.) [degrees]', color='cyan', weight='bold')
    axis.set_ylabel('latitude (lat.) [degrees]', color='magenta', weight='bold')
    for label in axis.get_xticklabels():
        label.set_fontweight('bold')
        label.set_color('cyan')
        label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])
    for label in axis.get_yticklabels():
        label.set_fontweight('bold')
        label.set_color('magenta')

def plot_transformed_sky_map(sky_map, longitude, latitude, long_offset, lat_offset, psi_offset):
    """
    Plots a rotated map of the nihgt sky with unlabeled graticules and a 
    highlighted coordinate.

    Parameters:
    -----------
    sky_map : numpy array
        The healpix array representing the night sky
    longitude : float
        The longitude to highlight
    latitude : float
        The latitude to highlight
    long_offset : float
        The longitude to rotate the map by
    lat_offset : float
        The latitude to rotate the map by
    psi_offset : float
        The psi coordinate normal to the projected plane of the sky to 
        rotate the map by
    """
    if longitude < 0 or longitude > 360:
        print(f'Make sure your right ascension value is between 0 and 360 before making your map of the sky! The right ascension value you provided is {longitude}')
        return
    if latitude < -90 or latitude > 90:
        print(f'Make sure your declination value is between -90 and 90 before making your map of the sky! The declination value you provided is {latitude}')
        return


    rotate = hp.Rotator(rot=(long_offset, lat_offset, psi_offset), deg=True)
    
    longitude_transform, latitude_transform = rotate(longitude, latitude, lonlat=True)
    
    hp.newvisufunc.projview(sky_map, cbar=None,
                            coord=('G','C'),
                            rot=[long_offset, lat_offset, psi_offset],
                            min=-.15,
                            max=2,
                            cmap = mpl.colormaps['bone']
                           )

    #graticules
    add_lat([-60,-30,0,30,60], rotate, coord='sky')
    add_lon([120, 60, 0, 300, 240, 180], rotate, coord='sky')
    
    # user-chosen coordinates
    add_lat([latitude], rotate, 'magenta', coord='sky')
    await_z_order = add_lon([longitude], rotate, 'cyan', coord='sky')
    
    zorder = await_z_order.get_zorder()
    
    add_coordinate(latitude, longitude, rotate, zorder, coord='sky')
    
    # plot settings
    axis = plt.gca()
    
    axis.set_xlabel('right ascension (RA) [degrees]', color='cyan', weight='bold')
    axis.set_ylabel('declination (Dec.) [degrees]', color='magenta', weight='bold')
    for label in axis.get_xticklabels():
        label.set_fontweight('bold')
        label.set_color('cyan')
        label.set_path_effects([patheffects.Stroke(linewidth=3, foreground='black'),
                                 patheffects.Normal()])
    for label in axis.get_yticklabels():
        label.set_fontweight('bold')
        label.set_color('magenta')