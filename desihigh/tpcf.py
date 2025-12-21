import yaml
import numpy as np
from numba import njit

def load_subway_data(city: str):
    """
    Load subway station data for a given city.

    Parameters
    ----------
    city : str
        The name of the city to load data for. Supported cities are 'boston', 'london', and 'paris'.

    Returns
    -------
    np.ndarray
        A 2D NumPy array containing the (longitude, latitude) coordinates of the subway stations.

    Raises
    ------
    ValueError
        If the city is not supported.
    """
    if city.lower() == 'boston':
        with open("../../../data/subway_boston.yaml", 'r') as stream: # source of MBTA https://erikdemaine.org/maps/mbta/
            boston_dict = yaml.safe_load(stream)

        # Dictionary to store unique stations
        # Key = station title, Value = (longitude, latitude)
        unique_stations = {}
    
        for line in boston_dict:
            stations = line.get('stations', [])
            for station in stations:
                title = station.get('title', 'Unknown')
                longitude = station.get('longitude')
                latitude = station.get('latitude')
                if title not in unique_stations:
                    unique_stations[title] = (longitude, latitude)
        
        coords_list = list(unique_stations.values()) # Extract the (longitude, latitude) pairs into a list
        data = np.array(coords_list) # Convert the list to a 2D NumPy array
    elif city.lower() == 'london':
        data = np.genfromtxt("../../../data/subway_london.csv",delimiter=',',skip_header=1,usecols=(1,2))[:, [1, 0]] # long first, lat second
    elif city.lower() == 'paris':
        data = np.genfromtxt('../../../data/subway_paris.csv', delimiter=';', skip_header=1, usecols=(5,4,-1))
    else:
        raise ValueError("City not supported.")
    return data

def geo_to_cartesian(lon, lat, R: float = 6371.0):
    """
    Convert geographic coordinates (longitude, latitude) to Cartesian coordinates (x, y) using a small angle approximation.

    Parameters
    ----------
    lon : np.ndarray
        Longitude coordinates in degrees.
    lat : np.ndarray
        Latitude coordinates in degrees.
    R : float
        Radius of the Earth in kilometers. Default is 6371.0 km.

    Returns
    -------
    np.ndarray
        Cartesian coordinates (x, y) in kilometers.
    """
    phi = np.radians(lon)
    theta = np.radians(lat)
    
    # Center to 0,0
    phi -= np.mean(phi)
    theta -= np.mean(theta)
    
    # Using small angle approximation :
    x = R * phi
    y = R * theta
    return np.column_stack((x, y))

@njit # This makes the function faster
def pair_count_2d(positions, edges):
    """
    Count the number of pairs of points with separations in given bins.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (N, 2) containing the (x, y) coordinates of N points.
    edges : np.ndarray
        Array of shape (M,) containing the bin edges for the histogram.

    Returns
    -------
    np.ndarray
        Array of shape (M-1,) containing the counts of pairs in each bin.
    """
    counts = np.zeros(len(edges) - 1)
    for i in range(positions.shape[0]):
        for j in range(i + 1, positions.shape[0]): # avoid double counting
            dx = positions[i, 0] - positions[j, 0]
            dy = positions[i, 1] - positions[j, 1]
            dist2 = dx * dx + dy * dy
            if dist2 < edges[-1]**2:  # Only count if within the maximum distance
                # Find the index in the edges array
                idx = int((np.sqrt(dist2) - edges[0]) / (edges[-1] - edges[0]) * len(counts))
                counts[idx] += 1
    return counts

def generate_randoms(data_positions, factor=10, seed=42, match_data=False):
    """
    Generate random points within the bounding box of the given data.

    Parameters
    ----------
    data : np.ndarray
        A NumPy array where each the two first columns represent the (longitude, latitude) coordinates of the data points.
    factor : int, optional
        The factor by which to increase the number of random points compared to the original data points. Default is 10.
    seed : int, optional
        The random seed for reproducibility. Default is 42.

    Returns
    -------
    np.ndarray
        A 2D NumPy array containing the generated random points.
    """
    rng = np.random.RandomState(seed=seed)
    
    nrandoms = len(data_positions) * factor
    # Generate random positions within the bounding box of the data
    data_min = np.min(data_positions, axis=0)
    data_max = np.max(data_positions, axis=0)
    randoms_positions = np.column_stack([rng.uniform(data_min[i], data_max[i], size=nrandoms) for i in [0, 1]])
    
    if match_data:
        # Create 2D data histogram
        hist, xedges, yedges = np.histogram2d(data_positions[:, 0], data_positions[:, 1], bins=20)
        hist = hist / hist.max()  # make it a probability
        # Bin randoms in the 2D space
        xidx = np.digitize(randoms_positions[:, 0], bins=xedges) - 1
        yidx = np.digitize(randoms_positions[:, 1], bins=yedges) - 1
        # Downsample randoms to match the amplitude of the 2D data hist
        prob = rng.uniform(0., 1., nrandoms)
        mask = prob < hist[xidx, yidx]
        randoms_positions = randoms_positions[mask]
    
    return randoms_positions

