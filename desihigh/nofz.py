from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import pickle
import numpy as np
import matplotlib.pyplot as plt

def generate_lss_nz(output_path = '../data/lss_catalogs_nz.pickle'):

    H0 = 67.4
    csm0 = FlatLambdaCDM(Om0=.315, H0=H0)
    
    bgs = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/BGS_ANY_clustering.dat.fits')
    elg = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/ELG_LOPnotqso_clustering.dat.fits')
    lrg = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/LRG_clustering.dat.fits')
    qso = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/QSO_clustering.dat.fits')
    
    bgs_Mlyr = csm0.comoving_distance(bgs['Z']).value * 3.26 # convert to Mlyr
    elg_Mlyr = csm0.comoving_distance(elg['Z']).value * 3.26 # convert to Mlyr
    lrg_Mlyr = csm0.comoving_distance(lrg['Z']).value * 3.26 # convert to Mlyr
    qso_Mlyr = csm0.comoving_distance(qso['Z']).value * 3.26 # convert to Mlyr
    
    bins = np.linspace(0, 15280/(H0/100), 100)
    
    bin_centers = bins[:-1] + np.diff(bins)/2
    
    r_min = bins[:-1]
    r_max = bins[1:]
    
    bgs_Mlyr_hist, _ = np.histogram(bgs_Mlyr, bins=bins)
    lrg_Mlyr_hist, _ = np.histogram(lrg_Mlyr, bins=bins)
    elg_Mlyr_hist, _ = np.histogram(elg_Mlyr, bins=bins)
    qso_Mlyr_hist, _ = np.histogram(qso_Mlyr, bins=bins)
    
    with open(output_path, 'wb') as file:
        pickle.dump((bin_centers, r_min, r_max, bgs_Mlyr_hist, lrg_Mlyr_hist, elg_Mlyr_hist, qso_Mlyr_hist), file)

def get_sky_coverage_fraction (tracer):
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

def plot_galaxy_distribution(bin_centers, histogram, label=None):

    plt.step(bin_centers, histogram, where='mid', linewidth = 1, label = label);
