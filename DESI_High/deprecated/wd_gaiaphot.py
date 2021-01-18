# Load the pure hydrogen atmosphere WD model (DA white dwarf)
# YouTubeVideo('sQbtjUgIVhg', width=800, height=400)
model_phot_da = np.genfromtxt('dat/WDphot/Table_DA.txt', skip_header=2)

# Eff. temperature.
T_da    = model_phot_da[:,0]

# Surface gravity.
logg_da = model_phot_da[:,1]

# Mass
mass_da = model_phot_da[:,2]

# absolute bolometric magnitude?
Mbol_da = model_phot_da[:,3]

# GAIA passbands. 
G_da    = model_phot_da[:,-6]
Gbp_da  = model_phot_da[:,-5]
Grp_da  = model_phot_da[:,-4]

# Interpolate between the models on the grid.
G_interpolator_da   = interpolate.CloughTocher2DInterpolator((T_da, logg_da), G_da,   rescale=True)
Gbp_interpolator_da = interpolate.CloughTocher2DInterpolator((T_da, logg_da), Gbp_da, rescale=True)
Grp_interpolator_da = interpolate.CloughTocher2DInterpolator((T_da, logg_da), Grp_da, rescale=True)

## -- Find the spectroscopic distance. 

# The expected absolute magnitude in the Gaia colours are
G_desi   = G_interpolator_da(T_desi, logg_desi)
Gbp_desi = Gbp_interpolator_da(T_desi, logg_desi)
Grp_desi = Grp_interpolator_da(T_desi, logg_desi)

# The distance modulus
dist_mod_G_desi   = G_obs   - G_desi
dist_mod_Gbp_desi = Gbp_obs - Gbp_desi
dist_mod_Grp_desi = Grp_obs - Grp_desi

# The medan distance modulus
dist_mod_mean_desi = np.mean((dist_mod_G_desi, dist_mod_Gbp_desi, dist_mod_Grp_desi))

# Converting the distance modulus back to distance
dist_spec = 10.**((dist_mod_mean_desi + 5.) / 5.) 

print('This Spectroscopic distance is %s pc' % dist_spec)

## -- Find the photometric distance.

# Define the least-square fucntion
def t_logg_fit(params, obs):
    # Unpack the input
    T, dist = params
    if (dist < 1.) or (T < 2500.) or (T > 1E+5):
        return np.inf
    logg, G_obs, Gbp_obs, Grp_obs, G_obs_err, Gbp_obs_err, Grp_obs_err = obs
    # Get the distance modulus
    dist_mod = 5. * np.log10(dist) - 5.
    # Get the model magnitude and apply the distance modulus
    G_model = G_interpolator_da(T, logg) + dist_mod
    Gbp_model = Gbp_interpolator_da(T, logg) + dist_mod
    Grp_model = Grp_interpolator_da(T, logg) + dist_mod
    # Check finite
    if not (np.isfinite(G_model) & np.isfinite(Gbp_model) & np.isfinite(Grp_model)):
        return np.inf
    else:
        return ((G_model - G_obs)/G_obs_err)**2. + ((Gbp_model - Gbp_obs)/Gbp_obs_err)**2. + ((Grp_model - Grp_obs)/Grp_obs_err)**2.

lsq_solution = optimize.minimize(
    t_logg_fit,
    [25000, 1000.],
    [logg_desi, G_obs, Gbp_obs, Grp_obs, G_obs_err, Gbp_obs_err, Grp_obs_err],
    method='Nelder-Mead')
T_phot = lsq_solution.x[0]
dist_phot = lsq_solution.x[1]

print('Photometric distance from Gaia photometry alone is %s pc' %dist_phot)
