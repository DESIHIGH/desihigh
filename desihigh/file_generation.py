# This file will contain all the functions created to generate the DESIHIGH files
# Some of those functions might have hardcoded paths or parameters.
import numpy as np

def gen_black_body(filename: str, save_to: str, shift: float = 0.3) -> None:
    """
    Generate a csv file containing the black body spectrum redshifted by a given value.
    Input file found at https://www.pveducation.org/pvcdrom/appendices/standard-solar-spectra

    Parameters
    ----------
    filename : str
        Path to the input file containing the black body spectrum.
    save_to : str
        Path to the output file where the redshifted spectrum will be saved.
    shift : float, optional
        Redshift value to apply to the spectrum, by default 0.3
    """
    black_body = np.genfromtxt(filename, delimiter=";", names=True)
    flux = black_body['Wm2nm1']
    wavelength = black_body['Wavelength_nm']
    shifted_wavelength = wavelength * np.sqrt(1+shift/(1-shift))
    
    mask = (wavelength < 4000) # Keep only wavelengths below 4000 nm
    wavelength = wavelength[mask]
    flux = flux[mask]
    shifted_wavelength = shifted_wavelength[mask]

    stack = np.column_stack((wavelength, shifted_wavelength, flux))
    headers = "wavelength_nm, shifted_wavelength_nm, flux_Wm2nm1"
    np.savetxt(save_to, stack, delimiter=",", header=headers)