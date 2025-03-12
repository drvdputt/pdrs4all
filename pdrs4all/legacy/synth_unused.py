"""Some parts of the original "synth.py" are used in the NIRSpec
postprocessing. Unused parts were moved here. We keep these archived in
case some of the more advanced steps by Ryan Chown are needed again."""

import numpy as np
from astropy.io import fits


def get_wave(fits_file):
    """
    Gets wavelength array from JWST fits file.
    Inputs:
    -------
    fits_file: 3D JWST Spectral cube with header.

    Outputs:
    -------
    wave_array : array with wavelenghts of the cube
    """
    # primary_hdr = fits.getheader(fits_file, "PRIMARY")
    science_hdr = fits.getheader(fits_file, "SCI")

    lambda0 = science_hdr["CRVAL3"] - science_hdr["CDELT3"] * (
        science_hdr["CRPIX3"] - 1
    )
    lambdas = np.arange(
        lambda0,
        lambda0 + (science_hdr["NAXIS3"] - 0.1) * science_hdr["CDELT3"],
        science_hdr["CDELT3"],
    )
    return lambdas
