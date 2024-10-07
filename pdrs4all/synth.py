from astropy.modeling.models import PowerLaw1D, BlackBody
import numpy as np
import matplotlib.pyplot as pl
import stitch_mrs
# import spectral_cube
from specutils.spectra import Spectrum1D
import astropy
from astroquery.svo_fps import SvoFps
from astropy import units as u
from jwst import datamodels
from reproject import reproject_exact
from astropy import wcs
from astropy.io import fits
import read_miri
import read_nircam
import pickle
import matplotlib
import os, sys

import astropy.visualization as astrovis
import matplotlib.colors as matcol

SEP4_GIT_DIR = '/Users/ryan/Documents/GitHub/sep4/'
path_to_data = '/Volumes/data1/jwst/sep4/'
DEFAULT_MIN_THROUGHPUT = 1.e-3 # 1e-4  # 5e-2
DEFAULT_MIN_COVERAGE = 0.95 # 1.


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
        science_hdr["CRPIX3"] - 1)
    lambdas = np.arange(
        lambda0,
        lambda0 + (science_hdr["NAXIS3"] - 0.1) * science_hdr["CDELT3"],
        science_hdr["CDELT3"])
    return lambdas


def get_2d_wcs_from_cube(fname):
    '''
    Gets 2D (spatial) WCS from IRS cube.
    For some reason, extracting WCS from cubism cubes doesn't work well
    (the spectral axis messes things up).

    '''
    fits_in = fits.open(fname)
    w_in = wcs.WCS(fits_in[1].header, fobj=fits_in, naxis=2)
    # Turn WCS into a header, and then back to WCS again (this cleans up garbage to do with the 3rd axis we don't want anyway)
    w_in = wcs.WCS(w_in.to_header())
    return w_in


def compute_colorcor(wave, bandpass, flux_ref, wave_ref, flux_source):
    """
    Compute the color correction K given the bandpass, reference spectrum,
    and source spectrum.  To use this color correction, divide the flux density
    for a band by K.  Such color corrections are needed to compute the correct
    flux density at the reference wavelength for a source with the flux_source
    spectral shape in the photometric convention that provides the flux density
    at a reference wavelength (convention B, see Gordon et al. 2022 for details).
    Parameters
    ----------
    wave : nd float array
       the wavelengths of the bandpass, flux_ref, and flux_source
    bandpass : nd float array
        end-to-end, total throughput bandpass of filter in fractional units
    flux_ref : nd float array
        reference flux density F(lambda) as a function of wave
    wave_ref : float
        reference wavelength
    flux_source : nd float array
        source flux density F(lambda) as a function of wave
    """
    # get the flux densities at the reference waveength
    flux_source_lambda_ref = np.interp(wave_ref, wave, flux_source)
    flux_ref_lambda_ref = np.interp(wave_ref, wave, flux_ref)

    # compute the top and bottom integrals
    inttop = np.trapz(wave * bandpass * flux_source / flux_source_lambda_ref,
                      wave)
    intbot = np.trapz(wave * bandpass * flux_ref / flux_ref_lambda_ref, wave)

    return inttop / intbot


def trapz_uncertainty(y_err, x):
    '''
    Uncertainty in np.trapz(y, x) using error propagation with y_err (uncertainty in y)
    '''
    dx = np.diff(x)
    res = 0.5 * np.sqrt(np.dot(dx**2, (y_err[1:]**2 + y_err[:-1]**2)))
    return res


def synthetic_image_jwst(waves,
                         spectrum,
                         bandpasses,
                         filter_key,
                         spectrum_unc=None,
                         min_throughput=DEFAULT_MIN_THROUGHPUT,
                         min_coverage=1.):
    '''
    Generic function to calculate JWST synthetic images.
    Will calculate uncertainties if spectrum_unc is provided.
    '''
    ref_shape = PowerLaw1D(amplitude=1.0, x_0=1.0, alpha=0.0)
    flux_ref = ref_shape(waves)

    print(filter_key)
    # 1. Effective wavelength
    # 2. Wavelength grid
    # 3. Total throughput of instrument over wavelength grid
    lambda_eff, lambdas_filt, tp_filt = bandpasses[filter_key]

    lambdas_filt = lambdas_filt.to(u.micron).value

    # Number of wavelength bins where throughput is > threshold
    good_filt = (tp_filt / tp_filt.max()) >= min_throughput
    lambdas_filt = lambdas_filt[good_filt]
    tp_filt = tp_filt[good_filt]

    lambda_min = lambdas_filt.min()
    lambda_max = lambdas_filt.max()

    # Calculate wavelength coverage (we don't want to consider spaxels without full coverage)
    good_unc = np.full(waves.size, fill_value=True, dtype=bool)
    if spectrum_unc is not None:
        good_unc = np.isfinite(spectrum_unc) & (spectrum_unc != 0.)

    all_spec = (waves >= lambda_min) & (waves <= lambda_max)
    good_spec = all_spec & good_unc
    waves_i = waves[good_spec]
    flux_source_i = spectrum[good_spec]
    flux_ref_i = flux_ref[good_spec]

    if spectrum_unc is not None:
        unc_flux_source_i = spectrum_unc[good_spec]

    # Total number of wavelengths over the bandpass
    n_waves = np.sum(all_spec)
    coverage = np.sum(spectrum[good_spec] != 0)
    if (coverage >= n_waves * min_coverage) and (n_waves >= 2):
        print(coverage, n_waves)
        print(f"{lambda_min} um to {lambda_max} um")
        # Interpolate throughput on wavelength grid of stitched cube
        tp_interp = np.interp(waves_i, lambdas_filt, tp_filt)

        colour_corr = compute_colorcor(waves_i, tp_interp, flux_ref_i,
                                       lambda_eff, flux_source_i)

        flux_a_numer = np.trapz(waves_i * flux_source_i * tp_interp, waves_i)
        flux_a_denom = np.trapz(waves_i * tp_interp, waves_i)

        flux_convention_a = flux_a_numer / flux_a_denom
        flux_convention_b = np.interp(lambda_eff, waves_i, flux_source_i)

        # Deal with uncertainties if needed
        if spectrum_unc is not None:
            ### Convention A
            # Calculate uncertainty in numerator and denominator
            unc_numer = trapz_uncertainty(
                waves_i * unc_flux_source_i * tp_interp, waves_i)
            unc_denom = 0.
            # Calculate uncertainty in numerator/denominator
            unc_flux_convention_a = np.abs(flux_convention_a) * np.sqrt(
                (unc_numer / flux_a_numer)**2 + (unc_denom / flux_a_denom)**2)

            ### Convention B
            unc_flux_convention_b = np.interp(lambda_eff, waves_i,
                                              unc_flux_source_i)

            return flux_convention_a, flux_convention_b, colour_corr, unc_flux_convention_a, unc_flux_convention_b
        else:
            return flux_convention_a, flux_convention_b, colour_corr
    else:
        if spectrum_unc is not None:
            return np.nan, np.nan, np.nan, np.nan, np.nan
        else:
            return np.nan, np.nan, np.nan


def make_synthetic_images_from_cube(comb_cube,
                                    waves,
                                    miri=True,
                                    nircam=False,
                                    filter_names=[],
                                    unc_comb_cube=None,
                                    min_throughput=DEFAULT_MIN_THROUGHPUT,
                                    min_coverage=DEFAULT_MIN_COVERAGE,
                                    path_to_throughputs=''):

    nw, nx, ny = comb_cube.shape

    ref_shape = PowerLaw1D(amplitude=1.0, x_0=1.0, alpha=0.0)
    flux_ref = ref_shape(waves)

    if miri:
        miri_nircam = 'MIRIm'
        mrs_nirspec = 'mrs'
        bandpasses = read_miri.read_miri()
        if len(filter_names) == 0:
            # Filters used in ERS 1288
            # filter_names = ['F1130W', 'F1500W', 'F2550W', 'F770W']
            filter_names = ['F1000W', 'F2100W', 'F1280W', 'F560W', 'F1130W', 'F1500W', 'F2550W', 'F770W']

        filter_names = [f.lower() for f in filter_names]
    elif nircam:
        miri_nircam = 'NIRCam'
        mrs_nirspec = 'nirspec'
        bandpasses = read_nircam.read_nircam(path_to_throughputs)
        if len(filter_names) == 0:
            # Filters used in ERS 1288
            filter_names = [
                'F140M', 'F182M', 'F187N', 'F210M', 'F212N', 'F277W', 'F300M',
                'F335M', 'F480M', 'F162M', 'F164N', 'F323N', 'F405N', 'F470N'
            ]

    print(f"Making synthetic {miri_nircam} images from {mrs_nirspec} cube.")
    print(f"Filters: {filter_names}")

    n_filt = len(filter_names)
    cc_dict = dict()
    synth_image_dict = dict()
    synth_image_dict['convention_a'] = dict()
    synth_image_dict['convention_b'] = dict()
    if unc_comb_cube is not None:
        synth_image_dict['unc_convention_a'] = dict()
        synth_image_dict['unc_convention_b'] = dict()

    for i in range(n_filt):
        cc_dict[filter_names[i]] = np.full((nx, ny), fill_value=np.nan)
        synth_image_dict['convention_a'][filter_names[i]] = np.full(
            (nx, ny), fill_value=np.nan)
        synth_image_dict['convention_b'][filter_names[i]] = np.full(
            (nx, ny), fill_value=np.nan)
        if unc_comb_cube is not None:
            synth_image_dict['unc_convention_a'][filter_names[i]] = np.full(
                (nx, ny), fill_value=np.nan)
            synth_image_dict['unc_convention_b'][filter_names[i]] = np.full(
                (nx, ny), fill_value=np.nan)

    for ix in range(nx):
        for iy in range(ny):
            print(f"ix, iy = {ix, iy}")
            flux_source = comb_cube[:, ix, iy]
            if unc_comb_cube is not None:
                unc_flux_source = unc_comb_cube[:, ix, iy]

            pwaves = []
            pcolor = []
            for filter_key in filter_names:
                min_throughput_temp = min_throughput
                min_coverage_temp = min_coverage

                if filter_key in ['F140M', 'F277M']:
                    min_throughput_temp = 1e-3
                    min_coverage_temp = 0.8

                if unc_comb_cube is not None:
                    flux_convention_a, flux_convention_b, colour_corr, unc_flux_convention_a, unc_flux_convention_b = synthetic_image_jwst(
                        waves,
                        flux_source,
                        bandpasses,
                        filter_key,
                        spectrum_unc=unc_flux_source,
                        min_throughput=min_throughput_temp,
                        min_coverage=min_coverage_temp)
                else:
                    flux_convention_a, flux_convention_b, colour_corr = synthetic_image_jwst(
                        waves,
                        flux_source,
                        bandpasses,
                        filter_key,
                        min_throughput=min_throughput_temp,
                        min_coverage=min_coverage_temp)

                cc_dict[filter_key][ix, iy] = colour_corr
                synth_image_dict['convention_a'][filter_key][
                    ix, iy] = flux_convention_a
                synth_image_dict['convention_b'][filter_key][
                    ix, iy] = flux_convention_b

                if unc_comb_cube is not None:
                    synth_image_dict['unc_convention_a'][filter_key][
                        ix, iy] = unc_flux_convention_a
                    synth_image_dict['unc_convention_b'][filter_key][
                        ix, iy] = unc_flux_convention_b

    #
    # # Pickle dictionary
    # fname_pickle = f'{output_path}{mrs_nirspec}_synthetic_images_both_conventions{pointing_group}.pk'
    # if nircam:
    #     fname_pickle = f'{output_path}synthetic_nircam_images_both_conventions_{pointing_group}.pk'
    # print(f"Writing output 1/2 to {fname_pickle}")
    # with open(fname_pickle, 'wb') as fl:
    #     pickle.dump(synth_image_dict, fl, pickle.HIGHEST_PROTOCOL)
    #
    # fname_pickle = f'{output_path}{mrs_nirspec}_synthetic_images_colour_corrections{pointing_group}.pk'
    # if nircam:
    #     fname_pickle = f'{output_path}synthetic_nircam_images_colour_corrections_{pointing_group}.pk'
    # print(f"Writing output 2/2 to {fname_pickle}")
    # with open(fname_pickle, 'wb') as fl:
    #     pickle.dump(cc_dict, fl, pickle.HIGHEST_PROTOCOL)

    return synth_image_dict, cc_dict


def write_synth_to_fits(nircam_synth_images, w, fname_out_synth,
                        fname_out_synth_unc):
    # New HDUList for synthetic images
    # Every fits file needs a PrimaryHDU. We'll make a blank one
    hdu0 = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list = [hdu0]
    # Same for uncertainties
    hdu0_unc = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list_unc = [hdu0_unc]

    filter_list = list(nircam_synth_images['convention_a'].keys())
    hdr0 = w.to_header()

    for filt_tmp in filter_list:
        arr = nircam_synth_images['convention_a'][filt_tmp]
        unc_arr = nircam_synth_images['unc_convention_a'][filt_tmp]
        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr['EXTNAME'] = (filt_tmp, 'Filter')
        q = 1 * u.MJy / u.sr
        hdr['UNIT'] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=arr, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr['EXTNAME'] = (filt_tmp, 'Filter')
        q = 1 * u.MJy / u.sr
        hdr['UNIT'] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=unc_arr, header=hdr)
        hdu_list_unc = hdu_list_unc + [hdu_i]

    hdul = fits.HDUList(hdus=hdu_list)
    hdul.writeto(fname_out_synth, output_verify='fix', overwrite=True)

    hdul_unc = fits.HDUList(hdus=hdu_list_unc)
    hdul_unc.writeto(fname_out_synth_unc, output_verify='fix', overwrite=True)


class Synth:
    """
    """

    def __init__(self,
                 path_to_throughputs='',
                 synth_image_dict=None,
                 cc_dict=None):
        """
        """
        self.path_to_throughputs = path_to_throughputs
        self.synth_image_dict = synth_image_dict
        self.cc_dict = cc_dict

    def make_save_miri_bandpasses(self):
        bandpasses = read_miri.read_miri()
        # Pickle dictionary
        fname_pickle = f'{path_to_data}MIRI_bandpasses.pk'
        with open(fname_pickle, 'wb') as fl:
            pickle.dump(bandpasses, fl, pickle.HIGHEST_PROTOCOL)

    def load_miri_bandpasses(self):
        fname_pickle = f'{path_to_data}MIRI_bandpasses.pk'
        with open(fname_pickle, 'rb') as fl:
            bandpasses = pickle.load(fl)
        return bandpasses

    def images(self, fname_stitched_cube, miri=False, nircam=True):
        comb_cube = fits.open(fname_stitched_cube)['CUBE'].data
        comb_cube_unc = fits.open(fname_stitched_cube)['ERR'].data
        waves = fits.open(fname_stitched_cube)['WAVE'].data

        # comb_cube, w = stitch_mrs.run_ers1288_nirspec(pointing_group=pointing_group, concat=False, cross_cut=False, force_zero_scale_factor=True)
        # waves = comb_cube[0, 0, 0, :]
        # comb_cube = comb_cube[:, :, 1, :].T
        if nircam:
            bandpasses = read_nircam.read_nircam(self.path_to_throughputs)
        if miri:
            bandpasses = read_miri.read_miri()

        synth_image_dict, cc_dict = make_synthetic_images_from_cube(
            comb_cube,
            waves,
            miri=miri,
            nircam=nircam,
            unc_comb_cube=comb_cube_unc,
            min_throughput=DEFAULT_MIN_THROUGHPUT,
            min_coverage=DEFAULT_MIN_COVERAGE,
            path_to_throughputs=self.path_to_throughputs)

        self.synth_image_dict = synth_image_dict
        self.cc_dict = cc_dict

        return synth_image_dict, cc_dict
