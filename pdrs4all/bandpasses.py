from astropy import units as u
import numpy as np

PATH_TO_THROUGHPUT = "/Users/driesvandeputte/Projects/pdrs4all/sep3_stitching/nirspec_code_ryan/nircam_throughputs"


def read_nircam(use_v0_filters=False):
    """
    Parameters
    ---------
        use_v0_filters: if True, use old filters. If False, use filters that became
            available as of late-November 2022 (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters#NIRCamFilters-filt_trans)

    Returns
    -------

    dict in the format: d[<filter name>] = (reference wavelength,
    wavelength array, efficiency array)

    """
    path_to_throughputs = PATH_TO_THROUGHPUT
    nircam_bandpasses = {}

    # Short wavelength filters
    filters_sw = [
        "F070W",
        "F090W",
        "F115W",
        "F140M",
        "F150W",
        "F162M",
        "F164N",
        "F182M",
        "F187N",
        "F200W",
        "F210M",
        "F212N",
    ]
    # Long wavelength filters
    filters_lw = [
        "F250M",
        "F277W",
        "F300M",
        "F323N",
        "F335M",
        "F356W",
        "F360M",
        "F405N",
        "F410M",
        "F430M",
        "F444W",
        "F460M",
        "F466N",
        "F470N",
        "F480M",
    ]
    filters_all = filters_sw + filters_lw

    for filtername in filters_all:

        if use_v0_filters:
            wave = np.loadtxt(
                f"{path_to_throughputs}nircam_throughputs/modAB_mean/nrc_plus_ote/{filtername}_NRC_and_OTE_ModAB_mean.txt",
                skiprows=1,
            )[:, 0]
            eff = np.loadtxt(
                f"{path_to_throughputs}nircam_throughputs/modAB_mean/nrc_plus_ote/{filtername}_NRC_and_OTE_ModAB_mean.txt",
                skiprows=1,
            )[:, 1]
        else:
            wave, eff = np.loadtxt(
                f"{path_to_throughputs}nircam_throughputs/mean_throughputs/{filtername}_mean_system_throughput.txt",
                skiprows=1,
                unpack=True,
            )

        # compute the reference wave
        # defined as the pivot wavelength in Gordon et al. (2022)
        inttop = np.trapz(wave * eff, wave)
        intbot = np.trapz(eff / wave, wave)
        ref_wave = np.sqrt(inttop / intbot)

        nircam_bandpasses[filtername] = (ref_wave, wave * u.micron, eff)

    return nircam_bandpasses