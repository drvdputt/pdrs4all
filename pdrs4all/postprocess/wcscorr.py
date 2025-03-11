from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.background import MedianBackground, Background2D
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.segmentation import detect_sources, SourceCatalog, deblend_sources

# hardcoded reference coordinate for the proplyd
C_PROPLYD_NIRSPEC = SkyCoord(
    83.83447199597 * u.deg, -5.41778904231 * u.deg, frame="icrs"
)
C_OTHER_PROPLYD_MRS = SkyCoord(
    83.83447199597 * u.deg, -5.41778904231 * u.deg, frame="icrs"
)


def apply_delta_ra_dec_to_wcs(current_wcs, delta_ra_dec):
    """Create modified copy of WCS with given shift in RA and Dec applied.

    Synergizes with the output of e.g.
    determine_offset_using_photutils_and_reference_coordinate

    Parameters
    ----------

    current_wcs: WCS

    delta_ra_dec: (quantity, quantity)
        The output of one of the functions that determines the RA and Dec
        offset.

    """
    w_new = current_wcs.deepcopy()
    dra, ddec = delta_ra_dec
    w_new.wcs.crval = w_new.wcs.crval + np.array([dra.value, ddec.value])
    return w_new


def delta_ra_dec_photometry_vs_ref_coord(c_true, reference_image, current_wcs):
    """Slightly more general. See nirspec version for explanation.

    Returns: dra, ddec
        Spherical offsets in RA and Dec units (check unit of the
        returned astropy quantity)"""
    npixels = 5

    # photutils source identification
    bkg_estimator = MedianBackground()
    xd, yd = reference_image.shape
    bkg = Background2D(
        reference_image, (10, 10), filter_size=(3, 3), bkg_estimator=bkg_estimator
    )
    threshold = 3.0 * bkg.background_rms
    fwhm = 2
    sigma = fwhm * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma)  # , xsize=1, ysize=1)
    kernel.normalize(mode="integral")
    convolved_data = convolve(reference_image - bkg.background, kernel)
    segm = detect_sources(convolved_data, threshold, npixels=npixels)
    segm_deblend = deblend_sources(
        convolved_data,
        segm,
        npixels=npixels,
        nlevels=32,
        contrast=0.001,
    )
    cat = SourceCatalog(reference_image, segm_deblend)
    tbl = cat.to_table()

    print("Found the following sources using photutils")
    print(tbl)

    good_tbl_rows = np.isfinite(tbl["xcentroid"]) & np.isfinite(tbl["ycentroid"])
    tbl = tbl[good_tbl_rows]

    n_apertures = len(tbl)
    sep = np.zeros(n_apertures)
    xc = tbl["xcentroid"]
    yc = tbl["ycentroid"]

    c_measured = SkyCoord.from_pixel(xc, yc, wcs=current_wcs, origin=0)
    sep = c_measured.separation(c_true).deg

    # Pixel indices of closest peak
    i_min = np.argmin(sep)
    xmin, ymin = xc[i_min], yc[i_min]
    sepmin = sep[i_min]
    print("Distance of measured peak to true coordinate.")
    print(f"xmin, ymin, sep = {xmin:.4f}, {ymin:.4f}, {sepmin * 3600.:.4f} arcsec")
    # Calculate offsets in RA and DEC needed to bring measured
    # center to the true center
    dra, ddec = c_measured[i_min].spherical_offsets_to(c_true)
    print(dra, ddec)
    return dra, ddec


def nirspec_wcscorr_using_proplyd(reference_image, current_wcs):
    """Use hardcoded proplyd position for WCS correction.

    Prerequisite: a slice in which the proplyd can be detected, or
    synthetic photometry.

    Technically we could correct the filters separately, but at the
    moment, this is used only for the "naive stitched" cube. But I
    expect the WCS differences between the pointings to be larger,
    rather than those between the filters. But the individual
    pointings don't contain any compact sources, so those cannot be
    corrected this way.

    1. Photutils source detection
       - estimate and subtract background around proplyd
       - detect and deblend segmentation map
       - make source catalog
       - clean catalog and get centroids (xy)

    2. Convert xy to RA, Dec, compute separations, compute RA and
    Dec offset for minimum.

    3. Copy WCS and add offsets to CRVAL

    Parameters
    ----------

    reference_image: 2d array
        data from which to derive correction. Needs to peak at the
        proplyd; F212N is used in Ryans approach.

    current_cwcs: WCS
        celestial wcs to be corrected

    Returns
    -------

    new_wcs: WCS
        This wcs should be applied to the stitched cube and the
        synthetic images.

    """
    dra, ddec = delta_ra_dec_photometry_vs_ref_coord(
        reference_image, C_PROPLYD_NIRSPEC, current_wcs
    )
    return apply_delta_ra_dec_to_wcs(current_wcs, dra, ddec)


def mrs_wcscorr_using_proplyd(images, current_wcss):
    """WCS correction method for the 12 MRS cubes.

    The special thing about this method is, that an offset is computed
    using reference images for A, B, and C (preferably channel 1) and
    then the same offsets are applied to channels 2,3,4, to A, B, and C
    respectively.

    Parameters
    ----------

    images: list of 2D array
        Images to perform photometry on to determine offset. Only the
        first 3 will be used.

    current_wcss: list of WCS
        The 12 WCS of the MRS cubes, in wavelenght order (ch1 A, B, C,
        ... ch4 A, B, C).

    """
    ref_images = images[:3]
    ref_wcss = current_wcss[:3]

    # calculate offsets for A, B, and C
    delta_ra_dec_abc = [
        delta_ra_dec_photometry_vs_ref_coord(C_OTHER_PROPLYD_MRS, ref_im, wcs)
        for ref_im, wcs in zip(ref_images, ref_wcss)
    ]

    # to the 12 bands, apply the offsets as 0 1 2, 0 1 2, ...
    new_wcss = []
    for ch in range(4):
        for abc in range(3):
            new_wcs = apply_delta_ra_dec_to_wcs(
                current_wcss[ch * 3 + abc], delta_ra_dec_abc[abc]
            )
            new_wcss.append(new_wcs)

    return new_wcss
