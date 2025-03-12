import matplotlib.pyplot as pl
import numpy as np

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u
from astroquery.gaia import Gaia
# Select early Data Release 3
Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"

from photutils import detect_sources
from photutils import SourceCatalog
from photutils import deblend_sources
from photutils import EllipticalAperture
from photutils import Background2D
from photutils import MedianBackground
from photutils.background import MMMBackground

import pickle
import compress_pickle
from reproject import reproject_exact

from scipy.stats import linregress
from scipy import stats


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


def run_linmix(x, y, xerr, yerr, parallelize=False, nchains=4, delta=None):
    from linmix import linmix
    pivot = np.mean(x)
    pivot_y = np.mean(y)
    lm_result = linmix.LinMix(x - pivot,
                              y - pivot_y,
                              xerr,
                              yerr,
                              K=3,
                              parallelize=parallelize,
                              nchains=nchains,
                              delta=delta)
    lm_result.run_mcmc(silent=True)
    alphas = lm_result.chain['alpha']
    betas = lm_result.chain['beta']
    chains = np.vstack([alphas, betas])
    result = dict()
    result['chains'] = chains
    result['intercept'] = np.average(chains[0])
    result['intercept_err'] = np.std(chains[0])
    result['slope'] = np.average(chains[1])
    result['slope_err'] = np.std(chains[1])
    result['pivot'] = pivot
    result['pivot_y'] = pivot_y
    return result


def plot_fit(t,
             a,
             b,
             a_err=0,
             b_err=0,
             inv=False,
             xin=None,
             yin=None,
             yin_err=None,
             s=None,
             pivot=0,
             pivot_y=0,
             ax=None,
             log=False,
             color='b',
             lw=2,
             alpha=0.5,
             text=True,
             fontsize=9,
             **kwargs):
    """
    alpha is used to shade the uncertainties from a_err and b_err
    **kwargs is passed to plt.plot() for the central line only
    the error band has zorder=-10
    """
    if log:
        if pivot == 0:
            pivot = 1
        y = lambda A, B: 10**A * (t / pivot)**B
    else:
        if inv:
            y = lambda A, B: A + B * t + pivot
        else:
            y = lambda A, B: A + B * (t - pivot) + pivot_y
    if ax is None:
        ax = pl
    # the length may vary depending on whether it's a default color
    # (e.g., 'r' or 'orange') or an rgb(a) color, etc, but as far as
    # I can think of none of these would have length 2.
    if len(color) != 2:
        color = (color, color)
    print('in lnr.plot: color =', color)
    ax.plot(t, y(a, b), ls='-', color=color[0], lw=lw, **kwargs)
    if a_err != 0 or b_err != 0:
        # to make it compatible with either one or two values
        a_err = np.array([a_err]).flatten()
        b_err = np.array([b_err]).flatten()
        if a_err.size == 1:
            a_err = [a_err, a_err]
        if b_err.size == 1:
            b_err = [b_err, b_err]
        err = [
            y(a - a_err[0], b - b_err[0]),
            y(a - a_err[0], b + b_err[1]),
            y(a + a_err[1], b - b_err[0]),
            y(a + a_err[1], b + b_err[1])
        ]
        ylo = np.min(err, axis=0)
        yhi = np.max(err, axis=0)
        ax.fill_between(t,
                        ylo,
                        yhi,
                        color=color[1],
                        alpha=alpha,
                        lw=0,
                        edgecolor='none',
                        zorder=-10)
    if s:
        if log:
            ax.plot(t, (1 + s) * y(a, b), ls='--', color=color[0], lw=lw)
            ax.plot(t, y(a, b) / (1 + s), ls='--', color=color[0], lw=lw)
        else:
            ax.plot(t, y(a, b) + s, ls='--', color=color[0], lw=lw)
            ax.plot(t, y(a, b) - s, ls='--', color=color[0], lw=lw)
    ax = pl.gca()
    print(a, a_err, b, b_err)
    if text:
        if pivot != 0:
            # string = '$y = a + b(x-x_0)$\n'
            string = '$y = a + bx$\n'
            # string += '$x_0 = %.3f$\n' % (pivot, )
            if inv:
                string += '$a=%.2f \pm %.2f$\n' % (a + pivot, a_err[0])
            else:
                string += '$a=%.2f \pm %.2f$\n' % (
                    a - b * pivot + pivot_y,
                    np.sqrt(a_err[0]**2 + (pivot * b_err[0])**2))
            string += '$b=%.2f \pm %.2f$\n' % (b, b_err[0])
        else:
            string = '$y = a + bx$\n'
            string += '$a=%.2f \pm %.2f$\n' % (a, a_err[0])
            string += '$b=%.2f \pm %.2f$\n' % (b, b_err[0])

        if b > 0:
            xt, yt = 0.05, 0.95
        else:
            xt, yt = 0.05, 0.4
        ax.text(xt,
                yt,
                string,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes,
                fontsize=fontsize)

        # txt = '${\\rm Spearman/Pearson}$\n'
        # r, p = stats.spearmanr(xin, yin)
        # txt += '$r=%.2g\, p=%.2g$\n' % (r, p)
        # r, p = stats.pearsonr(xin, yin)
        # txt += '$r=%.2g\, p=%.2g$\n' % (r, p)
        # txt = '${\\rm Pearson}$\n'
        r, p = stats.pearsonr(xin, yin)
        if yin_err is not None:
            # Calculate intrinsic scatter
            good_temp = (xin != 0) & (yin != 0) & (yin_err != 0) & np.isfinite(
                xin) & np.isfinite(yin) & np.isfinite(yin_err)
            sig_tot = np.sqrt(1. / (xin[good_temp].size - 2) * np.sum(
                (yin[good_temp] -
                 (a - b * pivot + pivot_y + b * xin[good_temp]))**2))
            sig_int = np.sqrt(sig_tot**2 - 1. / xin[good_temp].size *
                              np.sum(yin_err[good_temp]**2))
            txt = 'Pearson $r=%.2g$\n' % (r, )
            txt += '$\sigma_\mathrm{int}=%.2g$ dex' % (sig_int, )
        else:
            txt = 'Pearson $r=%.2g$' % (r, )

        if b > 0:
            xt, yt = 0.95, 0.2
        else:
            xt, yt = 0.95, 0.95
        ax.text(xt,
                yt,
                txt,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes,
                fontsize=fontsize)

    return


def get_radec_offsets(filename,
                      c_true,
                      make_plot=False,
                      npixels=5,
                      fwhm_pix=3.,
                      approx_radius_pix=2.):
    # cal_name : (str) name of calibration object e.g. '3C273'
    # filename : (str) full path to reduced calibration data fits file

    # ff = fits.open('/home/chownrj/18B_test/CRL618/850/DR1/CRL618-DR1-fix_cal.fits')
    ff = fits.open(filename)
    data = ff[0].data[0]
    var_data = ff[1].data[0]

    # bkg_estimator = MedianBackground()
    # bkg = Background2D(data, (50, 50), filter_size=(3, 3),
    #                    bkg_estimator=bkg_estimator)
    # threshold = bkg.background + (2. * bkg.background_rms)
    threshold = 2. * np.sqrt(var_data)

    sigma = fwhm_pix * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data,
                          threshold,
                          npixels=npixels,
                          filter_kernel=kernel)
    segm_deblend = deblend_sources(data,
                                   segm,
                                   npixels=npixels,
                                   filter_kernel=kernel,
                                   nlevels=32,
                                   contrast=0.001)

    cat = SourceCatalog(data, segm_deblend)
    tbl = cat.to_table()
    tbl['xcentroid'].info.format = '.2f'  # optional format
    tbl['ycentroid'].info.format = '.2f'
    tbl['cxx'].info.format = '.2f'
    tbl['cxy'].info.format = '.2f'
    tbl['cyy'].info.format = '.2f'
    print(tbl)

    # Make apertures
    # cat = SourceCatalog(data, segm_deblend)
    # approx_radius_pix = 2. # approximate isophotal extent
    apertures = []
    for obj in cat:
        position = (obj.xcentroid.value, obj.ycentroid.value)
        a = obj.semimajor_axis_sigma.value * approx_radius_pix
        b = obj.semiminor_axis_sigma.value * approx_radius_pix
        theta = obj.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))

    if make_plot:
        # Plot apertures on image
        # norm = ImageNormalize(stretch=SqrtStretch())
        # fig, (ax1, ax2) = pl.subplots(2, 1, figsize=(10, 12.5))
        fig, ax1 = pl.subplots(1, 1, figsize=(10, 10))
        # ax1.imshow(data, origin='lower', cmap='Greys_r', vmin=0, vmax=5, interpolation='none', norm=norm)
        ax1.imshow(data,
                   origin='lower',
                   cmap='Greys_r',
                   vmin=0,
                   vmax=5,
                   interpolation='none')
        # ax1.set_title('Data')
        # ax2.imshow(segm_deblend, origin='lower',
        #            cmap=segm_deblend.cmap(random_state=12345))
        # ax2.set_title('Segmentation Image')
        for aperture in apertures:
            aperture.plot(color='white', lw=1.5, ax=ax1)
            # aperture.plot(color='white', lw=1.5, ax=ax2)

    # Get WCS from reduced data
    # w = wcs.WCS(ff[0].header)
    # w = w.dropaxis(2)  # get rid of third axis
    w = get_2d_wcs_from_cube(filename)

    n_apertures = len(tbl)
    sep = np.zeros(n_apertures)
    dra = np.zeros(n_apertures)
    ddec = np.zeros(n_apertures)
    xc = np.zeros(n_apertures)
    yc = np.zeros(n_apertures)

    for i in range(0, n_apertures):
        xc[i] = tbl[i]['xcentroid'].value
        yc[i] = tbl[i]['ycentroid'].value

        # CRL618 true coordinates:
        # 04 42 53.6245215366 +36 06 53.397219192
        # c_true = SkyCoord("4h42m53.6245215366s +36d06m53.397219192s", ICRS, equinox="J2000")
        # c_true = true_coords[cal_name]

        # Sky coords of center of detected object
        c_measured = SkyCoord.from_pixel(xc[i], yc[i], wcs=w, origin=0).icrs

        # Calculate separation between measured center and true center
        sep[i] = c_measured.separation(c_true).arcsec
        # Calculate offsets in RA and DEC needed to bring measured center to the true center
        tmp = c_measured.spherical_offsets_to(c_true)
        dra[i], ddec[i] = tmp[0].arcsec, tmp[1].arcsec

    dra = dra[np.abs(sep) == np.min(np.abs(sep))]
    ddec = ddec[np.abs(sep) == np.min(np.abs(sep))]
    xc = xc[np.abs(sep) == np.min(np.abs(sep))]
    yc = yc[np.abs(sep) == np.min(np.abs(sep))]
    sep = sep[np.abs(sep) == np.min(np.abs(sep))]

    if make_plot:
        figure_name = filename.split('/')[-1].split('-DR1')[0]
        ax1.set_title('Data; offset = %3.3f arcsec (%3.3f RA, %3.3f dec)' %
                      (sep, dra, ddec))

        pl.savefig('/home/chownrj/SCUBA2-public-modified/offset_figures/' +
                   figure_name + '.pdf')

    return sep, dra, ddec


# Get center of FOV
# centers = np.zeros(2)
# pointing_group = '0_1_2'
# fwhm_pix = 3.
# npixels=5
# fwhm_pix=3.
# approx_radius_pix=2.
#
# output_path = '/Volumes/data1/jwst/data/ers1288/nirspec/NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/synthetic_nircam_images/'
# fname_synth_images = f'{output_path}synthetic_nircam_images_both_conventions_pointing_group_{pointing_group}.pk'
# filt = 'F470N'
#
# path_to_data = '/Volumes/data1/jwst/data/ers1288/nirspec/'
# fname_cube = f"{path_to_data}NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/stitched_cubes/NIRSpec_stitched_allchannels_pointing_group_{pointing_group}.fits"
# fname_out_cube = f"{path_to_data}NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/stitched_cubes/NIRSpec_stitched_allchannels_pointing_group_{pointing_group}_wcsfixed.fits"
# fname_out_synth = '/Volumes/data1/jwst/data/ers1288/nirspec/NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_group_0_1_2_wcsfixed.fits'

# In [44]: dra
# Out[44]: array([-0.00059043])
#
# In [45]: ddec
# Out[45]: array([4.98643004e-05])
do_pointing = False
outlier_off = True
if do_pointing:

    pointing_group = '5'
    fwhm_pix = 3.
    npixels = 5
    fwhm_pix = 3.
    approx_radius_pix = 2.

    output_path = '/Volumes/data1/jwst/data/ers1288/nirspec/NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/synthetic_nircam_images/'
    fname_synth_images = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}.pk'
    # filt = 'F470N'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F470N-F444W_i2d.fits'
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_CLEAR-F212N_i2d.fits'
    filt = 'F212N'
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F162M-F150W2_i2d.fits'
    filt = 'F162M'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F405N-F444W_i2d.fits'
    # filt = 'F405N'
    # 'Level3_F162M-F150W2_i2d.fits'

    path_to_data = '/Volumes/data1/jwst/data/ers1288/nirspec/'
    if outlier_off:
        fname_cube = f"/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/NIRSpec_stitched_allchannels_pointing_{pointing_group}.fits"
    else:
        fname_cube = f"/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/from_arcade/NIRSpec_stitched_allchannels_pointing_{pointing_group}.fits.gz"
    fname_out_cube = f"/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/NIRSpec_stitched_allchannels_pointing_{pointing_group}_wcsfixed.fits"
    fname_out_synth = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}_wcsfixed.fits'
    fname_out_synth_unc = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}_unc_wcsfixed.fits'

else:
    # pointing_group = '5'
    fwhm_pix = 3.
    npixels = 5
    fwhm_pix = 3.
    approx_radius_pix = 4.
    data_dir = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/'
    if outlier_off:
        # Fix this
        'NIRSpec_stitched_allchannels_all_pointings_wcsfixed_outlier_off'
    fname_cube = f'{data_dir}NIRSpec_stitched_allchannels_all_pointings.fits.gz'
    fname_out_cube = f'{data_dir}NIRSpec_stitched_allchannels_all_pointings_wcsfixed.fits'
    filt = 'F470N'

    fname_synth_images = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_all_pointings.pk'
    fname_out_synth = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_all_pointings_wcsfixed.fits'
    fname_out_synth_unc = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_all_pointings_unc_wcsfixed.fits'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F470N-F444W_i2d.fits'
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/OrionBar_NIRCAM_Final_mocaics_Boris_aligned_corrected_10Oct22/Level3_F470N-F444W-A_i2d.fits'

    # # ---
    # output_path = '/Volumes/data1/jwst/data/ers1288/nirspec/NIRSpec_all_filters_Orion_Bar_from_MAST_step2_15Sep22/synthetic_nircam_images/'
    # fname_synth_images = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}.pk'
    # # filt = 'F470N'
    # # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F470N-F444W_i2d.fits'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_CLEAR-F212N_i2d.fits'
    # filt = 'F212N'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F162M-F150W2_i2d.fits'
    # filt = 'F162M'
    # # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F405N-F444W_i2d.fits'
    # # filt = 'F405N'
    # # 'Level3_F162M-F150W2_i2d.fits'
    #
    # path_to_data = '/Volumes/data1/jwst/data/ers1288/nirspec/'
    # fname_cube = f"/Vo    lumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/NIRSpec_stitched_allchannels_pointing_{pointing_group}.fits"
    # fname_out_cube = f"/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/NIRSpec_stitched_allchannels_pointing_{pointing_group}_wcsfixed.fits"
    # fname_out_synth = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}_wcsfixed.fits'
    # fname_out_synth_unc = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/synthetic_nircam_images_both_conventions_pointing_{pointing_group}_unc_wcsfixed.fits'
    #

filt = 'F470N'


def realign_nircam(fname_nircam, w, c_true):
    # nircam_synth_images = pickle.load(open(fname_synth_images, 'rb'))
    # w = get_2d_wcs_from_cube(fname_cube)
    footprint = w.calc_footprint(fits.open(fname_nircam)[1].header)
    # ras, decs = np.meshgrid(ras, decs)
    ramin, ramax = np.min(footprint[:, 0]), np.max(footprint[:, 0])
    decmin, decmax = np.min(footprint[:, 1]), np.max(footprint[:, 1])
    ra0 = (ramin + ramax) / 2.
    dec0 = (decmin + decmax) / 2.
    search_radius = max([np.abs(ramax - ramin), np.abs(decmax - decmin)]) * 2

    coord_search = SkyCoord(ra=ra0,
                            dec=dec0,
                            unit=(u.degree, u.degree),
                            frame='icrs')
    width = u.Quantity(np.abs(ramax - ramin), u.deg)
    height = u.Quantity(np.abs(decmax - decmin), u.deg)
    # cat_gaia = Gaia.query_object_async(coordinate=coord_search,
    #                                    width=width,
    #                                                height=height)
    # cat_gaia.pprint(max_lines=12, max_width=130)

    # c_true = SkyCoord(gaia_row['ra'] * u.deg,
    #                   gaia_row['dec'] * u.deg,
    #                   frame='icrs')

    data = fits.open(fname_nircam)['SCI'].data
    data_unc = fits.open(fname_nircam)['ERR'].data
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (10, 10),
                       filter_size=(3, 3),
                       bkg_estimator=bkg_estimator)
    # threshold = 100
    threshold = bkg.background + (3. * bkg.background_rms)
    print(f"background = {bkg.background}")
    # threshold = bkg.background + (2. * nircam_synth_images['unc_convention_a'][filt])
    # threshold = bkg.background + 20 * nircam_synth_images['unc_convention_a'][filt]
    sigma = 2 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=20)  #, kernel=kernel)
    segm_deblend = deblend_sources(data,
                                   segm,
                                   npixels=npixels,
                                   kernel=kernel,
                                   nlevels=32,
                                   contrast=0.001)

    cat = SourceCatalog(data, segm_deblend)
    tbl = cat.to_table()
    print(tbl)
    apertures = []
    for obj in tbl:
        position = (obj['xcentroid'], obj['ycentroid'])
        a = obj['semimajor_sigma'] * approx_radius_pix
        b = obj['semiminor_sigma'] * approx_radius_pix
        theta = obj['orientation']
        apertures.append(
            EllipticalAperture(position, a.value, b.value, theta=theta))

    fig, ax1 = pl.subplots(1, 1, figsize=(10, 10))
    ax1.imshow(data, origin='lower', cmap='Greys_r', interpolation='none')
    for aperture in apertures:
        aperture.plot(color='b', lw=1.5)  #, ax=ax1)

    n_apertures = len(tbl)
    sep = np.zeros(n_apertures)
    dra = np.zeros(n_apertures)
    ddec = np.zeros(n_apertures)
    xc = np.zeros(n_apertures)
    yc = np.zeros(n_apertures)

    for i in range(0, n_apertures):
        xc[i] = tbl[i]['xcentroid']
        yc[i] = tbl[i]['ycentroid']

        # Sky coords of center of detected object
        c_measured = SkyCoord.from_pixel(xc[i], yc[i], wcs=w, origin=0).icrs

        # Calculate separation between measured center and true center
        # sep[i] = c_measured.separation(c_true).arcsec
        sep[i] = c_measured.separation(c_true).deg
        # Calculate offsets in RA and DEC needed to bring measured center to the true center
        tmp = c_measured.spherical_offsets_to(c_true)
        # dra[i], ddec[i] = tmp[0].arcsec, tmp[1].arcsec
        dra[i], ddec[i] = tmp[0].deg, tmp[1].deg

    # dx = np.abs(xc - 93.)
    # dy = np.abs(yc - 192.)
    # good2 = (dx == np.min(dx)) & (dy == np.min(dy))
    # dra = dra[good2]
    # ddec = ddec[good2]

    dra = dra[np.abs(sep) == np.min(np.abs(sep))]
    # ddec = ddec[np.abs(sep) == np.min(np.abs(sep))]

    # xc = xc[np.abs(sep) == np.min(np.abs(sep))]
    # yc = yc[np.abs(sep) == np.min(np.abs(sep))]
    # sep = sep[np.abs(sep) == np.min(np.abs(sep))]

    # update wcs
    w_new = w.deepcopy()
    # dra = -0.00059043
    # ddec = 4.98643004e-05
    w_new.wcs.crval = w_new.wcs.crval + np.array([dra[0], ddec[0]])
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra, ddec])
    # w_new = w.deepcopy()
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra[2], ddec[2]])
    # ['CRVAL1'] = tempHead['CRVAL1'] + (offsetShifts[galName]['RA']/3600.0)
    # tempHead['CRVAL2'] = tempHead['CRVAL2'] + (offsetShifts[galName]['DEC']/3600.0)
    return w_new


# fname_cube = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/NIRSpec_F290LP_Orion_Bar_Sept_18_2022/Pointing_5_Level3_g395h-f290lp_s3d.fits'
# fname_synth_images = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/27Oct22/synthetic_nircam_images_both_conventions_.pk'
# fname_out_cube = f"/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/stitched_cubes/NIRSpec_stitched_allchannels_pointing_5_wcsfixed.fits"
# fname_out_synth = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/27Oct22/synthetic_nircam_images_both_conventions_pointing_5_wcsfixed.fits'
# fname_out_synth_unc = f'/Volumes/data1/jwst/data/ers1288/from_ftp/NIRSPEC/synthetic_nircam_images/27Oct22/synthetic_nircam_images_both_conventions_pointing_5_unc_wcsfixed.fits'
# i_gaia = 0


def realign_nirspec(fname_synth_images,
                    fname_cube,
                    fname_out_cube,
                    fname_out_synth,
                    fname_out_synth_unc,
                    i_gaia=1):
    # i_gaia = 1
    nircam_synth_images = pickle.load(open(fname_synth_images, 'rb'))
    w = get_2d_wcs_from_cube(fname_cube)
    footprint = w.calc_footprint(fits.open(fname_cube)[1].header)
    # ras, decs = np.meshgrid(ras, decs)
    ramin, ramax = np.min(footprint[:, 0]), np.max(footprint[:, 0])
    decmin, decmax = np.min(footprint[:, 1]), np.max(footprint[:, 1])
    ra0 = (ramin + ramax) / 2.
    dec0 = (decmin + decmax) / 2.
    search_radius = max([np.abs(ramax - ramin), np.abs(decmax - decmin)]) * 2

    coord_search = SkyCoord(ra=ra0,
                            dec=dec0,
                            unit=(u.degree, u.degree),
                            frame='icrs')
    width = u.Quantity(np.abs(ramax - ramin), u.deg)
    height = u.Quantity(np.abs(decmax - decmin), u.deg)
    cat_gaia = Gaia.query_object_async(coordinate=coord_search,
                                       width=width,
                                       height=height)
    cat_gaia.pprint(max_lines=12, max_width=130)

    c_true = SkyCoord(cat_gaia[i_gaia]['ra'] * u.deg,
                      cat_gaia[i_gaia]['dec'] * u.deg,
                      frame='icrs')

    data = nircam_synth_images['convention_a'][filt]
    data_unc = nircam_synth_images['unc_convention_a'][filt]

    bkg_estimator = MedianBackground()
    xd, yd = data.shape
    nboxx = int(xd / 150)
    nboxy = int(yd / 150)
    mmm_bkg = MMMBackground()
    coverage_mask = (data == 0) | ~np.isnan(data)
    bkg = Background2D(data, (10, 10),
                       filter_size=(3, 3),
                       bkg_estimator=bkg_estimator
                       )  #, coverage_mask=coverage_mask, fill_value=0.)
    # threshold = bkg.background + (3. * bkg.background_rms)
    threshold = 3. * bkg.background_rms
    #
    #
    # bkg_estimator = MedianBackground()
    # # bkg = Background2D(data, (5, 5),
    # bkg = Background2D(data, (10, 10),
    #                    filter_size=(3, 3),
    #                    bkg_estimator=bkg_estimator)
    # print(f"background = {bkg.background}")
    # threshold = 20
    # threshold = bkg.background + (2. * bkg.background_rms)
    # threshold = 150 #bkg.background + (5. * nircam_synth_images['unc_convention_a'][filt])
    # threshold = bkg.background + 20 * nircam_synth_images['unc_convention_a'][filt]
    fwhm = 2
    sigma = fwhm * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma)  #, xsize=1, ysize=1)
    kernel.normalize(mode='integral')
    segm = detect_sources(data - bkg.background,
                          threshold,
                          npixels=npixels,
                          kernel=kernel)
    segm_deblend = deblend_sources(data - bkg.background,
                                   segm,
                                   npixels=npixels,
                                   kernel=kernel,
                                   nlevels=32,
                                   contrast=0.001)

    cat = SourceCatalog(data, segm_deblend)
    tbl = cat.to_table()
    print(tbl)
    apertures = []
    for obj in tbl:
        position = (obj['xcentroid'], obj['ycentroid'])
        a = obj['semimajor_sigma'] * approx_radius_pix
        b = obj['semiminor_sigma'] * approx_radius_pix
        theta = obj['orientation']
        apertures.append(
            EllipticalAperture(position, a.value, b.value, theta=theta))

    fig, ax1 = pl.subplots(1, 1, figsize=(10, 10))
    ax1.imshow(data, origin='lower', cmap='Greys_r', interpolation='none')
    for aperture in apertures:
        aperture.plot(color='b', lw=1.5)  #, ax=ax1)

    pl.savefig('/arc/home/rchown/segm.png')

    n_apertures = len(tbl)
    sep = np.zeros(n_apertures)
    dra = np.zeros(n_apertures)
    ddec = np.zeros(n_apertures)
    xc = np.zeros(n_apertures)
    yc = np.zeros(n_apertures)
    kflux = np.zeros(n_apertures)
    for i in range(0, n_apertures):
        xc[i] = tbl[i]['xcentroid']
        yc[i] = tbl[i]['ycentroid']
        kflux[i] = tbl[i]['kron_flux']

        # Sky coords of center of detected object
        c_measured = SkyCoord.from_pixel(xc[i], yc[i], wcs=w, origin=0).icrs

        # Calculate separation between measured center and true center
        # sep[i] = c_measured.separation(c_true).arcsec
        sep[i] = c_measured.separation(c_true).deg
        # Calculate offsets in RA and DEC needed to bring measured center to the true center
        tmp = c_measured.spherical_offsets_to(c_true)
        # dra[i], ddec[i] = tmp[0].arcsec, tmp[1].arcsec
        dra[i], ddec[i] = tmp[0].deg, tmp[1].deg

    dx = np.abs(xc - 95.)
    dy = np.abs(yc - 192.)
    good2 = (dx == np.min(dx[kflux > 1e4])) & (dy == np.min(dy[kflux > 1e4]))
    good2 = 4
    dra = dra[good2]
    ddec = ddec[good2]
    print(dra, ddec)

    # dra = dra[np.abs(sep) == np.min(np.abs(sep))]
    # ddec = ddec[np.abs(sep) == np.min(np.abs(sep))]

    # xc = xc[np.abs(sep) == np.min(np.abs(sep))]
    # yc = yc[np.abs(sep) == np.min(np.abs(sep))]
    # sep = sep[np.abs(sep) == np.min(np.abs(sep))]

    # update wcs
    w_new = w.deepcopy()
    # dra = -0.00059043
    # ddec = 4.98643004e-05
    # 25.080192 - 30.456068
    # 16.462075 - 15.155275
    # dra = -(30.456068 - 25.080192) * w.wcs.cdelt[1]
    # ddec = (16.462075 - 15.155275) * w.wcs.cdelt[1]
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra[0], ddec[0]])
    w_new.wcs.crval = w_new.wcs.crval + np.array([dra, ddec])
    # w_new = w.deepcopy()
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra[2], ddec[2]])
    # ['CRVAL1'] = tempHead['CRVAL1'] + (offsetShifts[galName]['RA']/3600.0)
    # tempHead['CRVAL2'] = tempHead['CRVAL2'] + (offsetShifts[galName]['DEC']/3600.0)

    hdul_cube = fits.open(fname_cube)
    hdul_cube[1].header = w_new.to_header()
    hdul = fits.HDUList(hdus=hdul_cube)
    hdul.writeto(fname_out_cube, output_verify='fix', overwrite=True)

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
    hdr0 = w_new.to_header()

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


def realign_nirspec_fits(filt,
                         fname_synth_images_fits,
                         fname_synth_images_unc_fits,
                         fname_cube,
                         fname_out_cube,
                         fname_out_synth,
                         fname_out_synth_unc,
                         fname_out_wcs,
                         i_gaia=1,
                         search_gaia=False):
    # i_gaia = 1
    # nircam_synth_images = pickle.load(open(fname_synth_images, 'rb'))
    nircam_synth_images = fits.open(fname_synth_images_fits)
    nircam_synth_images_unc = fits.open(fname_synth_images_unc_fits)

    w = get_2d_wcs_from_cube(fname_cube)

    if search_gaia:
        print("Searching Gaia catalog")
        footprint = w.calc_footprint(fits.open(fname_cube)[1].header)
        # ras, decs = np.meshgrid(ras, decs)
        ramin, ramax = np.min(footprint[:, 0]), np.max(footprint[:, 0])
        decmin, decmax = np.min(footprint[:, 1]), np.max(footprint[:, 1])
        ra0 = (ramin + ramax) / 2.
        dec0 = (decmin + decmax) / 2.
        search_radius = max([np.abs(ramax - ramin),
                             np.abs(decmax - decmin)]) * 2

        coord_search = SkyCoord(ra=ra0,
                                dec=dec0,
                                unit=(u.degree, u.degree),
                                frame='icrs')
        width = u.Quantity(np.abs(ramax - ramin), u.deg)
        height = u.Quantity(np.abs(decmax - decmin), u.deg)
        cat_gaia = Gaia.query_object_async(coordinate=coord_search,
                                           width=width,
                                           height=height)
        cat_gaia.pprint(max_lines=12, max_width=130)

        c_true = SkyCoord(cat_gaia[i_gaia]['ra'] * u.deg,
                          cat_gaia[i_gaia]['dec'] * u.deg,
                          frame='icrs')

    else:
        print("using hard-coded position of proplyd to align WCS")
        c_true = SkyCoord(83.83447199597 * u.deg,
                          -5.41778904231 * u.deg,
                          frame='icrs')

    data = nircam_synth_images[filt].data
    data_unc = nircam_synth_images_unc[filt].data

    bkg_estimator = MedianBackground()
    xd, yd = data.shape
    nboxx = int(xd / 150)
    nboxy = int(yd / 150)
    mmm_bkg = MMMBackground()
    coverage_mask = (data == 0) | ~np.isnan(data)
    bkg = Background2D(data, (10, 10),
                       filter_size=(3, 3),
                       bkg_estimator=bkg_estimator
                       )  #, coverage_mask=coverage_mask, fill_value=0.)
    # threshold = bkg.background + (3. * bkg.background_rms)
    threshold = 3. * bkg.background_rms
    #
    #
    # bkg_estimator = MedianBackground()
    # # bkg = Background2D(data, (5, 5),
    # bkg = Background2D(data, (10, 10),
    #                    filter_size=(3, 3),
    #                    bkg_estimator=bkg_estimator)
    # print(f"background = {bkg.background}")
    # threshold = 20
    # threshold = bkg.background + (2. * bkg.background_rms)
    # threshold = 150 #bkg.background + (5. * nircam_synth_images['unc_convention_a'][filt])
    # threshold = bkg.background + 20 * nircam_synth_images['unc_convention_a'][filt]
    fwhm = 2
    sigma = fwhm * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma)  #, xsize=1, ysize=1)
    kernel.normalize(mode='integral')
    segm = detect_sources(data - bkg.background,
                          threshold,
                          npixels=npixels,
                          kernel=kernel)
    segm_deblend = deblend_sources(data - bkg.background,
                                   segm,
                                   npixels=npixels,
                                   kernel=kernel,
                                   nlevels=32,
                                   contrast=0.001)

    cat = SourceCatalog(data, segm_deblend)
    tbl = cat.to_table()
    print(tbl)
    # apertures = []
    # for obj in tbl:
    #     position = (obj['xcentroid'], obj['ycentroid'])
    #     a = obj['semimajor_sigma'] * approx_radius_pix
    #     b = obj['semiminor_sigma'] * approx_radius_pix
    #     theta = obj['orientation']
    #     if np.isfinite(obj['xcentroid') and np.isfinite(obj['ycentroid']):
    #         apertures.append(
    #             EllipticalAperture(position, a.value, b.value, theta=theta))
    #     else:
    #         apertures.append([])
    #
    # fig, ax1 = pl.subplots(1, 1, figsize=(10, 10))
    # ax1.imshow(data, origin='lower', cmap='Greys_r', interpolation='none')
    # for aperture in apertures:
    #     aperture.plot(color='b', lw=1.5)  #, ax=ax1)
    #
    # pl.savefig('/arc/home/rchown/segm.png')

    good_tbl_rows = np.isfinite(tbl['xcentroid']) & np.isfinite(
        tbl['ycentroid'])
    tbl = tbl[good_tbl_rows]

    n_apertures = len(tbl)
    sep = np.zeros(n_apertures)
    dra = np.zeros(n_apertures)
    ddec = np.zeros(n_apertures)
    xc = np.zeros(n_apertures)
    yc = np.zeros(n_apertures)
    kflux = np.zeros(n_apertures)
    for i in range(0, n_apertures):
        xc[i] = tbl[i]['xcentroid']
        yc[i] = tbl[i]['ycentroid']
        kflux[i] = tbl[i]['kron_flux']

        # Sky coords of center of detected object
        c_measured = SkyCoord.from_pixel(xc[i], yc[i], wcs=w, origin=0).icrs

        # Calculate separation between measured center and true center
        # sep[i] = c_measured.separation(c_true).arcsec
        sep[i] = c_measured.separation(c_true).deg
        # Calculate offsets in RA and DEC needed to bring measured center to the true center
        tmp = c_measured.spherical_offsets_to(c_true)
        # dra[i], ddec[i] = tmp[0].arcsec, tmp[1].arcsec
        dra[i], ddec[i] = tmp[0].deg, tmp[1].deg

    # dx = np.abs(xc - 95.)
    # dy = np.abs(yc - 192.)
    # good2 = (dx == np.min(dx[kflux > 1e4])) & (dy == np.min(dy[kflux > 1e4]))
    # good2 = 4
    # dra = dra[good2]
    # ddec = ddec[good2]

    # Pixel indices of closest peak
    i_min = np.argmin(sep)
    xmin, ymin = xc[i_min], yc[i_min]
    print(
        f'xmin, ymin, sep = {xmin:.4f}, {ymin:.4f}, {sep[i_min] * 3600.:.4f} arcsec'
    )
    dra = dra[i_min]
    ddec = ddec[i_min]
    print(dra, ddec)

    # dra = dra[np.abs(sep) == np.min(np.abs(sep))]
    # ddec = ddec[np.abs(sep) == np.min(np.abs(sep))]

    # xc = xc[np.abs(sep) == np.min(np.abs(sep))]
    # yc = yc[np.abs(sep) == np.min(np.abs(sep))]
    # sep = sep[np.abs(sep) == np.min(np.abs(sep))]

    # update wcs
    w_new = w.deepcopy()
    # dra = -0.00059043
    # ddec = 4.98643004e-05
    # 25.080192 - 30.456068
    # 16.462075 - 15.155275
    # dra = -(30.456068 - 25.080192) * w.wcs.cdelt[1]
    # ddec = (16.462075 - 15.155275) * w.wcs.cdelt[1]
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra[0], ddec[0]])
    w_new.wcs.crval = w_new.wcs.crval + np.array([dra, ddec])
    # w_new = w.deepcopy()
    # w_new.wcs.crval = w_new.wcs.crval + np.array([dra[2], ddec[2]])
    # ['CRVAL1'] = tempHead['CRVAL1'] + (offsetShifts[galName]['RA']/3600.0)
    # tempHead['CRVAL2'] = tempHead['CRVAL2'] + (offsetShifts[galName]['DEC']/3600.0)
    print(f"Saving corrected WCS to {fname_out_wcs}")
    res_wcs = dict()
    res_wcs['optimal_wcs'] = w_new
    with open(fname_out_wcs, 'wb') as fl:
        compress_pickle.dump(res_wcs,
                             fl,
                             compression=None,
                             set_default_extension=False)

    hdul_cube = fits.open(fname_cube)
    hdul_cube[1].header = w_new.to_header()
    hdul = fits.HDUList(hdus=hdul_cube)
    hdul.writeto(fname_out_cube, output_verify='fix', overwrite=True)

    # New HDUList for synthetic images
    # Every fits file needs a PrimaryHDU. We'll make a blank one
    hdu0 = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list = [hdu0]
    # Same for uncertainties
    hdu0_unc = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list_unc = [hdu0_unc]

    # filter_list = list(nircam_synth_images['convention_a'].keys())
    filter_list = [
        nircam_synth_images[i].header['EXTNAME']
        for i in range(1, len(nircam_synth_images))
    ]
    hdr0 = w_new.to_header()

    for filt_tmp in filter_list:
        # arr = nircam_synth_images['convention_a'][filt_tmp]
        # unc_arr = nircam_synth_images['unc_convention_a'][filt_tmp]
        arr = nircam_synth_images[filt_tmp].data
        unc_arr = nircam_synth_images_unc[filt_tmp].data
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


def realign_nirspec_arcade(outlier_off=True, reduction_date=''):
    data_dir = '/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/'
    if outlier_off:
        tag_outlier = '_outlier_off'
        if reduction_date != '':
            tag_outlier = '_' + reduction_date
    else:
        tag_outlier = ''
    fname_cube = f'{data_dir}NIRSpec_stitched_allchannels_all_pointings{tag_outlier}.fits.gz'
    fname_out_cube = f'{data_dir}NIRSpec_stitched_allchannels_all_pointings_wcsfixed{tag_outlier}.fits'  #.gz'

    output_path = '/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/synthetic_nircam_images/'
    pgroup_str = 'all_pointings'
    if reduction_date != '':
        pgroup_str = pgroup_str + '_' + reduction_date
    fname_synth_images = f'{output_path}synthetic_nircam_images_both_conventions_{pgroup_str}.pk'

    fname_out_synth = f'{output_path}synthetic_nircam_images_both_conventions_{pgroup_str}_wcsfixed.fits'
    fname_out_synth_unc = f'{output_path}synthetic_nircam_images_both_conventions_{pgroup_str}_unc_wcsfixed.fits'

    print("###############")
    print(f"fname_cube = {fname_cube}")
    print(f"fname_out_cube = {fname_out_cube}")
    print(f"fname_synth_images = {fname_synth_images}")
    print(f"fname_out_synth = {fname_out_synth}")
    print(f"fname_out_synth_unc = {fname_out_synth_unc}")
    print("###############")
    i_gaia = 1
    if reduction_date != '':
        i_gaia = 1
    realign_nirspec(fname_synth_images,
                    fname_cube,
                    fname_out_cube,
                    fname_out_synth,
                    fname_out_synth_unc,
                    i_gaia=i_gaia)


def realign_nircam():
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F470N-F444W_i2d.fits'
    filt = 'F470N'
    # fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_CLEAR-F212N_i2d.fits'
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_F470N-F444W_i2d_wcsfixed.fits'
    fname_nircam = '/Volumes/data1/jwst/data/ers1288/from_ftp/NIRCAM_final_mosaics/Level3_CLEAR-F212N_i2d_wcsfixed.fits'
    filt = 'F212N'

    f_nrc = fits.open(fname_nircam)
    w_nrc = wcs.WCS(f_nrc['SCI'].header)
    #
    # gaia_row = cat_gaia[0]
    # w_nrc_new = realign_nircam(fname_nircam, f_nrc['SCI'].data, w_nrc, gaia_row)
    #
    # dx = (4472.13 - 4472.1866) * w_nrc.wcs.cdelt[1]
    # dy = (1784.81 - 1784.0198) * w_nrc.wcs.cdelt[1]
    # w_nrc.wcs.crval = w_nrc.wcs.crval + np.array([dx, dy])
    w_nrc_new = w_nrc

    nrc_reproj, _ = reproject_exact((f_nrc['SCI'].data, w_nrc_new),
                                    w_new,
                                    arr.shape,
                                    parallel=True)
    nrc_unc_reproj, _ = reproject_exact((f_nrc['ERR'].data**2, w_nrc_new),
                                        w_new,
                                        arr.shape,
                                        parallel=True)
    nrc_unc_reproj = np.sqrt(nrc_unc_reproj)

    pl.figure()
    im = pl.imshow(np.log10(hdul[filt].data), vmin=0, vmax=1)
    cb = pl.colorbar(mappable=im, orientation='horizontal', ticklocation='top')

    cb.set_label(f'{filt} ' + 'NIRCam / synthetic', labelpad=10)

    pl.figure()
    im = pl.imshow(nrc_reproj / hdul[filt].data, vmin=0, vmax=120)
    cb = pl.colorbar(mappable=im, orientation='horizontal', ticklocation='top')

    cb.set_label(f'{filt} ' + 'NIRCam / synthetic', labelpad=10)

    # filt = 'F212N'
    pl.figure()
    pl.subplot(121, projection=w_new)
    pl.imshow(nrc_reproj)
    pl.colorbar(label=f'{filt} [MJy sr$^{-1}$]')
    pl.xlabel('RA (J2000)')
    pl.ylabel('Dec (J2000)')
    pl.subplot(122, projection=w_new)
    pl.imshow(hdul[filt].data)
    pl.xlabel('RA (J2000)')
    pl.ylabel('Dec (J2000)')
    pl.colorbar(label=f'Synthetic {filt} [MJy sr$^{-1}$]')

    pl.figure()
    # pl.scatter(np.log10(hdul[filt].data.flatten()), np.log10(nrc_reproj.flatten()))
    x = hdul[filt].data.flatten()
    xerr = hdul_unc[filt].data.flatten()
    y = nrc_reproj.flatten()
    yerr = nrc_unc_reproj.flatten()
    snr_x = np.abs(x / xerr)
    snr_y = np.abs(y / yerr)
    snrflag = (snr_x > 150.) & (snr_y > 5.)
    good = (x != 0) & (y != 0) & np.isfinite(x) & np.isfinite(y) & np.isfinite(
        xerr) & np.isfinite(yerr) & (x > 0) & snrflag  #& (x < 500) & snrflag
    # pivot = np.mean(x[good])
    reg = linregress(x[good], y[good])
    # pl.errorbar(hdul[filt].data.flatten()[good], nrc_reproj.flatten()[good], xerr= xerr[good], yerr=yerr[good], linestyle=None)
    pl.errorbar(np.log10(hdul[filt].data.flatten()[good]),
                np.log10(nrc_reproj.flatten()[good]),
                xerr=0.434 * xerr[good] / x[good],
                yerr=0.434 * yerr[good] / y[good],
                marker='o',
                markeredgecolor='xkcd:blue',
                markerfacecolor='none',
                linewidth=1.5,
                alpha=1,
                capsize=0,
                markersize=3,
                ecolor="xkcd:black",
                elinewidth=0.25,
                markeredgewidth=0.25,
                linestyle="none")
    # pl.xlabel(r'Synthetic F470N / $10^{8.58}$ [MJy sr$^{-1}$]')
    pl.xlabel(f'Synthetic {filt}' + r' [MJy sr$^{-1}$]')
    pl.ylabel(f'True {filt}' + r' [MJy sr$^{-1}$]')

    xs, ys = hdul[filt].data.flatten()[good], nrc_reproj.flatten()[good]
    xs_err = xerr[good]
    ys_err = yerr[good]

    # xs, ys, xs_err, ys_err = np.loadtxt('blah.txt').T

    fit_result = run_linmix(xs, ys, xs_err,
                            ys_err)  #, parallelize=True)# , nchains=2)
    t = np.linspace(0, 400, 500)
    plot_fit(t,
             fit_result['intercept'],
             fit_result['slope'],
             a_err=fit_result['intercept_err'],
             b_err=fit_result['slope_err'],
             ax=pl.gca(),
             yin=xs,
             xin=ys,
             pivot=fit_result['pivot'],
             pivot_y=fit_result['pivot_y'],
             log=False,
             color='black',
             lw=1,
             alpha=0.5,
             text=True,
             fontsize=9)
    np.savetxt('blah.txt', np.vstack([xs, ys, xs_err, ys_err]).T)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        if sys.argv[1] == 'realign_nirspec_arcade':
            realign_nirspec_arcade(outlier_off=True, reduction_date='12Nov22')
