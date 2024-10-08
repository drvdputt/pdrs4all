import stitch
import glob
import xcal
import compress_pickle as pickle
import sys

# curl -T /arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/stitched_cubes/all_segments/NIRSpec_stitched_allchannels_all_pointings_calibrated_wcscorr.fits.gz ftp://jwst-1288-ftp.ias.u-psud.fr/full_archive/Backups_all_versions/NIRSPEC/stitched_cubes/crds1084_1overfcorr_25May23/not_final/ --user $FTP_USER:$FTP_PASS

snr_thresh = 1.0

# /arc/projects/PDRs4All/NIRSpec_reduction_crds_1027/Filter_100LP/
# /arc/projects/PDRs4All/NIRSpec_reduction_crds_1027/stitching/
path_to_throughputs = "/arc/home/rchown/"
fnames_f100lp = glob.glob("/arc/projects/PDRs4All/NIRSpec_1_f_corr/seg1/*")
fnames_f170lp = glob.glob("/arc/projects/PDRs4All/NIRSpec_1_f_corr/seg2/*")
fnames_f290lp = glob.glob("/arc/projects/PDRs4All/NIRSpec_1_f_corr/seg3/*")
output_dir_wcs = "/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/"
output_dir_reprojected_s3d = (
    "/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/reprojected_s3d/"
)
output_dir_stitched_cubes = (
    "/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/stitched_cubes/all_segments/"
)
fname_synth_images_fits = "/arc/projects/PDRs4All/NIRSpec_1_f_corr/synthetic_images/synthetic_nircam_images_convention_A_all_pointings.fits.gz"
fname_synth_images_unc_fits = "/arc/projects/PDRs4All/NIRSpec_1_f_corr/synthetic_images/synthetic_nircam_images_convention_A_all_pointings_unc.fits.gz"

fnames_nircam_bkgsub = glob.glob(
    "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRCAM/OrionBar_NIRCAM_final_mocaics_Boris_no1overf_bkgsubtracted_19Dec22/*-B*"
)
# 19 May 2023
# - Commented out line above
# - Files below are NOT background subtracted despite the variable name (wanted to keep it the same)
# - I want to see how different the slopes are from bkg sub to not bkg subbed
# - Also would allow us to calibrate first segment
fnames_nircam_bkgsub = glob.glob(
    "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRCAM/Final_products_Boris_aligned_corrected_10Oct22/Level3*-B*"
)

# Downloaded these from the FTP server on 1 June 2023
fnames_nircam_bkgsub = glob.glob(
    "/arc/home/rchown/sep4/lib/jwst-1288-ftp.ias.u-psud.fr/full_archive/PDRs4All_final_products/NIRCAM/Level3*-B*"
)


def fname_nircam_bkgsub_reproj(fname_nircam):
    fname_end = fname_nircam.split("/")[-1]
    fname_end = fname_end.replace(".fits", "_reproj.fits")
    return f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/NIRCam_bkgsub_reprojected_wcscorr/{fname_end}"


def fname_synth_images_fits_seg(seg):
    return f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/synthetic_images/synthetic_nircam_images_convention_A_all_pointings_{seg}.fits.gz"


def fname_synth_images_unc_fits_seg(seg):
    return f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/synthetic_images/synthetic_nircam_images_convention_A_all_pointings_unc_{seg}.fits.gz"


# fname_synth_images_fits_seg('f290lp')
# output_dir_stitched_cube = '/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/reprojected_s3d/f100lp/'

# create a Stitcher
stest = stitch.Stitcher(
    fnames_s3d_f100lp=fnames_f100lp,
    fnames_s3d_f170lp=fnames_f170lp,
    fnames_s3d_f290lp=fnames_f290lp,
    output_dir_wcs=output_dir_wcs,
    output_dir_reprojected_s3d=output_dir_reprojected_s3d,
    output_dir_stitched_cubes=output_dir_stitched_cubes,
)

# if __name__ == '__main__':
#     import sys
#     fun = sys.argv[1]
# Get WCS
# stest.calculate_optimal_wcs_for_nirspec_mosaic()
# # w, s = stest.load_optimal_wcs()
# # print(w, s)
#
# if __name__ == '__main__':
#     import sys
#     f = int(sys.argv[1])
#     unc = bool(int(sys.argv[2]))
#     if f == 1:
#         filt = 'f100lp'
#     if f == 2:
#         filt = 'f290lp'
#     if f == 3:
#         filt = 'f170lp'
#     print(f"filter = {filt}, uncertainty = {unc}")
#     stest.reproject_set_of_cubes_with_pattern([filt], unc)
#

# #
# #
# stest.reproject_set_of_cubes_with_pattern(['f100lp'], True)
# stest.reproject_set_of_cubes_with_pattern(['f100lp'], False)
# stest.reproject_set_of_cubes_with_pattern(['f290lp'], True)
# stest.reproject_set_of_cubes_with_pattern(['f290lp'], False)
# stest.reproject_set_of_cubes_with_pattern(['f170lp'], True)
# stest.reproject_set_of_cubes_with_pattern(['f170lp'], False)
#
# for lbl in ['1', '3', '5', '7', '9', 'b', 'd', 'f', 'h']:
#     stest.stitch_and_save_single_pointing(pointing_label=f'Pointing_{lbl}')
# # #
# plist = ['1', '3', '5', '7', '9', 'b', 'd', 'f', 'h']
# plist = ['Pointing_' + p for p in plist]
# stest.coadd_pointings(plist, seg='all')
# #
# #
# stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits,
#                                               fname_synth_images_unc_fits,
#                                               path_to_throughputs=path_to_throughputs)
#
# stest.realign_nirspec_fits('F212N', fname_synth_images_fits, fname_synth_images_unc_fits)
# # # #
# stest.coadd_pointings(plist, seg='f100lp')
# stest.coadd_pointings(plist, seg='f170lp')
# stest.coadd_pointings(plist, seg='f290lp')
# # # # #
# print("F100LP")
# stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f100lp'),
#                                               fname_synth_images_unc_fits_seg('f100lp'),
#                                               path_to_throughputs=path_to_throughputs, seg='f100lp')
# print("F170LP")
#
# stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f170lp'),
#                                               fname_synth_images_unc_fits_seg('f170lp'),
#                                               path_to_throughputs=path_to_throughputs, seg='f170lp')
# print("F290LP")
# stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f290lp'),
#                                               fname_synth_images_unc_fits_seg('f290lp'),
#                                               path_to_throughputs=path_to_throughputs, seg='f290lp')
# #
# # for fname_nircam in fnames_nircam_bkgsub:
# #     stest.reproject_nircam(fname_nircam, fname_nircam_bkgsub_reproj(fname_nircam))
# snr_thresh = 1.
# # #
# print("F405N")
# xc = stest.measure_xcal('F405N',
#                      fname_synth_images_fits_seg('f290lp'),
#                      fname_synth_images_unc_fits_seg('f290lp'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[-1]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
# # import stitch
# # import glob
# # import xcal
# # import compress_pickle as pickle
#
# # xc = pickle.load(open('/arc/home/rchown/xc_f405n_crds1027.pkl', 'rb'),
# #                            compression=None)
# #
# #
# # xcal.plot_xcal(xc)
#
# with open('/arc/home/rchown/xc_f405n_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
# #
# # #
# # # #
# print("F335M")
#
# xc = stest.measure_xcal('F335M',
#                      fname_synth_images_fits_seg('f290lp'),
#                      fname_synth_images_unc_fits_seg('f290lp'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[-2]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f335m_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
# # #
# # # # '''
# # # # x = c * x
# # # # sig_x = np.abs(c * x) * np.sqrt((sig_c/c)**2 + (sig_x/x)**2)
# # # # '''
# # # #
# xc = stest.measure_xcal('F210M',
#                      fname_synth_images_fits_seg('f170lp'),
#                      fname_synth_images_unc_fits_seg('f170lp'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[1]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f210m_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
# #
# # #
# xc = stest.measure_xcal('F182M',
#                      fname_synth_images_fits_seg('f170lp'),
#                      fname_synth_images_unc_fits_seg('f170lp'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[0]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f182m_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
#
#
# fnames_xcal = ['/arc/home/rchown/xc_f182m_crds1027.pkl', '/arc/home/rchown/xc_f210m_crds1027.pkl']
# # stest.coadd_pointings(plist, seg='f100lp')
# stest.coadd_pointings(plist, seg='f170lp', fnames_xcal=fnames_xcal)
# fnames_xcal = ['/arc/home/rchown/xc_f405n_crds1027.pkl', '/arc/home/rchown/xc_f335m_crds1027.pkl']
# stest.coadd_pointings(plist, seg='f290lp', fnames_xcal=fnames_xcal)

################
# Same as above, but with calibrations applied
###############

# Pa gamma
# 1.09411 um
# Br gamma
# 2.16612 um

# with open(
#         f'/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{filt.lower()}_not_cal_not_scaled_13May23.pkl',
#         'wb') as fl:

# fnames_xcal_f100lp = []
# fnames_xcal_f170lp = [
#     '/arc/home/rchown/xc_f182m_crds1027.pkl',
#     '/arc/home/rchown/xc_f210m_crds1027.pkl'
# ]
# fnames_xcal_f290lp = [
#     '/arc/home/rchown/xc_f405n_crds1027.pkl',
#     '/arc/home/rchown/xc_f335m_crds1027.pkl'
# ]


# # version_tag = '13May23'
# # version_tag = '30May23'
# version_tag = 'corner_only_mask_9June23'
def fname_xcal_1overf(version_tag_xcal, filt):
    return f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{filt.lower()}_not_cal_not_scaled_{version_tag_xcal}.pkl"


def load_cal_dict(version_tag_xcal):
    """
    TODO (9 June 2023):
    - change this to use a different set of filters depending on how things look with non bkg subbed data
    """

    if version_tag_xcal == "corner_only_mask_9June23":
        # 10 Aug 2023:
        #   Verified (Slack msg https://ers-pdrs.slack.com/archives/C046838V9A7/p1686668506907289)
        #   that the mosaic from 13 June 2023
        #   /arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/stitched_cubes/all_segments/
        #        ...NIRSpec_stitched_allchannels_all_pointings_calibrated_wcscorr.fits.gz
        #   used the following filters to calibrate the three segments:
        #       F100LP: same cal factor as F170LP
        #       F170LP: average of 'f187n', 'f210m', 'f212n'
        #       F290LP: average of 'f323n', 'f405n', 'f470n', 'f480m'

        fnames_xcal_f100lp = []

        flist = ["f187n", "f210m", "f212n"]
        fnames_xcal_f170lp = [
            fname_xcal_1overf(version_tag_xcal, filt) for filt in flist
        ]
        # fnames_xcal_f170lp = [fname_xcal_1overf(version_tag_xcal, 'f182m'), fname_xcal_1overf(version_tag_xcal, 'f210m')]

        flist = ["f323n", "f405n", "f470n", "f480m"]

        fnames_xcal_f290lp = [
            fname_xcal_1overf(version_tag_xcal, filt) for filt in flist
        ]

        # fnames_xcal_f290lp = [
        #     fname_xcal_1overf(version_tag_xcal, 'f405n'),
        #     fname_xcal_1overf(version_tag_xcal, 'f335m')
        # ]
    else:
        fnames_xcal_f100lp = []
        fnames_xcal_f170lp = [
            fname_xcal_1overf(version_tag_xcal, "f182m"),
            fname_xcal_1overf(version_tag_xcal, "f210m"),
        ]
        fnames_xcal_f290lp = [
            fname_xcal_1overf(version_tag_xcal, "f405n"),
            fname_xcal_1overf(version_tag_xcal, "f335m"),
        ]

    cal_dict = xcal.make_cal_dict(
        fnames_xcal_f100lp, fnames_xcal_f170lp, fnames_xcal_f290lp
    )
    print(cal_dict)
    return cal_dict, fnames_xcal_f100lp, fnames_xcal_f170lp, fnames_xcal_f290lp


def plot_all_xcal_results(version_tag_xcal, plot_dir):
    import matplotlib.pyplot as pl

    # filters = ['f182m', 'f210m', 'f405n', 'f335m']
    filt_combos = [
        "CLEAR-F182M",
        "CLEAR-F140M",
        "F164N-F150W2",
        "CLEAR-F300M",
        "CLEAR-F210M",
        "CLEAR-F335M",
        "CLEAR-F277W",
        "CLEAR-F480M",
        "F162M-F150W2",
        "CLEAR-F212N",
        "F323N-F322W2",
        "F470N-F444W",
        "CLEAR-F187N",
        "F405N-F444W",
    ]
    filters = [
        fc.split("-")[-1] if "CLEAR" in fc else fc.split("-")[0] for fc in filt_combos
    ]

    for filt in sorted(filters):
        cal_dict = pickle.load(fname_xcal_1overf(version_tag_xcal, filt))
        xcal.plot_xcal(cal_dict, filt)
        pl.title(filt.upper())
        pl.savefig(plot_dir + f"xcal_{filt}.pdf")
        # try:
        #     cal_dict = pickle.load(fname_xcal_1overf(version_tag_xcal, filt))
        #     xcal.plot_xcal(cal_dict, filt)
        #     pl.title(filt.upper())
        #     pl.savefig(plot_dir + 'xcal_{filt}.pdf')
        # except:
        #     print(f"skipped {filt}")


if __name__ == "__main__":
    # import argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-c", "--corners_only", action="store_true")

    # parser.add_argument("-c", "--corners_only", action="store_true")
    # parsed_args = parser.parse_args()
    # print(parsed_args)
    # if parsed_args.corners_only:
    #     print("corners_only flag given")
    # else:
    #     print("corners_only flag not given")
    mask_top_left_bottom_right = False
    version_tag_xcal = "13May23"  # Using backgorund subtracted nircam
    # version_tag_xcal = '30May23'  # Using non-backgorund subtracted nircam

    if "--corners_only" in sys.argv:
        print("corners")
        mask_top_left_bottom_right = True
        version_tag_xcal = (
            "corner_only_mask_9June23"  # non-background subtracted nircam
        )

    fun = sys.argv[1]
    """
    # 10 May 23:
    
    - Copied this script from stitch_arcade_crds1027.py
    - Adjusted filenames for 1/f-corrected cubes, except for xcal files. Still from previous version.
    - Need to wait for seg2 reduction to finish running before i start this

    # To do:
    1. Add functions for uncalibrated data
    2. Run as many steps as possible on uncalibrated data
    3. Make new xcal dictionary, compare with previous
    4. Make all the versions of stitched cubes once xcal looks okay


    """
    if fun == "step0_mkdir":
        import os

        os.system(
            "mkdir /arc/projects/PDRs4All/NIRSpec_1_f_corr/NIRCam_bkgsub_reprojected_wcscorr"
        )
        os.system("mkdir -p /arc/projects/PDRs4All/NIRSpec_1_f_corr/synthetic_images/")
        os.system(
            "mkdir -p /arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/stitched_cubes/all_segments/"
        )
        os.system(
            "mkdir -p /arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/reprojected_s3d"
        )

    if fun == "step1_wcs":
        # Run once
        # Get WCS
        stest.calculate_optimal_wcs_for_nirspec_mosaic()

    if fun == "step2_reproj":
        # Reproject
        # e.g. arguments:
        #   step2_reproj f100lp is_unc
        # or
        #   step2_reproj f100lp not_unc
        filt = sys.argv[2].lower()
        unc = sys.argv[3]

        # Make sure arguments are valid:
        valid_filters = ["f100lp", "f170lp", "f290lp"]
        valid_unc = ["is_unc", "not_unc", "both"]
        if filt not in valid_filters:
            print(f"Filter ({filt}) must be one of {valid_filters}")
        if unc not in valid_unc:
            print(f"Uncertainty flag ({unc}) must be one of {valid_unc}")

        if unc == "is_unc":
            stest.reproject_set_of_cubes_with_pattern([filt], True)
        if unc == "not_unc":
            stest.reproject_set_of_cubes_with_pattern([filt], False)
        if unc == "both":
            stest.reproject_set_of_cubes_with_pattern([filt], True)
            stest.reproject_set_of_cubes_with_pattern([filt], False)

    # step3_coadd_full_mosaic
    # step3_realign_full_mosaic
    # step3_coadd_single_seg_mosaics
    # step3_coadd_single_seg_mosaics_cal

    # Make single-segment mosaics
    # - No scaling, not calibrated
    # - These are to be used for measuring xcal factors
    if fun == "step3_stitch_pointings":
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        for lbl in plist:
            stest.stitch_and_save_single_pointing(
                pointing_label=f"Pointing_{lbl}",
                cal_dict=None,
                snr_cutoff=0.1,
                force_zero_scale_factor=True,
                do_not_scale=True,
            )

    if fun == "step3_coadd_full_mosaic":
        # Make single-segment mosaics
        # - No scaling, not calibrated
        # - These are to be used for measuring xcal factors
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]

        plist = ["Pointing_" + p for p in plist]
        coadd_masks = stest.coadd_pointings(
            plist,
            seg="all",
            wcscorr=False,
            calibrated=False,
            scaled=False,
            mask_top_left_bottom_right=mask_top_left_bottom_right,
        )
        with open(
            "/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/test_coadd_masks_8June23/coadd_masks_topleft_bottomright.pkl",
            "wb",
        ) as fl:
            pickle.dump(coadd_masks, fl, compression=None, set_default_extension=False)

    if fun == "step3_synth_full_mosaic":
        stest.make_and_save_synthetic_images_mosaic(
            fname_synth_images_fits,
            fname_synth_images_unc_fits,
            path_to_throughputs=path_to_throughputs,
            wcscorr=False,
            calibrated=False,
            scaled=False,
        )

    if fun == "step3_realign_full_mosaic":
        stest.realign_nirspec_fits(
            "F212N",
            fname_synth_images_fits,
            fname_synth_images_unc_fits,
            calibrated=False,
            scaled=False,
        )

    if fun == "step3_coadd_single_seg_mosaics":
        # Make single-segment mosaics
        # - No scaling, not calibrated
        # - These are to be used for measuring xcal factors
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        plist = ["Pointing_" + p for p in plist]
        for seg in ["f100lp", "f170lp", "f290lp"]:
            _ = stest.coadd_pointings(
                plist,
                seg=seg,
                wcscorr=True,
                calibrated=False,
                scaled=False,
                mask_top_left_bottom_right=mask_top_left_bottom_right,
            )

    if fun == "step3_coadd_single_seg_mosaics_cal":
        # Make single-segment mosaics
        # - No scaling, not calibrated
        # - These are to be used for measuring xcal factors
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        # for seg in ['f100lp', 'f170lp', 'f290lp']:
        plist = ["Pointing_" + p for p in plist]
        cal_dict, fnames_xcal_f100lp, fnames_xcal_f170lp, fnames_xcal_f290lp = (
            load_cal_dict(version_tag_xcal)
        )
        fname_calfac_lists = [
            fnames_xcal_f100lp,
            fnames_xcal_f170lp,
            fnames_xcal_f290lp,
        ]
        # Note: need to pass fname_calfac_lists for single seg calibrated mosaics, but not
        for j, seg in enumerate(["f100lp", "f170lp", "f290lp"]):
            _ = stest.coadd_pointings(
                plist,
                seg=seg,
                fnames_xcal=fname_calfac_lists[j],
                wcscorr=True,
                calibrated=True,
                scaled=False,
                mask_top_left_bottom_right=mask_top_left_bottom_right,
            )

    if fun == "step4_prepare_nircam_images":
        for fname_nircam in fnames_nircam_bkgsub:
            stest.reproject_nircam(
                fname_nircam, fname_nircam_bkgsub_reproj(fname_nircam)
            )

    if fun == "step5_synth_single_seg":
        # Not needed? xcal is done on stitched cubes
        for seg in ["f100lp", "f170lp", "f290lp"]:
            stest.make_and_save_synthetic_images_mosaic(
                fname_synth_images_fits_seg(seg),
                fname_synth_images_unc_fits_seg(seg),
                path_to_throughputs=path_to_throughputs,
                seg=seg,
                wcscorr=True,
                scaled=False,
                calibrated=False,
            )
            # stest.realign_nirspec_fits('F212N', fname_synth_images_fits,
            #                            fname_synth_images_unc_fits)

    if fun == "step6_xcal_all_ptg_not_cal_not_scaled":
        # Cross calibration, using mosaic from
        #   all pointings, NOT calibrated, segments NOT scaled to match
        filt = sys.argv[2]
        snr_thresh = 1.0

        # if filt == 'F210M':
        #     i = 1
        # if filt == 'F182M':
        #     i = 0
        # if filt == 'F335M':
        #     i = -2
        # if filt == 'F405N':
        #     i = -1

        # stitching script not using right indices for filters when using backgorund vs not bkg subbed images
        filts_without_clear = ["F164N", "F162M", "F323N", "F470N", "F405N"]

        def search_list_of_filenames_for_pattern(list_of_fnames, pattern):
            for fname in list_of_fnames:
                if pattern in fname:
                    return fname
            return ""

        filt_combos = [
            "CLEAR-F182M",
            "CLEAR-F140M",
            "F164N-F150W2",
            "CLEAR-F300M",
            "CLEAR-F210M",
            "CLEAR-F335M",
            "CLEAR-F277W",
            "CLEAR-F480M",
            "F162M-F150W2",
            "CLEAR-F212N",
            "F323N-F322W2",
            "F470N-F444W",
            "CLEAR-F187N",
            "F405N-F444W",
        ]

        if filt == "all":
            # ilist = [1, 0, -2, -1]
            # filters = ['F210M', 'F182M', 'F335M', 'F405N']
            filters = sorted(
                [
                    fc.split("-")[-1] if "CLEAR" in fc else fc.split("-")[0]
                    for fc in filt_combos
                ]
            )[:3]
        else:
            filters = [filt]

        fnames_nircam_temp = []
        for ftemp in filters:
            if ftemp in filts_without_clear:
                look_for = f"_{ftemp}"
            else:
                look_for = f"CLEAR-{ftemp}"

            fname = search_list_of_filenames_for_pattern(fnames_nircam_bkgsub, look_for)
            fnames_nircam_temp.append(fname)

        print(f"These are the filters: {filters}")
        print(f"These are the filenames: {fnames_nircam_temp}")

        skipped_filters = []
        for fname_nircam_temp, ftemp in zip(fnames_nircam_temp, filters):
            try:
                xc = stest.measure_xcal(
                    ftemp,
                    fname_synth_images_fits,
                    fname_synth_images_unc_fits,
                    fname_nircam_bkgsub_reproj(fname_nircam_temp),
                    return_xc=True,
                    edge_dist_thresh_pix=9,
                    region_mask=None,
                    snr_thresh=snr_thresh,
                )

                # with open(
                #         f'/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{filtname.lower()}_not_cal_not_scaled_13May23.pkl',
                #         'wb') as fl:
                #     pickle.dump(xc, fl, compression=None, set_default_extension=False)
                with open(
                    f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{ftemp.lower()}_not_cal_not_scaled_{version_tag_xcal}.pkl",
                    "wb",
                ) as fl:
                    pickle.dump(xc, fl, compression=None, set_default_extension=False)
            except:
                print(f"Skipped filter {ftemp}")
                skipped_filters.append(ftemp)
        print(f"Done. Skipped filters: {skipped_filters}")
    # print("F100LP")
    # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f100lp'),
    #                                               fname_synth_images_unc_fits_seg('f100lp'),
    #                                               path_to_throughputs=path_to_throughputs, seg='f100lp')
    # print("F170LP")
    #
    # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f170lp'),
    #                                               fname_synth_images_unc_fits_seg('f170lp'),
    #                                               path_to_throughputs=path_to_throughputs, seg='f170lp')
    # print("F290LP")
    # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f290lp'),
    #                                               fname_synth_images_unc_fits_seg('f290lp'),
    #                                               path_to_throughputs=path_to_throughputs, seg='f290lp')
    # #
    # # for fname_nircam in fnames_nircam_bkgsub:
    # #     stest.reproject_nircam(fname_nircam, fname_nircam_bkgsub_reproj(fname_nircam))

    if fun == "step7_stitch_pointings_calibrated_and_scaled":
        # One pointing, calibrated, scale segments to match
        do_not_scale = False
        # lbl = sys.argv[2]
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        cal_dict, fnames_xcal_f100lp, fnames_xcal_f170lp, fnames_xcal_f290lp = (
            load_cal_dict(version_tag_xcal)
        )

        for lbl in plist:
            stest.stitch_and_save_single_pointing(
                pointing_label=f"Pointing_{lbl}",
                cal_dict=cal_dict,
                snr_cutoff=0.1,
                force_zero_scale_factor=True,
            )
            scale_factor_dict = stest.scale_factor_maps
            with open(
                f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/stitching/scale_factors/scale_factors_{lbl}.pkl",
                "wb",
            ) as fl:
                pickle.dump(
                    scale_factor_dict, fl, compression=None, set_default_extension=False
                )

    if fun == "step7_stitch_pointings_calibrated_not_scaled":
        # One pointing, calibrated, scale segments to match
        do_not_scale = True
        # lbl = sys.argv[2]
        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        cal_dict, fnames_xcal_f100lp, fnames_xcal_f170lp, fnames_xcal_f290lp = (
            load_cal_dict(version_tag_xcal)
        )

        for lbl in plist:
            stest.stitch_and_save_single_pointing(
                pointing_label=f"Pointing_{lbl}",
                cal_dict=cal_dict,
                snr_cutoff=0.1,
                force_zero_scale_factor=True,
                do_not_scale=do_not_scale,
            )

    if fun == "step8_make_mosaic_calibrated_and_scaled":
        # All pointings, calibrated, scale segments to match
        # - Use corrected WCS (doesn't recompute it)
        do_not_scale = False

        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        plist = ["Pointing_" + p for p in plist]
        stest.coadd_pointings(
            plist,
            seg="all",
            wcscorr=True,
            calibrated=True,
            scaled=(not do_not_scale),
            mask_top_left_bottom_right=mask_top_left_bottom_right,
        )

        # stest.make_and_save_synthetic_images_mosaic(
        #     fname_synth_images_fits.replace('.fits',
        #                                     '_calibrated_scaled.fits'),
        #     fname_synth_images_unc_fits.replace('.fits',
        #                                         '_calibrated_scaled.fits'),
        #     path_to_throughputs=path_to_throughputs,
        #     calibrated=True,
        #     scaled=(not do_not_scale))

        # stest.realign_nirspec_fits(
        #     'F212N',
        #     fname_synth_images_fits.replace('.fits',
        #                                     '_calibrated_scaled.fits'),
        #     fname_synth_images_unc_fits.replace('.fits',
        #                                         '_calibrated_scaled.fits'),
        #     calibrated=True,
        #     scaled=(not do_not_scale))

    if fun == "step8_make_mosaic_calibrated_not_scaled":
        # All pointings, calibrated, segments NOT scaled to match
        # - Use corrected WCS (doesn't recompute it)
        do_not_scale = True

        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        plist = ["Pointing_" + p for p in plist]
        stest.coadd_pointings(
            plist,
            seg="all",
            wcscorr=True,
            calibrated=True,
            scaled=(not do_not_scale),
            mask_top_left_bottom_right=mask_top_left_bottom_right,
        )

    if fun == "all_ptg_cal_not_scaled":
        # All pointings, calibrated, DON'T scale segments to match
        do_not_scale = True
        for lbl in ["1", "3", "5", "7", "9", "b", "d", "f", "h"]:
            stest.stitch_and_save_single_pointing(
                pointing_label=f"Pointing_{lbl}",
                cal_dict=cal_dict,
                snr_cutoff=1.1,
                force_zero_scale_factor=True,
                do_not_scale=do_not_scale,
            )

        plist = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        plist = ["Pointing_" + p for p in plist]
        stest.coadd_pointings(
            plist, seg="all", calibrated=True, scaled=(not do_not_scale)
        )
        stest.make_and_save_synthetic_images_mosaic(
            fname_synth_images_fits.replace(".fits", "_calibrated.fits"),
            fname_synth_images_unc_fits.replace(".fits", "_calibrated.fits"),
            path_to_throughputs=path_to_throughputs,
            calibrated=True,
            scaled=(not do_not_scale),
        )

        stest.realign_nirspec_fits(
            "F212N",
            fname_synth_images_fits.replace(".fits", "_calibrated.fits"),
            fname_synth_images_unc_fits.replace(".fits", "_calibrated.fits"),
            calibrated=True,
            scaled=(not do_not_scale),
        )

    if fun == "xcal_all_ptg_cal_not_scaled":
        # Cross calibration, using mosaic from
        #   all pointings, calibrated, segments NOT scaled to match
        filt = sys.argv[2]
        snr_thresh = 1.0

        if filt == "F210M":
            i = 1
        if filt == "F182M":
            i = 0
        if filt == "F335M":
            i = -2
        if filt == "F405N":
            i = -1

        xc = stest.measure_xcal(
            filt,
            fname_synth_images_fits.replace(".fits", "_calibrated.fits"),
            fname_synth_images_unc_fits.replace(".fits", "_calibrated.fits"),
            fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[i]),
            return_xc=True,
            edge_dist_thresh_pix=9,
            region_mask=None,
            snr_thresh=snr_thresh,
        )

        with open(
            f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{filt.lower()}_cal_not_scaled_crds1027_27March23.pkl",
            "wb",
        ) as fl:
            pickle.dump(xc, fl, compression=None, set_default_extension=False)

    if fun == "xcal_all_ptg_cal_scaled":
        # Cross calibration, using mosaic from
        #   all pointings, calibrated, segments scaled to match
        filt = sys.argv[2]
        snr_thresh = 1.0

        if filt == "F210M":
            i = 1
        if filt == "F182M":
            i = 0
        if filt == "F335M":
            i = -2
        if filt == "F405N":
            i = -1

        xc = stest.measure_xcal(
            filt,
            fname_synth_images_fits.replace(".fits", "_calibrated_scaled.fits"),
            fname_synth_images_unc_fits.replace(".fits", "_calibrated_scaled.fits"),
            fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[i]),
            return_xc=True,
            edge_dist_thresh_pix=9,
            region_mask=None,
            snr_thresh=snr_thresh,
        )

        with open(
            f"/arc/projects/PDRs4All/NIRSpec_1_f_corr/xcal/xc_{filt.lower()}_cal_scaled_crds1027_27March23.pkl",
            "wb",
        ) as fl:
            pickle.dump(xc, fl, compression=None, set_default_extension=False)

# import os
#
# os.system(
#     f"mv {stest.fname_stitched_cube_all_pointings_calibrated} {stest.fname_stitched_cube_all_pointings_calibrated.replace('.fits', '_no_scaling.fits')}"
# )
# #
# for lbl in ['1', '3', '5', '7', '9', 'b', 'd', 'f', 'h']:
#     # stest.stitch_and_save_single_pointing(pointing_label=f'Pointing_{lbl}', cal_dict=cal_dict, snr_cutoff=1.1)
#     stest.stitch_and_save_single_pointing(pointing_label=f'Pointing_{lbl}',
#                                           cal_dict=cal_dict,
#                                           snr_cutoff=0.1,
#                                           force_zero_scale_factor=False)
#
# plist = ['1', '3', '5', '7', '9', 'b', 'd', 'f', 'h']
# plist = ['Pointing_' + p for p in plist]
# stest.coadd_pointings(plist, seg='all', calibrated=True)
# import os
# os.system(f"mv {stest.fname_stitched_cube_all_pointings_calibrated} {stest.fname_stitched_cube_all_pointings_calibrated.replace('.fits', '_no_scaling.fits')}")

# # #
# stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                                               fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                                               path_to_throughputs=path_to_throughputs, calibrated=True)
#
# stest.realign_nirspec_fits('F212N', fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                                               fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'), calibrated=True)
# # # # #
# # stest.coadd_pointings(plist, seg='f100lp')
# # stest.coadd_pointings(plist, seg='f170lp')
# # stest.coadd_pointings(plist, seg='f290lp')
# # # # # # #
# # # print("F100LP")
# # # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f100lp'),
# # #                                               fname_synth_images_unc_fits_seg('f100lp'),
# # #                                               path_to_throughputs=path_to_throughputs, seg='f100lp')
# # # print("F170LP")
# # #
# # # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f170lp'),
# # #                                               fname_synth_images_unc_fits_seg('f170lp'),
# # #                                               path_to_throughputs=path_to_throughputs, seg='f170lp')
# # # print("F290LP")
# # # stest.make_and_save_synthetic_images_mosaic(fname_synth_images_fits_seg('f290lp'),
# # #                                               fname_synth_images_unc_fits_seg('f290lp'),
# # #                                               path_to_throughputs=path_to_throughputs, seg='f290lp')
# # # #
# for fname_nircam in fnames_nircam_bkgsub:
#     stest.reproject_nircam(fname_nircam, fname_nircam_bkgsub_reproj(fname_nircam))
# #
# print("F405N")
# xc = stest.measure_xcal('F405N',
#                      fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                      fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[-1]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
# # import stitch
# # import glob
# # import xcal
# # import compress_pickle as pickle
#
# # xc = pickle.load(open('/arc/home/rchown/xc_f405n_crds1027.pkl', 'rb'),
# #                            compression=None)
# #
# #
# # xcal.plot_xcal(xc)
#
# with open('/arc/home/rchown/xc_f405n_cal_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
#
# #
# # #
# print("F335M")
#
# xc = stest.measure_xcal('F335M',
#                      fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                      fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[-2]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f335m_cal_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
# #
# # # '''
# # # x = c * x
# # # sig_x = np.abs(c * x) * np.sqrt((sig_c/c)**2 + (sig_x/x)**2)
# # # '''
# # #
# xc = stest.measure_xcal('F210M',
#                      fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                      fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[1]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f210m_cal_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
# #
# # #
# xc = stest.measure_xcal('F182M',
#                      fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                      fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                      fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[0]),
#                      return_xc=True,
#                      edge_dist_thresh_pix=9,
#                      region_mask=None, snr_thresh=snr_thresh)
#
# with open('/arc/home/rchown/xc_f182m_cal_crds1027.pkl', 'wb') as fl:
#     pickle.dump(xc,
#                 fl,
#                 compression=None,
#                 set_default_extension=False)
#
# import sys
# if __name__ == '__main__':
#     filt = sys.argv[1]
#
#     if filt == 'F210M':
#         i = 1
#     if filt == 'F182M':
#         i = 0
#     if filt == 'F335M':
#         i = -2
#     if filt == 'F405N':
#         i = -1
#
#     xc = stest.measure_xcal(filt,
#                          fname_synth_images_fits.replace('.fits', '_calibrated.fits'),
#                          fname_synth_images_unc_fits.replace('.fits', '_calibrated.fits'),
#                          fname_nircam_bkgsub_reproj(fnames_nircam_bkgsub[i]),
#                          return_xc=True,
#                          edge_dist_thresh_pix=9,
#                          region_mask=None, snr_thresh=snr_thresh)
#
#     with open(f'/arc/home/rchown/xc_{filt.lower()}_cal_crds1027_25March23.pkl', 'wb') as fl:
#         pickle.dump(xc,
#                     fl,
#                     compression=None,
#                     set_default_extension=False)
