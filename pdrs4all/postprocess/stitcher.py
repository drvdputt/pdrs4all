import glob
import numpy as np
import stitch_utils
import synth
from astropy.io import fits
import astropy.units as u
import astropy.wcs as wcs
import matplotlib.pyplot as pl
import compress_pickle as pickle


def fname_fits_add_wcscorr(fname):
    if "wcscorr" in fname:
        raise ValueError(f"Found 'wcscorr' in filename {fname}. Aborting.")
    return fname.replace(".fits", "_wcscorr.fits")


def reshape_cube(cube):
    nx, ny, nw = cube.shape
    cube_out = np.zeros((nw, nx, ny))
    for i in range(nw):
        cube_out[i] = cube[:, :, i]
    return cube_out


class Stitcher:
    """Summary of class here.

    Longer class information...
    Longer class information...

    Attributes:
        likes_spam: A boolean indicating if we like SPAM or not.
        eggs: An integer count of the eggs we have laid.
    """

    def __init__(
        self,
        fnames_s3d_f100lp,
        fnames_s3d_f170lp,
        fnames_s3d_f290lp,
        output_dir_wcs,
        output_dir_reprojected_s3d,
        output_dir_stitched_cubes,
        ifu_name="NIRSpec",
        xcal_dict=dict(),
    ):
        """Inits Stitcher"""
        # segments
        # path to s3d files
        # filenames of s3d files (take as input pointing and segment, return filename)
        # path to reprojected pickle files (same as above)
        # path to output cubes
        self.fnames_s3d_f100lp = fnames_s3d_f100lp
        self.fnames_s3d_f170lp = fnames_s3d_f170lp
        self.fnames_s3d_f290lp = fnames_s3d_f290lp

        self.ifu_name = ifu_name

        self.output_dir_wcs = output_dir_wcs
        self.fname_wcs = f"{self.output_dir_wcs}optimal_wcs_nirspec_mosaic.pkl"
        self.fname_wcs_corrected = (
            f"{self.output_dir_wcs}optimal_wcs_nirspec_mosaic_corrected.pkl"
        )

        self.output_dir_reprojected_s3d = output_dir_reprojected_s3d
        self.output_dir_stitched_cubes = output_dir_stitched_cubes
        self.fname_stitched_cube_all_pointings = f"{self.output_dir_stitched_cubes}{ifu_name}_stitched_allchannels_all_pointings.fits.gz"
        self.fname_stitched_cube_all_pointings_calibrated = f"{self.output_dir_stitched_cubes}{ifu_name}_stitched_allchannels_all_pointings_calibrated.fits.gz"
        self.fname_stitched_cube_all_pointings_wcs_corrected = fname_fits_add_wcscorr(
            self.fname_stitched_cube_all_pointings
        )
        self.fname_stitched_cube_all_pointings_wcs_corrected_calibrated = (
            fname_fits_add_wcscorr(self.fname_stitched_cube_all_pointings_calibrated)
        )
        self.fname_mask_stitched_cube_all_pointings = (
            lambda edge_dist_thresh_pix: f"{self.output_dir_stitched_cubes}minimal_edge_mask_{edge_dist_thresh_pix}pix.fits.gz"
        )

        self.xcal_dict = xcal_dict

        self.scale_factor_maps = {"seg1": dict(), "seg2": dict()}

        # If making synthetic images
        # self.path_to_throughputs = path_to_throughputs

    def load_optimal_wcs(self, wcscorr=False):
        if not wcscorr:
            with open(self.fname_wcs, "rb") as fl:
                optimal_wcs = pickle.load(fl, compression=None)
                w = optimal_wcs["optimal_wcs"]
                shape_out = optimal_wcs["shape"]
        else:
            # Look for wcs-corrected data
            with open(self.fname_wcs_corrected, "rb") as fl:
                optimal_wcs = pickle.load(fl, compression=None)
                w = optimal_wcs["optimal_wcs"]

            with open(self.fname_wcs, "rb") as fl:
                optimal_wcs = pickle.load(fl, compression=None)
                shape_out = optimal_wcs["shape"]
        return w, shape_out

    # Some filename functions
    def fname_reprojected_s3d(self, fname_cube, uncertainty):
        unc_str = ""
        if uncertainty:
            unc_str = "_unc"
        fname = (
            self.output_dir_reprojected_s3d
            + fname_cube.split("/")[-1].split(".")[0]
            + unc_str
            + ".pk.gz"
        )
        return fname

    def fnames_s3d_with_pattern(self, patterns):
        """
        Gets filenames of s3d files that match a pattern or set of patterns.
        patterns is a list of strings
        """
        fnames_s3d_all_pointings = (
            self.fnames_s3d_f100lp + self.fnames_s3d_f170lp + self.fnames_s3d_f290lp
        )
        res = [
            fname
            for fname in fnames_s3d_all_pointings
            if (sum([pat in fname for pat in patterns]) == len(patterns))
        ]
        # print("Filenames with matches:")
        # print(res)
        return res

    def fname_stitched_cube_single_pointing(
        self, pointing_label, calibrated=False, scaled=True
    ):
        cal_str = ""
        scale_str = ""
        if calibrated:
            cal_str = "_calibrated"
        if scaled:
            scale_str = "_scaled"

        return f"{self.output_dir_stitched_cubes}{self.ifu_name}_stitched_allchannels_pointing_{pointing_label}{cal_str}{scale_str}.fits.gz"

    def fname_mosaic_single_segment(
        self, seg, wcscorr=True, calibrated=False, scaled=True
    ):
        wcs_str = ""
        cal_str = ""
        scale_str = ""
        if wcscorr:
            wcs_str = "_wcscorr"
        if calibrated:
            cal_str = "_calibrated"
        if scaled:
            scale_str = "_scaled"
        return f"{self.output_dir_stitched_cubes}{self.ifu_name}_mosaic_{seg.lower()}{wcs_str}{cal_str}{scale_str}.fits.gz"

    # End filename functions

    def calculate_optimal_wcs_for_nirspec_mosaic(self, save=True):
        """
        Get optimal WCS and shape of mosaic for NIRSpec data

        TODO:
        - test that this works (done)
        """
        print("Getting optimal WCS and shape of mosaic...")
        fnames_s3d_all_pointings = (
            self.fnames_s3d_f100lp + self.fnames_s3d_f170lp + self.fnames_s3d_f290lp
        )
        print(f"Input s3d files: {fnames_s3d_all_pointings}")
        ra, dec, w = stitch_utils.get_ra_dec_grid_optimal_wcs(fnames_s3d_all_pointings)
        ras, decs = np.meshgrid(ra, dec)
        shape_out = ras.shape
        if save:
            print(f"Saving to {self.fname_wcs}")
            res = dict()
            res["optimal_wcs"] = w
            res["shape"] = shape_out
            with open(self.fname_wcs, "wb") as fl:
                pickle.dump(res, fl, compression=None, set_default_extension=False)
        else:
            return w, shape_out

    # reproject
    # one wavelength slice
    def reproject_cube(self, fname_cube, uncertainty, clip_spikes=True):
        """
        TODO:
            - define variables:
                fname_cube
                output_dir
                fname_optimal_wcs
            - add a check for optimal wcs
        """
        print("reproject_nirspec_cube:")
        print(f"   fname_cube = {fname_cube}")
        # print(f"   output_dir = {output_dir}")
        # path_to_stage3 = f"/Volumes/data1/jwst/data/ers1288/nirspec/pointing_group_{pointing_group}/stage3/"
        # path_to_data = '/Volumes/data1/jwst/data/ers1288/from_ftp/'
        # if arcade:
        #     path_to_data = '/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/'

        w, shape_out = self.load_optimal_wcs()

        cube = stitch_utils.projection_field(
            fname_cube, w, shape_out, uncertainty=uncertainty, clip_spikes=clip_spikes
        )
        # Save
        res = dict()
        res["cube"] = cube

        fname_out = self.fname_reprojected_s3d(fname_cube, uncertainty)

        print(f"Saving to {fname_out}")
        with open(fname_out, "wb") as fl:
            pickle.dump(res, fl, compression="gzip", set_default_extension=False)

    def reproject_set_of_cubes(self, list_of_fnames, uncertainty, clip_spikes=True):
        """
        e.g. self.reproject_set_of_cubes(self.fnames_s3d_f100lp, False)
        """
        for fname_cube in list_of_fnames:
            self.reproject_cube(fname_cube, uncertainty, clip_spikes=clip_spikes)

    def reproject_set_of_cubes_with_pattern(
        self, patterns, uncertainty, clip_spikes=True
    ):
        """
        e.g. self.reproject_set_of_cubes_with_pattern(['Pointing_1', 'F100LP'], False)
        """
        list_of_fnames = self.fnames_s3d_with_pattern(patterns)
        for fname_cube in list_of_fnames:
            self.reproject_cube(fname_cube, uncertainty, clip_spikes=clip_spikes)

    def reproject_all_cubes(self, uncertainty, clip_spikes=True):
        fnames_s3d_all_pointings = (
            self.fnames_s3d_f100lp + self.fnames_s3d_f170lp + self.fnames_s3d_f290lp
        )
        for fname_cube in fnames_s3d_all_pointings:
            self.reproject_cube(fname_cube, uncertainty, clip_spikes=clip_spikes)

    # make a stitched cube
    def run_ers1288_nirspec(
        self,
        pointing_label,
        cal_dict=None,
        snr_cutoff=None,
        force_zero_scale_factor=False,
        do_not_scale=False,
    ):
        """
        Args:
            (optional) cal_dict: dictionary containing calibration factors
                and their uncertainties. Should have the form
                    {'f100lp': {'calfac': 1.0, 'calfac_err': 1.0},
                    'f170lp': {...},
                    'f290lp': {...}}

        TODO:
            3Feb23:
            - Finish implementing calibration factors
            -
        """

        w, shape_out = self.load_optimal_wcs()

        list_of_fnames = (
            self.fnames_s3d_with_pattern([pointing_label, "100lp"])
            + self.fnames_s3d_with_pattern([pointing_label, "170lp"])
            + self.fnames_s3d_with_pattern([pointing_label, "290lp"])
        )

        print(list_of_fnames)
        # Make wavelength array
        waves_3d = [stitch_utils.get_wave(x) for x in list_of_fnames]
        # waves_3d = np.array([stitch_utils.get_wave(x) for x in list_of_fnames])
        waves3d = np.empty(0)

        if cal_dict is not None:
            print("Found calibration dict")
        else:
            cal_dict = {
                "f100lp": {"calfac": 1.0, "calfac_err": 0.0},
                "f170lp": {"calfac": 1.0, "calfac_err": 0.0},
                "f290lp": {"calfac": 1.0, "calfac_err": 0.0},
            }

        cal_keys = list(cal_dict.keys())

        # Load the first cube
        fname_cube = list_of_fnames[0]
        fname_cube_reproj = self.fname_reprojected_s3d(fname_cube, False)
        cube = pickle.load(open(fname_cube_reproj, "rb"), compression="gzip")["cube"]

        # Uncertainty
        fname_cube_unc_reproj = self.fname_reprojected_s3d(fname_cube, True)
        cube_unc = pickle.load(open(fname_cube_unc_reproj, "rb"), compression="gzip")[
            "cube"
        ]

        # Apply calibration factor
        def apply_calfac(cube, cube_unc, calfac, calfac_err):
            cube_unc = np.abs(cube * calfac) * np.sqrt(
                (cube_unc / cube) ** 2 + (calfac_err * calfac) ** 2
            )
            cube = cube * calfac
            return cube, cube_unc

        i_cal = 0
        if cal_dict[cal_keys[i_cal]]["calfac"] == 1.0:  # and (do_not_scale == False):
            i_cal = 1
            print(
                f"First cal factor is 1.0. Applying cal factor for 2nd seg instead ({cal_dict[cal_keys[1]]['calfac']})"
            )

        cube, cube_unc = apply_calfac(
            cube,
            cube_unc,
            cal_dict[cal_keys[i_cal]]["calfac"],
            cal_dict[cal_keys[i_cal]]["calfac_err"],
        )

        # Load the rest of the cubes and their uncertainties
        for i in range(1, len(list_of_fnames)):
            fname_cube = list_of_fnames[i]
            fname_cube_reproj = self.fname_reprojected_s3d(fname_cube, False)
            cube_i = pickle.load(open(fname_cube_reproj, "rb"), compression="gzip")[
                "cube"
            ]

            fname_cube_unc_reproj = self.fname_reprojected_s3d(fname_cube, True)
            cube_unc_i = pickle.load(
                open(fname_cube_unc_reproj, "rb"), compression="gzip"
            )["cube"]

            cube_i, cube_unc_i = apply_calfac(
                cube_i,
                cube_unc_i,
                cal_dict[cal_keys[i]]["calfac"],
                cal_dict[cal_keys[i]]["calfac_err"],
            )

            # Accumulate cubes and uncertainties
            cube = np.concatenate([cube, cube_i], axis=-1)
            cube_unc = np.concatenate([cube_unc, cube_unc_i], axis=-1)

        for i in range(len(list_of_fnames)):
            waves3d = np.concatenate([waves3d, waves_3d[i]])

        # method = "cross cut"
        # if concat:
        #     method = "concat"
        # else:
        #     method = "cross cut"
        method = "cross cut"

        nx, ny = cube.shape[0], cube.shape[1]
        spec_cube = []
        # waves_all = []
        # flux_all = []
        # flux_all_err = []

        ind_1 = 0
        ind_2 = waves_3d[0].size
        ind_3 = ind_2 + waves_3d[1].size
        ind_4 = ind_3 + waves_3d[2].size
        indices = [ind_1, ind_2, ind_3, ind_4]

        # Data for segment 1
        j = 0
        wave1 = waves_3d[j]
        cube1 = cube[:, :, : indices[j + 1]]
        cube_unc1 = cube_unc[:, :, : indices[j + 1]]

        # Data for segment 2
        j = 1
        wave2 = waves_3d[j]
        cube2 = cube[:, :, indices[j] : indices[j + 1]]
        cube_unc2 = cube_unc[:, :, indices[j] : indices[j + 1]]

        # Data for segment 3
        j = 2
        wave3 = waves_3d[j]
        cube3 = cube[:, :, indices[j] : indices[j + 1]]
        cube_unc3 = cube_unc[:, :, indices[j] : indices[j + 1]]

        if do_not_scale:
            scale_factor_map_seg1 = np.ones((nx, ny))
            scale_factor_map_seg2 = np.ones((nx, ny))
            where_to_calc_scale_factors = np.zeros((nx, ny))
        else:
            scale_factor_map_seg1 = np.zeros((nx, ny))
            scale_factor_map_seg2 = np.zeros((nx, ny))
            where_to_calc_scale_factors = np.zeros((nx, ny))
            # 1 = seg1
            # 2 = seg2
            # 3 = both
            # Loop for scaling factors
            cube1_mask = np.sum((cube1 != 0) & (~np.isnan(cube1)), axis=2)
            cube2_mask = np.sum((cube2 != 0) & (~np.isnan(cube2)), axis=2)
            cube3_mask = np.sum((cube3 != 0) & (~np.isnan(cube3)), axis=2)

            pix_radius = 0
            for x in range(0, nx):
                for y in range(ny):
                    result = 0
                    if cube1_mask[x, y] + cube2_mask[x, y] + cube3_mask[x, y] > 0:
                        print("x, y = ", x, y)
                        res_seg1 = stitch_utils.should_we_calc_scale_factor(
                            x,
                            y,
                            wave1,
                            cube1,
                            cube_unc1,
                            wave2,
                            cube2,
                            cube_unc2,
                            pix_radius=pix_radius,
                            n_min_overlap=10,
                        )

                        res_seg2 = stitch_utils.should_we_calc_scale_factor(
                            x,
                            y,
                            wave2,
                            cube2,
                            cube_unc2,
                            wave3,
                            cube3,
                            cube_unc3,
                            pix_radius=pix_radius,
                            n_min_overlap=10,
                        )

                        if res_seg1 and (not res_seg2):
                            result = 1
                        if (not res_seg1) and res_seg2:
                            result = 2
                        if res_seg1 and res_seg2:
                            result = 3
                    where_to_calc_scale_factors[x, y] = result

            # Loop for scaling factors
            for x in range(nx):
                for y in range(ny):
                    # Note:
                    # where_to_calc_scale_factors = 1 --> calculate scaling factor for seg 1
                    # where_to_calc_scale_factors = 2 --> calculate scaling factor for seg 2
                    # where_to_calc_scale_factors = 3 --> calculate scaling factor for seg 1 and 2
                    # where_to_calc_scale_factors = 0 --> don't calculate

                    print("x, y = ", x, y)
                    scale_factor_seg2_seg3 = 1.0
                    scale_factor_seg1_seg2 = 1.0
                    if cube1_mask[x, y] + cube2_mask[x, y] + cube3_mask[x, y] > 0:
                        # for k, wv in enumerate([wave1, wave2, wave3]):
                        #     print(
                        #         f"Segment {k}: wavelengths {wv.min()} um to {wv.max()} um"
                        #     )

                        if where_to_calc_scale_factors[x, y] == 3:
                            print("Calculating scaling factor for seg 1 and 2")
                            # Scale factor for segments 2 and 3
                            scale_factor_seg2_seg3 = (
                                stitch_utils.measure_scale_factor_ix_iy(
                                    x,
                                    y,
                                    wave2,
                                    cube2,
                                    cube_unc2,
                                    wave3,
                                    cube3,
                                    cube_unc3,
                                    pix_radius=pix_radius,
                                )
                            )

                            scale_factor_seg1_seg2 = (
                                stitch_utils.measure_scale_factor_ix_iy(
                                    x,
                                    y,
                                    wave1,
                                    cube1,
                                    cube_unc1,
                                    wave2,
                                    cube2 * scale_factor_seg2_seg3,
                                    cube_unc2 * scale_factor_seg2_seg3,
                                    pix_radius=pix_radius,
                                )
                            )

                        if where_to_calc_scale_factors[x, y] == 2:
                            print("Calculating scaling factor for seg 2")
                            # Scale factor for segments 2 and 3
                            scale_factor_seg2_seg3 = (
                                stitch_utils.measure_scale_factor_ix_iy(
                                    x,
                                    y,
                                    wave2,
                                    cube2,
                                    cube_unc2,
                                    wave3,
                                    cube3,
                                    cube_unc3,
                                    pix_radius=pix_radius,
                                )
                            )

                        if where_to_calc_scale_factors[x, y] == 1:
                            print("Calculating scaling factor for seg 1")
                            # Scale factor for segments 1 and 2
                            # Apply scale factor to segment 2 first
                            scale_factor_seg1_seg2 = (
                                stitch_utils.measure_scale_factor_ix_iy(
                                    x,
                                    y,
                                    wave1,
                                    cube1,
                                    cube_unc1,
                                    wave2,
                                    cube2 * scale_factor_seg2_seg3,
                                    cube_unc2 * scale_factor_seg2_seg3,
                                    pix_radius=pix_radius,
                                )
                            )

                        if scale_factor_seg1_seg2 == 1.0:
                            print(
                                f"Scale factor for seg 1 is 1.0. Applying scale factor for 2nd segment instead ({scale_factor_seg2_seg3})"
                            )
                            scale_factor_seg1_seg2 = scale_factor_seg2_seg3

                        # # Apply scale factor to segment 1
                        # cube1 = cube1 * scale_factor_seg1_seg2
                        # cube_unc1 = cube_unc1 * scale_factor_seg1_seg2

                    # Fill in entries of scale factor maps to use in the following loop
                    scale_factor_map_seg1[x, y] = scale_factor_seg1_seg2
                    scale_factor_map_seg2[x, y] = scale_factor_seg2_seg3

        self.scale_factor_maps["seg1"][pointing_label] = {
            "scale_factor": scale_factor_map_seg1,
            "mask": where_to_calc_scale_factors,
        }
        self.scale_factor_maps["seg2"][pointing_label] = {
            "scale_factor": scale_factor_map_seg2,
            "mask": where_to_calc_scale_factors,
        }

        # Loop for stitching
        for x in range(nx):
            for y in range(ny):
                print("x, y = ", x, y)

                # # for k, wv in enumerate([wave1, wave2, wave3]):
                # #     print(
                # #         f"Segment {k}: wavelengths {wv.min()} um to {wv.max()} um"
                # #     )
                #
                # # Scale factor for segments 2 and 3
                # scale_factor_seg2_seg3 = stitch_utils.measure_scale_factor_ix_iy(
                #     x, y, wave2, cube2, cube_unc2, wave3, cube3, cube_unc3)
                #
                # # Scale factor for segments 1 and 2
                #
                # # Apply scale factor to segment 2 first
                # cube2 = cube2 * scale_factor_seg2_seg3
                # cube_unc2 = cube_unc2 * scale_factor_seg2_seg3
                #
                # scale_factor_seg1_seg2 = stitch_utils.measure_scale_factor_ix_iy(
                #     x, y, wave1, cube1, cube_unc1, wave2, cube2, cube_unc2)
                #
                # # Apply scale factor to segment 1
                # cube1 = cube1 * scale_factor_seg1_seg2
                # cube_unc1 = cube_unc1 * scale_factor_seg1_seg2

                # Stitch segments 1 and 2
                out3d = stitch_utils.stitch(
                    wave1,
                    cube1[x, y] * scale_factor_map_seg1[x, y],
                    cube_unc1[x, y] * scale_factor_map_seg1[x, y],
                    wave2,
                    cube2[x, y] * scale_factor_map_seg2[x, y],
                    cube_unc2[x, y] * scale_factor_map_seg2[x, y],
                    method=method,
                    snr_cutoff=snr_cutoff,
                    force_zero_scale_factor=force_zero_scale_factor,
                )

                # Stitch segments 2 and 3
                i = 2
                out3d = stitch_utils.stitch(
                    wave1=out3d[0],
                    flux1=out3d[1],
                    eflux1=out3d[2],
                    wave2=wave3,
                    flux2=cube3[x, y],
                    eflux2=cube_unc3[x, y],
                    method=method,
                    snr_cutoff=snr_cutoff,
                    force_zero_scale_factor=force_zero_scale_factor,
                )

                spec_cube.append(out3d)
                # waves_all.append(waves_i)
                # flux_all.append(flux_i)
                # flux_all_err.append(flux_i_err)

        # print(f"len(spec_cube) = {len(spec_cube)}")
        # print(f"spec_cube[1].shape = {spec_cube[1].shape}")
        spec_cube = np.array(spec_cube)  # , dtype=object)
        spec_cube = np.reshape(spec_cube, (nx, ny, 3, spec_cube.shape[-1]))

        # # Plot one spaxel of combined cube
        # pl.figure(figsize=(6, 3))
        # ix = int(spec_cube.shape[0] / 2)
        # iy = int(spec_cube.shape[1] / 2)
        # pl.errorbar(spec_cube[ix, iy, 0, :],
        #             spec_cube[ix, iy, 1, :],
        #             yerr=spec_cube[ix, iy, 2, :])
        # pl.xlabel(r"$\lambda$ [$\mu$m]", fontsize=14)
        # pl.ylabel(r"$I_\lambda$ [MJy sr$^{-1}$]", fontsize=14)
        # pl.tight_layout()
        # waves = waves_i
        # print(waves)
        # # for l in [waves_all, flux_all, flux_all_err]:
        # print(f"len(waves_all) = {len(waves_all)}")
        # print(f"waves_all[1].shape = {waves_all[1].shape}")
        #
        # print(f"len(flux_all) = {len(flux_all)}")
        # print(f"flux_all[0] = {flux_all[0]}")
        # print(f"flux_all[1].shape = {flux_all[1].shape}")
        #
        # print(f"len(flux_all_err) = {len(flux_all_err)}")
        # print(f"flux_all_err[0] = {flux_all_err[0]}")
        # print(f"flux_all_err[1].shape = {flux_all_err[1].shape}")
        #
        # cube = np.concatenate(flux_all)
        # print(nx*ny*waves.size)
        # cube = cube.reshape((nx, ny, waves.size))
        # cube_unc = np.concatenate(flux_all_err)
        # cube_unc = cube.reshape((nx, ny, waves.size))
        return spec_cube, w
        # return waves, cube, cube_unc, w

    def stitch_and_save_single_pointing(
        self,
        pointing_label,
        cal_dict=None,
        snr_cutoff=None,
        force_zero_scale_factor=False,
        do_not_scale=False,
    ):
        """ """
        spec_cube, w_in = self.run_ers1288_nirspec(
            pointing_label=pointing_label,
            cal_dict=cal_dict,
            snr_cutoff=snr_cutoff,
            force_zero_scale_factor=force_zero_scale_factor,
            do_not_scale=do_not_scale,
        )
        # waves, comb_cube_temp, comb_cube_var_temp, w_in = self.run_ers1288_nirspec(
        #     pointing_label=pointing_label,
        #     cal_dict=cal_dict,
        #     snr_cutoff=snr_cutoff)

        if cal_dict is not None:
            calibrated = True
        else:
            calibrated = False

        fname_out = self.fname_stitched_cube_single_pointing(
            pointing_label, calibrated=calibrated, scaled=(not do_not_scale)
        )

        waves = spec_cube[0, 0, 0, :]
        comb_cube_temp = spec_cube[:, :, 1, :]
        comb_cube_var_temp = spec_cube[:, :, 2, :]
        nx, ny, nw = comb_cube_temp.shape
        comb_cube = np.zeros((nw, nx, ny))
        comb_cube_unc = np.zeros((nw, nx, ny))
        for i in range(nw):
            comb_cube[i] = comb_cube_temp[:, :, i]
            comb_cube_unc[i] = comb_cube_var_temp[:, :, i]

        hdr0 = w_in.to_header()
        # Every fits file needs a PrimaryHDU. We'll make a blank one
        hdu0 = fits.PrimaryHDU()
        # Start an HDUList
        hdu_list = [hdu0]

        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("CUBE", "Stitched NIRSpec datacube")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("ERR", "NIRSpec datacube uncertainty")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube_unc, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One for wavelength array
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("WAVE", "Wavelength array")
        q = 1 * u.um
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=waves, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # Make an HDUList from all HDUs we want to combine into a single FITS file
        hdul = fits.HDUList(hdus=hdu_list)

        print(f"Saving to {fname_out}")
        hdul.writeto(fname_out, overwrite=True)

    def coadd_pointings(
        self,
        pointing_labels,
        seg="all",
        edge_dist_thresh_pix=6,
        fname_out_mask="",
        fnames_xcal=[],
        wcscorr=True,
        calibrated=False,
        scaled=True,
        mask_top_left_bottom_right=False,
    ):
        """
        Coadd all single-pointing stitched cubes

        Args:
            seg (str):
                'all' -> 3-segment cubes
                'f100lp'
                'f170lp'
                'f290lp'
            # wcscorr (bool): if True, will look for wcs-corrected WCS and save
            #     output with that. Otherwise, will use uncorrected WCS.
        """
        from scipy import ndimage as ndi

        make_save_minimal_edge_mask = False

        calfac = 1.0
        calfac_err = 0.0
        # calibrated = False

        if len(fnames_xcal) > 0:
            calfac = 0.0
            calfac_err = 0.0
            sum_weights = 0.0
            for fname in fnames_xcal:
                with open(fname, "rb") as fl:
                    res_xcal = pickle.load(fl, compression=None)
                slope = res_xcal["slope"]
                slope_err = res_xcal["slope_err"]
                print(f"fname={fname}")
                print(f"slope= {slope}")
                print(f"slope_err= {slope_err}")
                wt = 1.0 / slope_err**2
                calfac += slope * wt
                sum_weights += wt
                calfac_err += (slope_err * wt) ** 2

            calfac = calfac / sum_weights
            calfac_err = np.sqrt(calfac_err / (sum_weights**2))
            print(f"calfac= {calfac}")
            print(f"calfac_err= {calfac_err}")

            calibrated = True

        if seg == "all":
            print("All segments")
            if wcscorr:
                fname_out = self.fname_stitched_cube_all_pointings_wcs_corrected
                if calibrated:
                    fname_out = (
                        self.fname_stitched_cube_all_pointings_wcs_corrected_calibrated
                    )
            else:
                fname_out = self.fname_stitched_cube_all_pointings
                if calibrated:
                    fname_out = self.fname_stitched_cube_all_pointings_calibrated

            if scaled:
                fname_out = fname_out.replace(".fits", "_scaled.fits")
            make_save_minimal_edge_mask = True
            fname_out_mask = self.fname_mask_stitched_cube_all_pointings(
                edge_dist_thresh_pix
            )

            pointing_label = pointing_labels[0]
            fname_cube = self.fname_stitched_cube_single_pointing(
                pointing_label, calibrated=calibrated, scaled=scaled
            )

            with fits.open(fname_cube) as f:
                # w_in = wcs.WCS(f['CUBE'].header)
                w_in, _ = self.load_optimal_wcs(wcscorr=wcscorr)
                comb_cube = np.zeros(f["CUBE"].data.shape)
                sum_weights = np.zeros(f["CUBE"].data.shape)
                comb_cube_unc = np.zeros(f["CUBE"].data.shape)
                wave = f["WAVE"].data
        else:
            print(f"Single segment: {seg}")
            # Use one reprojected single segment cube to initialize the mosaic
            list_of_fnames = self.fnames_s3d_with_pattern(
                [pointing_labels[0], seg.replace("f", "")]
            )
            print(list_of_fnames)
            # Load the first cube
            fname_cube = list_of_fnames[0]
            fname_cube_reproj = self.fname_reprojected_s3d(fname_cube, False)
            cube = pickle.load(open(fname_cube_reproj, "rb"), compression="gzip")[
                "cube"
            ]
            cube = reshape_cube(cube)

            w_in, _ = self.load_optimal_wcs(wcscorr=wcscorr)

            # Initialize the mosaic
            comb_cube = np.zeros(cube.shape)
            sum_weights = np.zeros(cube.shape)
            comb_cube_unc = np.zeros(cube.shape)

            wave = stitch_utils.get_wave(fname_cube)
            fname_out = self.fname_mosaic_single_segment(
                seg, wcscorr=wcscorr, calibrated=calibrated, scaled=scaled
            )

        print(f"fname_out = {fname_out}")

        # pointing_labels = ['1', '3', '5', '7', '9', 'b', 'd', 'f', 'h']
        n_pointings = len(pointing_labels)
        coverage_mask = np.zeros(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2]), dtype=int
        )
        coadd_mask = np.ones(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2]), dtype=int
        )
        distance_transforms = np.zeros(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2])
        )

        print("Calculating coverage array and distance transform for each pointing...")
        for i, pointing_label in enumerate(pointing_labels):
            fname_cube = self.fname_stitched_cube_single_pointing(
                pointing_label, calibrated=calibrated, scaled=scaled
            )
            # fname_cube = f'{data_dir}NIRSpec_stitched_allchannels_pointing_{pointing_label}{tag_outlier}.fits{gz}'
            print(f"    {fname_cube}")
            with fits.open(fname_cube) as f:
                # Set up inverse-variance weights
                coverage_i = np.nansum(f["CUBE"].data, axis=0)

                coverage_i[coverage_i != 0] = 1
                dist_i = ndi.distance_transform_edt(coverage_i)
                # set 0 distance to NaN
                dist_i[dist_i == 0] = np.nan
                distance_transforms[i] = dist_i
                coverage_mask[i] = coverage_i

        print(
            "Done. Now making decisions on which data to include in regions of overlap..."
        )
        # Need a pixel index grid
        xs, ys = np.arange(coverage_mask[0].shape[0]), np.arange(
            coverage_mask[0].shape[1]
        )
        xs, ys = np.meshgrid(xs, ys, indexing="ij")
        print(coverage_mask[0].shape)

        # Want to see how coadd_mask[i] change at fixed i
        coadd_mask_dict = dict()
        for i in range(n_pointings - 1):
            coadd_mask_dict[f"iloop{i}"] = dict()
            for stage in range(0, 8):
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"] = dict()

        notch_mask = np.ones(coverage_mask[0].shape)

        for i in range(n_pointings - 1):
            print(f"    Pointing {i}")
            # Now go through again, focus on places with multiple coverage
            # overlap = (coverage_mask[i] == 1) & (coverage_mask[i + 1] == 1)
            # Worst case: both close to border
            if mask_top_left_bottom_right:
                # Just mask top left and bottom right of each pointing

                # Update coadd_mask dict
                # Stage 0
                stage = 0
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                not_already_masked = coadd_mask[i] != 0

                both_covered = (coverage_mask[i] == 1) & (coverage_mask[i + 1] == 1)

                both_close_to_edge = (
                    distance_transforms[i] <= edge_dist_thresh_pix
                ) & (distance_transforms[i + 1] <= edge_dist_thresh_pix)

                idcs_worst_case = both_covered & both_close_to_edge

                # Mask the upper left edge of i + 1
                ymin_i = np.min(ys[idcs_worst_case])
                xmax_i = np.max(xs[idcs_worst_case])
                idcs_worst_case_fix = (
                    not_already_masked
                    & both_covered
                    & both_close_to_edge
                    & (ys >= ymin_i)
                    & (ys <= ymin_i + edge_dist_thresh_pix)
                    & (xs <= xmax_i)
                    & (xs >= xmax_i - edge_dist_thresh_pix)
                )

                coadd_mask[i][idcs_worst_case_fix] = 1
                coadd_mask[i + 1][idcs_worst_case_fix] = 0
                # Update coadd_mask dict
                # Stage 1
                stage = 1
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                # Mask the bottom right edge of i
                xmin_i = np.min(xs[idcs_worst_case])
                ymax_i = np.max(ys[idcs_worst_case])
                idcs_worst_case_fix_lower = (
                    both_covered
                    & (distance_transforms[i] <= edge_dist_thresh_pix)
                    & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                    & (xs >= xmin_i)
                    & (xs <= xmin_i + edge_dist_thresh_pix)
                    & (ys <= ymax_i)
                    & (ys >= ymax_i - edge_dist_thresh_pix)
                )

                coadd_mask[i][idcs_worst_case_fix_lower] = 0
                coadd_mask[i + 1][idcs_worst_case_fix_lower] = 1
                # Update coadd_mask dict
                # Stage 2
                stage = 2
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                # Upper notch
                notch_mask[idcs_worst_case_fix] = 0
                # Lower notch
                notch_mask[idcs_worst_case_fix_lower] = 0

                # # Next case: i + 1 closer to border but i is not
                # idcs_case2 = not_already_masked & both_covered & (
                #     distance_transforms[i] >= edge_dist_thresh_pix) & (
                #         distance_transforms[i + 1] <= edge_dist_thresh_pix)
                # # Favour i
                # coadd_mask[i][idcs_case2] = 1
                # coadd_mask[i + 1][idcs_case2] = 0
                # # Update coadd_mask dict
                # # Stage 3
                # stage = 3
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}'][
                #     'i'] = coadd_mask[i].copy()
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}'][
                #     'i+1'] = coadd_mask[i + 1].copy()

                # # Clean up at the end
                # idcs_case6 = not_already_masked & both_covered & (
                #     coadd_mask[i] == 0) & (coadd_mask[i + 1] == 0) & (
                #         distance_transforms[i] <= edge_dist_thresh_pix) & (
                #             distance_transforms[i + 1] <= edge_dist_thresh_pix)
                # # Favour neither, coadd them
                # coadd_mask[i][idcs_case6] = 1
                # coadd_mask[i + 1][idcs_case6] = 1
                # # Update coadd_mask dict
                # # Stage 4
                # stage = 4
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}'][
                #     'i'] = coadd_mask[i].copy()
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}'][
                #     'i+1'] = coadd_mask[i + 1].copy()
            else:

                # Update coadd_mask dict
                # Stage 0
                stage = 0
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                not_already_masked = coadd_mask[i] != 0

                both_covered = (coverage_mask[i] == 1) & (coverage_mask[i + 1] == 1)

                both_close_to_edge = (
                    distance_transforms[i] <= edge_dist_thresh_pix
                ) & (distance_transforms[i + 1] <= edge_dist_thresh_pix)

                idcs_worst_case = both_covered & both_close_to_edge
                # If both are close to the border, favour i + 1
                coadd_mask[i][idcs_worst_case] = 0
                coadd_mask[i + 1][idcs_worst_case] = 1
                # Update coadd_mask dict
                # Stage 1
                stage = 1
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                # Except for the upper left edge of i + 1
                ymin_i = np.min(ys[idcs_worst_case])
                xmax_i = np.max(xs[idcs_worst_case])
                idcs_worst_case_fix = (
                    not_already_masked
                    & both_covered
                    & both_close_to_edge
                    & (ys >= ymin_i)
                    & (ys <= ymin_i + edge_dist_thresh_pix)
                    & (xs <= xmax_i)
                    & (xs >= xmax_i - edge_dist_thresh_pix)
                )

                coadd_mask[i][idcs_worst_case_fix] = 1
                coadd_mask[i + 1][idcs_worst_case_fix] = 0
                # Update coadd_mask dict
                # Stage 2
                stage = 2
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                # Bottom right edge of i
                xmin_i = np.min(xs[idcs_worst_case])
                ymax_i = np.max(ys[idcs_worst_case])
                idcs_worst_case_fix_lower = (
                    both_covered
                    & (distance_transforms[i] <= edge_dist_thresh_pix)
                    & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                    & (xs >= xmin_i)
                    & (xs <= xmin_i + edge_dist_thresh_pix)
                    & (ys <= ymax_i)
                    & (ys >= ymax_i - edge_dist_thresh_pix)
                )

                # Upper notch
                notch_mask[idcs_worst_case_fix] = 0
                # Lower notch
                notch_mask[idcs_worst_case_fix_lower] = 0

                # Next case: i + 1 closer to border but i is not
                idcs_case2 = (
                    not_already_masked
                    & both_covered
                    & (distance_transforms[i] >= edge_dist_thresh_pix)
                    & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                )
                # Favour i
                coadd_mask[i][idcs_case2] = 1
                coadd_mask[i + 1][idcs_case2] = 0
                # Update coadd_mask dict
                # Stage 3
                stage = 3
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

                # # Next case: i closer to border but i + 1 is not
                # idcs_case3 = not_already_masked & both_covered & (
                #     distance_transforms[i] <= edge_dist_thresh_pix) & (
                #         distance_transforms[i + 1] > edge_dist_thresh_pix)
                # # Favour i + 1
                # coadd_mask[i][idcs_case3] = 0
                # coadd_mask[i + 1][idcs_case3] = 1
                # # Update coadd_mask dict
                # # Stage 4
                # stage = 4
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i'] = coadd_mask[i].copy()
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i+1'] = coadd_mask[
                #     i + 1].copy()

                # # Next case: both close to border
                # idcs_case4 = not_already_masked & both_covered & (
                #     distance_transforms[i] <= edge_dist_thresh_pix) & (
                #         distance_transforms[i + 1] <= edge_dist_thresh_pix)
                # # Favour i + 1
                # coadd_mask[i][idcs_case4] = 0
                # coadd_mask[i + 1][idcs_case4] = 1
                # # Update coadd_mask dict
                # # Stage 5
                # stage = 5
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i'] = coadd_mask[i].copy()
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i+1'] = coadd_mask[
                #     i + 1].copy()

                # # Next case: i and i+1 far from border --> safe to coadd them
                # idcs_case5 = not_already_masked & both_covered & (
                #     distance_transforms[i] >= edge_dist_thresh_pix) & (
                #         distance_transforms[i + 1] >= edge_dist_thresh_pix)
                # # Favour neither, coadd them
                # coadd_mask[i][idcs_case5] = 1
                # coadd_mask[i + 1][idcs_case5] = 1
                # # Update coadd_mask dict
                # # Stage 6
                # stage = 6
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i'] = coadd_mask[i].copy()
                # coadd_mask_dict[f'iloop{i}'][f'stage_{stage}']['i+1'] = coadd_mask[
                #     i + 1].copy()

                # Clean up at the end
                idcs_case6 = (
                    not_already_masked
                    & both_covered
                    & (coadd_mask[i] == 0)
                    & (coadd_mask[i + 1] == 0)
                    & (distance_transforms[i] <= edge_dist_thresh_pix)
                    & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                )
                # Favour neither, coadd them
                coadd_mask[i][idcs_case6] = 1
                coadd_mask[i + 1][idcs_case6] = 1
                # Update coadd_mask dict
                # Stage 5
                stage = 4
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i"] = coadd_mask[
                    i
                ].copy()
                coadd_mask_dict[f"iloop{i}"][f"stage_{stage}"]["i+1"] = coadd_mask[
                    i + 1
                ].copy()

        if make_save_minimal_edge_mask:
            print("Making minimal edge mask")
            # Now make a minimal mask
            minimal_edge_mask = np.zeros(coadd_mask[0].shape)
            for i in range(n_pointings):
                good = (coadd_mask[i] == 1) & (coverage_mask[i] == 1)
                minimal_edge_mask[good] = 1
            minimal_edge_mask *= notch_mask

            pl.figure()
            pl.imshow(minimal_edge_mask, origin="lower")
            pl.savefig("/arc/home/rchown/mask.png")
            pl.close()

            print(f"Done. Saving to {fname_out_mask}")
            hdr0 = w_in.to_header()
            # Every fits file needs a PrimaryHDU. We'll make a blank one
            hdu0 = fits.PrimaryHDU()
            # Start an HDUList
            hdu_list = [hdu0]
            # One ImageHDU per image
            hdr = hdr0.copy()
            hdr["EXTNAME"] = ("MASK", "Mask for edge pixels")
            q = 1 * u.m / u.m  # dimensionless
            hdr["UNIT"] = f"{q.unit:FITS}"
            hdu_i = fits.ImageHDU(data=minimal_edge_mask, header=hdr)
            hdu_list = hdu_list + [hdu_i]
            hdul = fits.HDUList(hdus=hdu_list)
            hdul.writeto(fname_out_mask, overwrite=True)

        print("Done. Now coadding...")
        # Now do the coadd
        for i, pointing_label in enumerate(pointing_labels):
            if seg == "all":
                fname_cube = self.fname_stitched_cube_single_pointing(
                    pointing_label, calibrated=calibrated, scaled=scaled
                )
                print(f"    {fname_cube}")
                with fits.open(fname_cube) as f:
                    # Set up inverse-variance weights
                    # Note: multiplying a 3d arr by a 2d arr gives a 3d array with each plane multiplied by the 2d arr
                    wt_i = f["ERR"].data * coadd_mask[i]
                    cube_i = f["CUBE"].data * coadd_mask[i]
                    good = (
                        (wt_i != 0)
                        & np.isfinite(wt_i)
                        & (cube_i != 0)
                        & np.isfinite(cube_i)
                    )
                    wt_i[good] = 1.0 / wt_i[good] ** 2
                    sum_weights[good] += wt_i[good]
                    comb_cube[good] += wt_i[good] * cube_i[good]

            else:
                # Grab single segment data
                list_of_fnames = self.fnames_s3d_with_pattern(
                    [pointing_label, seg.replace("f", "")]
                )
                print(list_of_fnames)
                # Load the first cube
                fname_cube = list_of_fnames[0]
                fname_cube_reproj = self.fname_reprojected_s3d(fname_cube, False)
                cube = pickle.load(open(fname_cube_reproj, "rb"), compression="gzip")[
                    "cube"
                ]
                cube = reshape_cube(cube)

                fname_cube_unc_reproj = self.fname_reprojected_s3d(fname_cube, True)
                cube_unc = pickle.load(
                    open(fname_cube_unc_reproj, "rb"), compression="gzip"
                )["cube"]
                cube_unc = reshape_cube(cube_unc)

                # unc(a * b) = |a * b| * sqrt((unc_a / a)**2 + (unc_b / b)**2)
                cube_unc = np.abs(cube * calfac) * np.sqrt(
                    (cube_unc / cube) ** 2 + (calfac_err * calfac) ** 2
                )
                cube = cube * calfac

                wt_i = cube_unc * coadd_mask[i]
                cube_i = cube * coadd_mask[i]
                good = (
                    (wt_i != 0)
                    & np.isfinite(wt_i)
                    & (cube_i != 0)
                    & np.isfinite(cube_i)
                )
                wt_i[good] = 1.0 / wt_i[good] ** 2
                sum_weights[good] += wt_i[good]
                comb_cube[good] += wt_i[good] * cube_i[good]

        good = (sum_weights != 0) & np.isfinite(sum_weights)
        comb_cube[good] /= sum_weights[good]
        comb_cube_unc[good] = 1.0 / np.sqrt(sum_weights[good])

        print(f"Done. Now writing to {fname_out}...")
        # Write to fits
        hdr0 = w_in.to_header()
        # Every fits file needs a PrimaryHDU. We'll make a blank one
        hdu0 = fits.PrimaryHDU()
        # Start an HDUList
        hdu_list = [hdu0]
        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("CUBE", "Stitched NIRSpec datacube")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("ERR", "NIRSpec datacube uncertainty")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube_unc, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One for wavelength array
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("WAVE", "Wavelength array")
        q = 1 * u.um
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=wave, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # Make an HDUList from all HDUs we want to combine into a single FITS file
        hdul = fits.HDUList(hdus=hdu_list)

        hdul.writeto(fname_out, overwrite=True)
        print("Done.")

        return coadd_mask_dict

    def make_and_save_synthetic_images(
        self,
        fname_cube,
        fname_synth_images_fits,
        fname_synth_images_unc_fits,
        path_to_throughputs,
    ):
        """
        For NIRSpec only at the moment
        """
        sim = synth.Synth(path_to_throughputs)
        synth_im_dict, cc_dict = sim.images(fname_cube)

        w = wcs.WCS(fits.open(fname_cube)["CUBE"].header).dropaxis(2)
        synth.write_synth_to_fits(
            synth_im_dict, w, fname_synth_images_fits, fname_synth_images_unc_fits
        )
        # synth.write_synth_to_pickle(fname_synth_images_pickle)

        return synth_im_dict, cc_dict

    def make_and_save_synthetic_images_mosaic(
        self,
        fname_synth_images_fits,
        fname_synth_images_unc_fits,
        path_to_throughputs,
        seg="all",
        wcscorr=True,
        calibrated=False,
        scaled=True,
    ):
        """
        Specific case of make_and_save_synthetic_images but for
        stitched cube of all pointings.
        """
        # if path_to_throughputs is not None:
        #     self.path_to_throughputs = path_to_throughputs
        # if (self.path_to_throughputs is None) and (path_to_throughputs is
        #                                            None):
        #     print("Need to define a path to throughputs to do anything")
        # pass
        if seg == "all":
            fname_cube = self.fname_stitched_cube_all_pointings
            if calibrated:
                fname_cube = self.fname_stitched_cube_all_pointings_calibrated
            if scaled:
                fname_cube = fname_cube.replace(".fits", "_scaled.fits")
        else:
            fname_cube = self.fname_mosaic_single_segment(
                seg, wcscorr=wcscorr, calibrated=calibrated, scaled=scaled
            )

        return self.make_and_save_synthetic_images(
            fname_cube,
            fname_synth_images_fits,
            fname_synth_images_unc_fits,
            path_to_throughputs,
        )

    # Calculate wcs offsets
    def realign_nirspec_fits(
        self,
        filt,
        fname_synth_images_fits,
        fname_synth_images_unc_fits,
        search_gaia=False,
        i_gaia=1,
        calibrated=False,
        scaled=True,
    ):
        import wcs_offsets

        fname_cube = self.fname_stitched_cube_all_pointings
        if calibrated:
            fname_cube = self.fname_stitched_cube_all_pointings_calibrated
        if scaled:
            fname_cube = fname_cube.replace(".fits", "_scaled.fits")

        # Output filenames
        fname_out_cube = self.fname_stitched_cube_all_pointings_wcs_corrected
        if calibrated:
            fname_out_cube = (
                self.fname_stitched_cube_all_pointings_wcs_corrected_calibrated
            )
        if scaled:
            fname_out_cube = fname_out_cube.replace(".fits", "_scaled.fits")

        fname_out_synth = fname_fits_add_wcscorr(fname_synth_images_fits)
        fname_out_synth_unc = fname_fits_add_wcscorr(
            fname_synth_images_unc_fits
        ).replace(".fits", "_calibrated.fits")

        if scaled:
            fname_out_synth_unc = fname_out_synth_unc.replace(".fits", "_scaled.fits")

        fname_out_wcs = self.fname_wcs_corrected
        wcs_offsets.realign_nirspec_fits(
            filt,
            fname_synth_images_fits,
            fname_synth_images_unc_fits,
            fname_cube,
            fname_out_cube,
            fname_out_synth,
            fname_out_synth_unc,
            fname_out_wcs,
            i_gaia=i_gaia,
            search_gaia=search_gaia,
        )

    def reproject_nircam(self, fname_nircam, fname_out):
        wcs_out, shape_out = self.load_optimal_wcs(wcscorr=True)

        stitch_utils.reproject_nircam(fname_nircam, wcs_out, shape_out, fname_out)

    def coadd_pointings_one_segment(
        self,
        segment_label="",
        edge_dist_thresh_pix=6,
        arcade=True,
        outlier_off=False,
        which_input_dir="",
    ):
        from scipy import ndimage as ndi

        # files3d = [
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F100LP_2{pointing_group.upper()}_g140h-f100lp_s3d.fits',
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F170LP_1{pointing_group.upper()}_g235h-f170lp_s3d.fits',
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F290LP_0{pointing_group.upper()}_g395h-f290lp_s3d.fits'
        # ]

        data_dir = "/Volumes/data1/jwst/data/ers1288/from_ftp/"
        if arcade:
            output_dir = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/single_segments_12Nov22/"
            data_dir = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/"
        if outlier_off:
            # tag_outlier = '_outlier_off'
            tag_outlier = "_12Nov22"
        else:
            tag_outlier = ""

        tag1 = segment_label.split("-")[1].upper()
        if which_input_dir == "":
            fname_out = (
                f"{output_dir}Level3_{tag1}_all_pointings_{segment_label}_s3d.fits.gz"
            )
            fname_out_wcsfixed = f"{output_dir}Level3_{tag1}_all_pointings_{segment_label}_s3d_wcsfixed.fits.gz"
            path_to_stitched_cube_data = "/arc/home/rchown/stitched_cube_data/"
        else:
            # fname_out = f'/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/single_segments_7Nov22/Level3_{tag1}_all_pointings_{segment_label}_s3d.fits.gz'
            # fname_out_wcsfixed = f'/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/single_segments_7Nov22/Level3_{tag1}_all_pointings_{segment_label}_s3d_wcsfixed.fits.gz'
            # path_to_stitched_cube_data = '/arc/home/rchown/stitched_cube_data_7Nov22/'
            fname_out = f"/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/single_segments_12Nov22/Level3_{tag1}_all_pointings_{segment_label}_s3d.fits.gz"
            fname_out_wcsfixed = f"/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/single_segments_12Nov22/Level3_{tag1}_all_pointings_{segment_label}_s3d_wcsfixed.fits.gz"
            path_to_stitched_cube_data = "/arc/home/rchown/stitched_cube_data_12Nov22/"

        print(f"fname_out = {fname_out}")

        pointing_label = "1"
        fname_cube = f"{data_dir}NIRSpec_stitched_allchannels_pointing_{pointing_label}{tag_outlier}.fits"

        with fits.open(fname_cube) as f:
            w_in = wcs.WCS(f["CUBE"].header)
            comb_cube = np.zeros(f["CUBE"].data.shape)
            sum_weights = np.zeros(f["CUBE"].data.shape)
            comb_cube_unc = np.zeros(f["CUBE"].data.shape)
            wave = f["WAVE"].data

        path_to_s3d = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/"
        # files3d = glob.glob(
        #     f'{path_to_s3d}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_{tag1}_*s3d.fits'
        # )
        files3d_f100 = glob.glob(
            f"{path_to_s3d}NIRSPEC/NIRSpec_F100LP_Orion_Bar_Nov_12_2022/Pointing_*_Level3_g140h-f100lp_s3d.fits"
        )
        files3d_f170 = glob.glob(
            f"{path_to_s3d}NIRSPEC/NIRSpec_F170LP_Orion_Bar_Nov_12_2022/Pointing_*_Level3_g235h-f170lp_s3d.fits"
        )
        files3d_f290 = glob.glob(
            f"{path_to_s3d}NIRSPEC/NIRSpec_F290LP_Orion_Bar_Nov_12_2022/Pointing_*_Level3_g395h-f290lp_s3d.fits"
        )
        files3d = files3d_f100 + files3d_f170 + files3d_f290
        # files3d = glob.glob(
        #     f'{path_to_s3d}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_{tag1}_*s3d.fits'
        # )

        pointing_labels = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        n_pointings = len(pointing_labels)
        coverage_mask = np.zeros(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2]), dtype=int
        )
        coadd_mask = np.ones(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2]), dtype=int
        )
        distance_transforms = np.zeros(
            (n_pointings, comb_cube.shape[1], comb_cube.shape[2])
        )

        print("Calculating coverage array and distance transform for each pointing...")
        for i, pointing_label in enumerate(pointing_labels):
            fname_cube = f"{data_dir}NIRSpec_stitched_allchannels_pointing_{pointing_label}{tag_outlier}.fits"
            print(f"    {fname_cube}")
            with fits.open(fname_cube) as f:
                # Set up inverse-variance weights
                coverage_i = np.nansum(f["CUBE"].data, axis=0)

                coverage_i[coverage_i != 0] = 1
                dist_i = ndi.distance_transform_edt(coverage_i)
                # set 0 distance to NaN
                dist_i[dist_i == 0] = np.nan
                distance_transforms[i] = dist_i
                coverage_mask[i] = coverage_i

        print(
            "Done. Now making decisions on which data to include in regions of overlap..."
        )
        # Need a pixel index grid
        xs, ys = np.arange(coverage_mask[0].shape[0]), np.arange(
            coverage_mask[0].shape[1]
        )
        xs, ys = np.meshgrid(xs, ys, indexing="ij")
        print(coverage_mask[0].shape)

        notch_mask = np.ones(coverage_mask[0].shape)

        for i in range(n_pointings - 1):
            print(f"    Pointing {i}")
            # Now go through again, focus on places with multiple coverage
            # overlap = (coverage_mask[i] == 1) & (coverage_mask[i + 1] == 1)
            # Worst case: both close to border
            idcs_worst_case = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] <= edge_dist_thresh_pix)
                & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
            )
            # Favour i + 1
            coadd_mask[i][idcs_worst_case] = 0
            coadd_mask[i + 1][idcs_worst_case] = 1
            # Except for the upper left edge of i + 1
            ymin_i = np.min(ys[idcs_worst_case])
            xmax_i = np.max(xs[idcs_worst_case])
            idcs_worst_case_fix = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] <= edge_dist_thresh_pix)
                & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                & (ys >= ymin_i)
                & (ys <= ymin_i + edge_dist_thresh_pix)
                & (xs <= xmax_i)
                & (xs >= xmax_i - edge_dist_thresh_pix)
            )

            coadd_mask[i][idcs_worst_case_fix] = 1
            coadd_mask[i + 1][idcs_worst_case_fix] = 0

            # Bottom right edge of i
            xmin_i = np.min(xs[idcs_worst_case])
            ymax_i = np.max(ys[idcs_worst_case])
            idcs_worst_case_fix_lower = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] <= edge_dist_thresh_pix)
                & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
                & (xs >= xmin_i)
                & (xs <= xmin_i + edge_dist_thresh_pix)
                & (ys <= ymax_i)
                & (ys >= ymax_i - edge_dist_thresh_pix)
            )

            # Upper notch
            notch_mask[idcs_worst_case_fix] = 0
            # Lower notch
            notch_mask[idcs_worst_case_fix_lower] = 0

            # Next case: i + 1 closer to border but i is not
            idcs_case2 = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] > edge_dist_thresh_pix)
                & (distance_transforms[i + 1] <= edge_dist_thresh_pix)
            )
            # Favour i
            coadd_mask[i][idcs_case2] = 1
            coadd_mask[i + 1][idcs_case2] = 0

            # Next case: i closer to border but i + 1 is not
            idcs_case3 = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] <= edge_dist_thresh_pix)
                & (distance_transforms[i + 1] > edge_dist_thresh_pix)
            )
            # Favour i + 1
            coadd_mask[i][idcs_case3] = 0
            coadd_mask[i + 1][idcs_case3] = 1

            # Next case: i and i+1 far from border --> safe to coadd them
            idcs_case4 = (
                (coverage_mask[i] == 1)
                & (coverage_mask[i + 1] == 1)
                & (distance_transforms[i] >= edge_dist_thresh_pix)
                & (distance_transforms[i + 1] >= edge_dist_thresh_pix)
            )
            # Favour neither, coadd them
            coadd_mask[i][idcs_case4] = 1
            coadd_mask[i + 1][idcs_case4] = 1

        pointing_labels = ["1", "3", "5", "7", "9", "b", "d", "f", "h"]
        n_pointings = len(pointing_labels)
        if "f100lp" in segment_label:
            files3d = files3d_f100
        if "f170lp" in segment_label:
            files3d = files3d_f170
        if "f290lp" in segment_label:
            files3d = files3d_f290

        # if which_input_dir != '':
        #     path_to_s3d = '/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/'
        #     files3d = glob.glob(
        #         f'{path_to_s3d}NIRSPEC/NIRSpec_{filter_label}_Orion_Bar_{date}/Pointing_*s3d.fits'
        #     )

        wave = stitch_utils.get_wave(files3d[0])

        # fnames_cubes = [
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F100LP_2{which_pointing.upper()}_g140h-f100lp_s3d.fits',
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F170LP_1{which_pointing.upper()}_g235h-f170lp_s3d.fits',
        #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F290LP_0{which_pointing.upper()}_g395h-f290lp_s3d.fits'
        # ]
        unc_str = ""
        name_cube_pickle = files3d[0].split("/")[-1].split(".")[0] + unc_str + ".pk.gz"
        cube = pickle.load(open(path_to_stitched_cube_data + name_cube_pickle, "rb"))[
            "cube"
        ]
        nx, ny, nw = cube.shape
        comb_cube = np.zeros((nw, nx, ny))
        comb_cube_unc = np.zeros((nw, nx, ny))
        sum_weights = np.zeros((nw, nx, ny))

        # with fits.open(fname_cube) as f:
        #     w_in = wcs.WCS(f['CUBE'].header)
        #     comb_cube = np.zeros(f['CUBE'].data.shape)
        #     sum_weights = np.zeros(f['CUBE'].data.shape)
        #     comb_cube_unc = np.zeros(f['CUBE'].data.shape)
        #     wave = f['WAVE'].data

        # Now do the coadd
        for i, pointing_label in enumerate(pointing_labels):
            # fnames_cubes = [
            #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F100LP_2{which_pointing.upper()}_g140h-f100lp_s3d.fits',
            #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F170LP_1{which_pointing.upper()}_g235h-f170lp_s3d.fits',
            #     f'{path_to_data}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_F290LP_0{which_pointing.upper()}_g395h-f290lp_s3d.fits'
            # ]
            if which_input_dir == "":
                path_to_s3d = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/"
                if tag1 == "F100LP":
                    j = 2
                if tag1 == "F170LP":
                    j = 1
                if tag1 == "F290LP":
                    j = 0

                # files3d = glob.glob(f'{path_to_s3d}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_{tag1}_{j}{pointing_label}_*s3d.fits')
                # fname_cube = f'{data_dir}NIRSpec_stitched_allchannels_pointing_{pointing_label}{tag_outlier}.fits.gz'
                fname_cube = glob.glob(
                    f"{path_to_s3d}NIRSPEC/NIRSpec_all_Orion_Bar_Oct_14_2022/Level3_{tag1}_{j}{pointing_label.upper()}_*s3d.fits"
                )[0]
                print(f"    {fname_cube}")
                unc_str = ""
                name_cube_pickle = (
                    fname_cube.split("/")[-1].split(".")[0] + unc_str + ".pk"
                )
                unc_str = "_unc"
                name_cube_unc_pickle = (
                    fname_cube.split("/")[-1].split(".")[0] + unc_str + ".pk"
                )
            else:
                # Pointing_1_Level3_g140h-f100lp_s3d.fits
                path_to_s3d = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/"

                if "f100lp" in segment_label:
                    fname_cube = glob.glob(
                        f"{path_to_s3d}NIRSPEC/NIRSpec_F100LP_Orion_Bar_Nov_12_2022/Pointing_{pointing_label}_Level3_g140h-f100lp_s3d.fits"
                    )[0]
                if "f170lp" in segment_label:
                    fname_cube = glob.glob(
                        f"{path_to_s3d}NIRSPEC/NIRSpec_F170LP_Orion_Bar_Nov_12_2022/Pointing_{pointing_label}_Level3_g235h-f170lp_s3d.fits"
                    )[0]
                if "f290lp" in segment_label:
                    fname_cube = glob.glob(
                        f"{path_to_s3d}NIRSPEC/NIRSpec_F290LP_Orion_Bar_Nov_12_2022/Pointing_{pointing_label}_Level3_g395h-f290lp_s3d.fits"
                    )[0]

                # fname_cube = glob.glob(
                #     f'{path_to_s3d}NIRSPEC/NIRSpec_{filter_label}_Orion_Bar_{date}/Pointing_{pointing_label.lower()}_*s3d.fits'
                # )[0]
                print(f"    {fname_cube}")
                unc_str = ""
                name_cube_pickle = (
                    fname_cube.split("/")[-1].split(".")[0] + unc_str + ".pk"
                )
                unc_str = "_unc"
                name_cube_unc_pickle = (
                    fname_cube.split("/")[-1].split(".")[0] + unc_str + ".pk"
                )

            # Set up inverse-variance weights
            # Note: multiplying a 3d arr by a 2d arr gives a 3d array with each plane multiplied by the 2d arr
            # Fluxes

            cube = pickle.load(
                open(path_to_stitched_cube_data + name_cube_pickle, "rb"),
                compression="gzip",
            )["cube"]

            cube = reshape_cube(cube)
            # Uncertainties
            cube_unc = pickle.load(
                open(path_to_stitched_cube_data + name_cube_unc_pickle, "rb"),
                compression="gzip",
            )["cube"]

            cube_unc = reshape_cube(cube_unc)

            wt_i = cube_unc * coadd_mask[i]
            cube_i = cube * coadd_mask[i]
            good = (wt_i != 0) & np.isfinite(wt_i) & (cube_i != 0) & np.isfinite(cube_i)
            wt_i[good] = 1.0 / wt_i[good] ** 2
            sum_weights[good] += wt_i[good]
            comb_cube[good] += wt_i[good] * cube_i[good]

        good = (sum_weights != 0) & np.isfinite(sum_weights)
        comb_cube[good] /= sum_weights[good]
        comb_cube_unc[good] = 1.0 / np.sqrt(sum_weights[good])

        print(f"Done. Now writing to {fname_out}...")
        # Write to fits
        hdr0 = w_in.to_header()
        # Every fits file needs a PrimaryHDU. We'll make a blank one
        hdu0 = fits.PrimaryHDU()
        # Start an HDUList
        hdu_list = [hdu0]
        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("CUBE", "Stitched NIRSpec datacube")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("ERR", "NIRSpec datacube uncertainty")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube_unc, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One for wavelength array
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("WAVE", "Wavelength array")
        q = 1 * u.um
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=wave, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # Make an HDUList from all HDUs we want to combine into a single FITS file
        hdul = fits.HDUList(hdus=hdu_list)
        hdul.writeto(fname_out, overwrite=True)

        # Now write out the version with WCS fixed
        data_dir = "/arc/home/rchown/jwst-1288-ftp.ias.u-psud.fr/data/jwst/full_archive/OrionBar/NIRSPEC/stitched_cubes/"
        # if outlier_off:
        #     tag_outlier = '_outlier_off'
        # else:
        #     tag_outlier = ''
        fname_cube_ref = f"{data_dir}NIRSpec_stitched_allchannels_all_pointings_wcsfixed{tag_outlier}.fits.gz"
        print("Using reference cube to save wcs-corrected cube:")
        print(fname_cube_ref)
        print("Saving to")
        print(fname_out_wcsfixed)
        w_fixed = wcs.WCS(fits.open(fname_cube_ref)[1].header).dropaxis(2)

        # Write to fits
        hdr0 = w_fixed.to_header()
        # Every fits file needs a PrimaryHDU. We'll make a blank one
        hdu0 = fits.PrimaryHDU()
        # Start an HDUList
        hdu_list = [hdu0]
        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("CUBE", "Stitched NIRSpec datacube")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("ERR", "NIRSpec datacube uncertainty")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=comb_cube_unc, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One for wavelength array
        hdr = hdr0.copy()
        hdr["EXTNAME"] = ("WAVE", "Wavelength array")
        q = 1 * u.um
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=wave, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # Make an HDUList from all HDUs we want to combine into a single FITS file
        hdul = fits.HDUList(hdus=hdu_list)
        hdul.writeto(fname_out_wcsfixed, overwrite=True)

    def measure_xcal(
        self,
        filt,
        fname_synth_images_fits,
        fname_synth_images_unc_fits,
        fname_nircam,
        return_xc=True,
        edge_dist_thresh_pix=9,
        region_mask=None,
        snr_thresh=0.0,
    ):
        """
        filter
        synth im
        synth im unc
        nircam im
        nircam im unc
        """
        import xcal

        fname_mask = self.fname_mask_stitched_cube_all_pointings(edge_dist_thresh_pix)

        mask = fits.open(fname_mask)["MASK"].data

        # fnames_nircam
        # fname_cal_dict

        synth_im = fits.open(fname_synth_images_fits)[filt].data
        synth_im_unc = fits.open(fname_synth_images_unc_fits)[filt].data
        w = wcs.WCS(fits.open(fname_synth_images_fits)[filt].header)
        nircam_im = fits.open(fname_nircam)["SCI"].data
        nircam_im_unc = fits.open(fname_nircam)["ERR"].data

        xmax_limit = 5000
        if filt == "F182M":
            xmax_limit = 500
        if filt == "F335M":
            xmax_limit = 1000
        if filt == "F210M":
            xmax_limit = 1000
        res_xc = xcal.measure_cal_from_nircam_nirspec(
            synth_im,
            synth_im_unc,
            w,
            nircam_im,
            nircam_im_unc,
            mask,
            region_mask=region_mask,
            mask_star=True,
            snr_thresh=snr_thresh,
            xmax_limit=xmax_limit,
        )

        # xc = xcal.XCal(fname_mosaic, fnames_nircam, fname_cal_dict)
        # res_xc = xc.measure_xcal(fname_synth_images_fits,
        #                          fname_synth_images_unc_fits, filt, fname_mask,
        #                          region_mask)

        self.xcal_dict = res_xc

        if return_xc:
            return res_xc
