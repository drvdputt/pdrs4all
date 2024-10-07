"""Correct NIRSpec WCS and calibration based on NIRCam data.

This implements the worflow by Ryan to improve the NIRSpec data. By
comparing wavelength-integrated NIRSpec cubes to reprojected NIRCam
images, the WCS and flux calibration are improved. At the end, the
improved products are sitched. Here are main steps, with some
interdependency.

0. Naive stitch.

   To make the synthetic photometry in the next step work, we will need
   to make an initial stitched cube first, so that the integration can
   happen without problems.

1. Create comparable maps

   a. NIRSpec cubes x NIRCam filter curves -> synthetic photometry (synphot)
   b. NIRCam images x NIRSpec WCS -> NIRCam maps on same grid (ncrpj)

2. WCS correction

   a. Measure centroid of proplyd in synphot -> RA, Dec offset
   b. Apply to cubes and synphot -> wcscorr cubes

3. Calibration factor

   a. ncrpj x synphot -> factors
   b. wcscorr cubes x factors -> calcubes (See Peeters 2024 to see which ones to use)

4. Stitch end results

Ryans has a lot of workarounds (e.g. masking certain pixels in the cube)
and variations (single pointing vs mosaic vs single segment vs stitched)
in his code that make it hard to understand what works with what. I
think it will be better to reimplement some things, make sure that all
the output is in a format that I understand, and copy only the parts
that I need.

"""

from pdrs4all.stitcher import Stitcher
from pdrs4all.synth import Synth

def step_1a_make_and_save_synthetic_images(fname_cube,
                                       fname_synth_images_fits,
                                       fname_synth_images_unc_fits,
                                       path_to_throughputs):
        '''
        For NIRSpec only at the moment
        '''
        sim = Synth(path_to_throughputs)
        synth_im_dict, cc_dict = sim.images(fname_cube)

        w = wcs.WCS(fits.open(fname_cube)['CUBE'].header).dropaxis(2)
        synth.write_synth_to_fits(synth_im_dict, w, fname_synth_images_fits,
                                  fname_synth_images_unc_fits)
        # synth.write_synth_to_pickle(fname_synth_images_pickle)

        return synth_im_dict, cc_dict


if __name__ == "__main__":
    stest = Stitcher(fnames_s3d_f100lp=fnames_f100lp,
                        fnames_s3d_f170lp=fnames_f170lp,
                        fnames_s3d_f290lp=fnames_f290lp,
                        output_dir_wcs=output_dir_wcs,
                        output_dir_reprojected_s3d=output_dir_reprojected_s3d,
                        output_dir_stitched_cubes=output_dir_stitched_cubes)

    # run stitch_arcade_1overf_corr_10May.py step3_synth_full_mosaic
    stest.make_and_save_synthetic_images_mosaic(
        fname_synth_images_fits,
        fname_synth_images_unc_fits,
        path_to_throughputs=path_to_throughputs,
        wcscorr=False,
        calibrated=False,
        scaled=False)

