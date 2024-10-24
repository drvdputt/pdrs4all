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

from specutils import Spectrum1D
from argparse import ArgumentParser
from pdrs4all.postprocess import spectral_segments, write_cube, synth, wcscorr
from astropy.wcs import WCS
from astropy.io import fits
from pathlib import Path
from jwst import datamodels


def parse_args(argv=None):
    ap = ArgumentParser()
    ap.add_argument(
        "nirspec_cubes",
        nargs=3,
        help="nirspec _s3d.fits of the same spatial dimensions, in wavelength order, i.e. 100, 170, 290",
    )
    ap.add_argument("--output_dir", default="./")
    args = ap.parse_args(argv)
    return args


def main(args):
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)

    print("Step 0: load cubes and make naive stitched cube")
    s3ds = [Spectrum1D.read(fn) for fn in args.nirspec_cubes]
    # merge algorithm needs sorted list of Spectrum1D cubes
    s3ds.sort(key=lambda x: x.spectral_axis.value[0])
    nirspec_cwcs = WCS(s3ds[0].meta["header"]).celestial

    s3dm = spectral_segments.merge_nd_memfriendly(s3ds)
    naive_stitch_fits = output_path / "nirspec_naive_stitch.fits"
    write_cube.write_cube_s1d_wavetab_jwst_s3d_format(
        naive_stitch_fits, s3dm, nirspec_cwcs
    )

    print("Step 1a: nirspec synthetic photometry to nircam fiters")
    synth_image_dict, cc_dict = synth.synthesize_nircam_images(s3dm)
    # write them out to see if its working
    synth_fits = output_path / "nirspec_synth.fits"
    synth_unc_fits = output_path / "nirspec_synth_unc.fits"
    synth.write_synth_to_fits(
        synth_image_dict,
        nirspec_cwcs,
        synth_fits,
        synth_unc_fits,
    )

    print("Step 2a: compute offset")

    # open the synth image fits file we just wrote out, get one image to
    # use for nirspec wcscorr. Alternatively, a single slice or
    # collapsed range of the cube could be used.
    with fits.open(synth_fits) as hdul:
        synth_F212N = hdul["F212N"].data

    new_cwcs = wcscorr.nirspec_wcscorr_using_proplyd(synth_F212N, nirspec_cwcs)

    print("Step 2b: apply offset to cube")

    # tested this in iPython, and I think it's going to work.
    naive_cube_dm = datamodels.open(naive_stitch_fits)
    naive_cube_dm.meta.wcsinfo.crval1 = new_cwcs.wcs.crval[0]
    naive_cube_dm.meta.wcsinfo.crval2 = new_cwcs.wcs.crval[1]
    naive_cube_dm.write(output_path / "nirspec_naive_stitch_wcscorr.fits")

    # alternatively, we could open the individual cubes, and correct all
    # three of their WCS


if __name__ == "__main__":
    main(parse_args())
