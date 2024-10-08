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
from pdrs4all import synth
from pdrs4all.postprocess import spectral_segments, cube
from astropy.wcs import WCS
from pathlib import Path


def parse_args():
    ap = ArgumentParser()
    ap.add_argument(
        "nirspec_cubes",
        nargs=3,
        help="nirspec _s3d.fits of the same spatial dimensions, in wavelength order, i.e. 100, 170, 290",
    )
    ap.add_argument("--output_dir", default="./")
    args = ap.parse_args()
    return args


def main(args):
    # step 0: load individual cubes and make naive stitched cube to be
    # used for synthetic photometry
    output_path = Path(args.output_dir)

    # sort cubes by shortest wavelength
    s3ds = [Spectrum1D.read(fn) for fn in args.nirspec_cubes]
    s3ds.sort(key=lambda x: x.spectral_axis.value[0])

    # use algorithm from Dries' package. We can include a copy later.
    nirspec_cwcs = WCS(s3ds[0].meta["header"]).celestial
    s3dm = spectral_segments.merge_nd(s3ds)
    cube.write_cube_s1d_wavetab_jwst_s3d_format(
        output_path / "nirspec_naive_stitch.fits",
        s3dm,
        nirspec_cwcs,
    )

    # step 1a: nirspec synthetic photometry
    synth_image_dict, cc_dict = synth.synthesize_nircam_images(s3dm)
    # write them out to see if its working
    synth.write_synth_to_fits(
        synth_image_dict,
        nirspec_cwcs,
        output_path / "nirspec_synth_nircam.fits",
        output_path / "nirspec_synth_nircam_unc.fits",
    )


if __name__ == "__main__":
    main(parse_args())
