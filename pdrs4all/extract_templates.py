"""Script to do the extraction, given the template file and a list of
cubes.

"""

import argparse
from astropy.table import Table
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from astropy.nddata import StdDevUncertainty
from regions import Regions, SkyRegion
from specutils import Spectrum1D
from myastro import spectral_segments
import numpy as np
from matplotlib import pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "region_file",
        help="File containing DS9 regions representing template apertures.",
    )
    ap.add_argument(
        "cube_files", nargs="+", help="Data cubes to extract from and merge."
    )
    ap.add_argument(
        "--apply_offsets",
        action="store_true",
        help="""Apply additive offsets to make the spectral stitching
        smoother. Resulting continuum might be unrealistic.""",
    )
    ap.add_argument(
        "--reference_segment",
        type=int,
        help="""Index of the cube to use as a reference for spectral stitching""",
    )
    ap.add_argument(
        "--template_names",
        nargs="+",
        help="""Optional list of names to give to the template spectra.
        Number of arguments must equal the number of apertures in the given
        region file. E.g. "HII", "Atomic", "DF3", "DF2", "DF1" for the Orion
        templates.""",
    )
    args = ap.parse_args()

    # load the cubes
    cubes = spectral_segments.sort(
        [Spectrum1D.read(cf, format="JWST s3d") for cf in args.cube_files]
    )

    # set up template apertures and names
    regions = Regions.read(args.region_file)
    # apertures = [regionhacks.skyregion_to_aperture_auto(r) for r in regions]

    # determine template names
    if args.template_names is None:
        template_names = [f"T{i}" for i in range(1, len(regions) + 1)]
    else:
        template_names = args.template_names

    t = extract_templates_table(
        cubes, regions, template_names, args.apply_offsets, args.reference_segment
    )
    fname = "templates.ecsv"
    print(f"Writing extracted spectra to {fname}")

    # add some info about which files these templates were generated from
    t.meta["stitch_method"] = "additive" if args.apply_offsets else "none"
    t.meta["stitch_reference_segment"] = (
        args.reference_segment if args.apply_offsets else "none"
    )
    t.meta["cubes"] = args.cube_files

    t.write(fname, overwrite=True)


# define this local utility function
def extract_and_merge(cubes, aperture, apply_offsets, offset_reference=0):
    """Steps that need to happen for every aperture.

    1. extract from every given cube
    2. apply stitching corrections
    3. return a single merged spectrum"""
    specs = [cube_sky_aperture_extraction_v3(s, aperture) for s in cubes]

    if apply_offsets:
        shifts = spectral_segments.overlap_shifts(specs)
        # choose a segment in the middle as the reference for the stitching
        offsets = spectral_segments.shifts_to_offsets(shifts, offset_reference)
        specs_to_merge = [s + o for s, o in zip(specs, offsets)]
    else:
        specs_to_merge = specs

    return spectral_segments.merge_1d(specs_to_merge)


def extract_templates_table(
    cubes, apertures, template_names, apply_offsets=False, offset_reference=0
):
    templates = {
        k: extract_and_merge(cubes, a, apply_offsets, offset_reference)
        for k, a in zip(template_names, apertures)
    }

    # Construct astropy table and save as ECSV
    columns = {
        "wavelength": templates[template_names[0]].spectral_axis.to(u.micron),
    }
    for k, v in templates.items():
        columns[f"flux_{k}"] = v.flux
        columns[f"unc_{k}"] = v.uncertainty.array * v.flux.unit

    t = Table(columns)
    return t


def cube_sky_aperture_extraction_v3(
    cube_spec1d, sky_region: SkyRegion, average_per_sr=True, wcs_2d=None
):
    """Extract spectrum cube using aperture in sky coordinates.

    v3 uses a regions.Skyregion as a parameter, instead of a
    photutils.SkyAperture. The newer regions package has duplicated some
    of the internal functionality of photutils. And it turns out it
    better suits my needs, for the following reasons.

    - There is an example of how to multiply the masks in the
      documentation, meaning that this use case is supported.

    - I do not have to convert between SkyRegion and SkyAperture. So I
      no longer need this annoying boilerplate (per type of region) and
      will be safer from bugs (e.g. the definition of the position angle
      of a region was different in some cases).

    Parameters
    ----------
    cube_spec1d: Spectrum1D
        Typically an object loaded from a cube file using
        Spectrum1D.read. Can be a custom Spectrum1D object, but needs to
        have the right header stored in cube_spec1d.meta["header"] so
        that a WCS can be derived. Alternatively, the wcs_2d parameter
        should be used to pass a celestial WCS manually.

    sky_region: SkyRegion
        The region to use as an aperture.

    Returns
    -------
    spectrum: Spectrum1D
        The collapsed spectrum.

    """
    # make 2D wcs of cube if needed
    if wcs_2d is None:
        the_wcs_2d = WCS(cube_spec1d.meta["header"]).celestial
    else:
        the_wcs_2d = wcs_2d

    nx, ny = cube_spec1d.shape[:2]

    pixel_region = sky_region.to_pixel(the_wcs_2d)
    aperture_mask = pixel_region.to_mask()

    slices_large, slices_small = aperture_mask.get_overlap_slices((ny, nx))
    yx_slc = slices_large
    weights = aperture_mask.data
    if slices_small is None:
        print("No overlap between aperture and data!")

    # Cut out the relevant data from the cube, and set up masked array
    # to ignore some things.. We also create a handy 3D mask that we can
    # multiply the other arrays with, to make sure we ignore the same
    # pixels for them.
    cube_ma_yx = np.ma.masked_invalid(np.swapaxes(cube_spec1d.flux.value, 0, 1))
    cube_cutout = cube_ma_yx[yx_slc]
    cube_cutout_used = np.where(cube_cutout.mask, 0, 1)

    # plt.subplot(311)
    # plt.imshow(aperture_mask.to_image((ny, nx)))
    # plt.subplot(312)
    # plt.imshow(cube_cutout.sum(axis=-1))
    # plt.subplot(313)
    # plt.imshow(aperture_mask.data)
    # plt.show()

    # # Do the sum, broadcasting the mask over the slices. Einstein
    # summation for clarity. Basically a vectorized version of aperture_mask.multiply
    sum_per_slice = np.einsum(
        "xyw,xy->w", np.where(cube_cutout_used, cube_cutout, 0), weights
    )

    # Add units again, and convert to MJy
    one_px_area = (proj_plane_pixel_area(the_wcs_2d) * u.deg**2).to(u.sr)
    flux = sum_per_slice * cube_spec1d.flux.unit * one_px_area

    # Apply this to the variance array when uncertainty is provided
    if cube_spec1d.uncertainty is not None:
        uncertainty_cutout = np.swapaxes(cube_spec1d.uncertainty.array, 0, 1)[yx_slc]
        variance_cutout = np.square(np.where(cube_cutout_used, uncertainty_cutout, 0))
        varsum_per_slice = np.einsum("xyw,xy->w", variance_cutout, weights)
        sigma_flux = np.sqrt(varsum_per_slice) * cube_spec1d.flux.unit * one_px_area
    else:
        sigma_flux = None

    if average_per_sr:
        # To compute an average, we need the total weight per slice.
        # Weight * pixel size = total area. Num pixels used = total
        # weight but without bad pixels -> slightly different per slice
        weight_per_slice = np.einsum("xyw,xy->w", cube_cutout_used, weights)
        area_per_slice = weight_per_slice * one_px_area
        flux = (flux / area_per_slice).to(cube_spec1d.flux.unit)

        if sigma_flux is not None:
            sigma_flux = (sigma_flux / area_per_slice).to(cube_spec1d.flux.unit)

    if sigma_flux is not None:
        sigma_flux = StdDevUncertainty(sigma_flux)

    # make a spectrum1d object for convenience
    s1d = Spectrum1D(
        spectral_axis=cube_spec1d.spectral_axis,
        flux=flux,
        uncertainty=sigma_flux,
    )
    return s1d


if __name__ == "__main__":
    main()
