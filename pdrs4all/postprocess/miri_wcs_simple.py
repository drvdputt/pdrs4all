from argparse import ArgumentParser
from pathlib import Path
from specutils import Spectrum1D
from astropy.wcs import WCS
from pdrs4all.postprocess import custom_io, wcscorr

def main():
    ap = ArgumentParser()
    ap.add_argument(
        "miri_cubes",
        nargs=12,
        help="All 12 MRS _s3d.fits files (order doesn't matter)",
    )
    ap.add_argument("--output_dir", default="./")
    args = ap.parse_args()
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)

    # every channel has some PAH or dust emission in it. Let's just
    # integrate them and see if the result is usable to determine the
    # proplyd centroid

    # load cubes and sort by wavelength
    s3ds = [Spectrum1D.read(fn) for fn in args.miri_cubes]
    s3ds.sort(key=lambda s: s.spectral_axis.value[0])

    # get all celestial WCS
    cwcss = [WCS(s.meta["header"]).celestial for s in s3ds]

    # collapse cubes along wavelength direction. (Not sure what the unit
    # of this is. Proper way would be to convert flux from MJy/sr to a
    # per-wavelength unit, then integrate)
    images = [s.collapse('sum', axis='spectral') for s in s3ds]
    
    print("Writing collapsed MRS cubes to image files")
    for i in range(len(images)):
        custom_io.write_i2d(f"mrs_collapsed{i}.fits", images[i], cwcss[i])

    # problem from these written images: the proplyd is not visible in
    # the two reddest channels. IDEA: I only need to find the offsets
    # for the first channel (SHORT MEDIUM and LONG), since all four
    # channels are observed simultaneously.
    new_cwcss = wcscorr.mrs_wcscorr_using_proplyd(images, cwcss)

    print("Writing wcs corrected images as test")
    for i in range(len(images)):
        custom_io.write_i2d(f"mrs_collapsed{i}_wcscorr.fits", images[i], new_cwcss[i])
    
if __name__ == "__main__":
    main()
    
