"""Command line tool to run NSClean with the right settings

Takes 1 file and 1 output name on the command line. This makes it easy
to parallellize the NSClean runs using GNU parallel.

The mask required by NSClean is automatically selected depending on the
grating. For medium resolution, custom masks are used. For high
resolution, those included in NSClean are used.

"""
import argparse
import numpy as np
from astropy.io import fits
import nsclean
from importlib.resources import files
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Wrapper around NSClean")
    parser.add_argument("rate_in")
    parser.add_argument("rate_out")
    args = parser.parse_args()

    # keep the input file open, since we want to write everything out
    # again after doing our edits
    all_hdus = fits.open(args.rate_in)
    H0 = all_hdus[0].header
    D1 = np.array(all_hdus[1].data)

    # get some info about the exposure
    grating = H0["GRATING"].lower()
    if "nrs1" in args.rate_in:
        detector = "nrs1"
    elif "nrs2" in args.rate_in:
        detector = "nrs2"
    else:
        print("detector nrs1 or nrs2 not found in file name")
        exit()

    # choose right mask here. For medium gratings, use masks by Karl M.
    # (nrs1 only)
    if grating[-1] == "m":
        mask_fn = f"nrs1_ifu_mask_thorough_update_01192_{grating}.fits.gz"
    elif grating[-1] == "h":
        mask_fn = f"{detector}_ifu_mask_thorough.fits.gz"
    else:
        print("grating can only be m or h")
        exit()

    mask_path = str(files("pdr_reduction.nsclean_masks").joinpath(mask_fn))
    with fits.open(mask_path) as hdul:
        M = np.array(hdul[0].data, dtype=np.bool_)

    # nsclean hates nans
    where_nan = np.where(np.isnan(D1))

    # replace the reference pixel nans at the edges
    D1[:3, :] = 0
    D1[2044:, :] = 0
    D1[:, :3] = 0
    D1[:, 2044:] = 0

    # Substitute the remaining nans with local medians. I've seen other
    # people just replace them with 0.
    b = 3  # box size
    for i, j in zip(*np.where(np.isnan(D1))):
        D1[i, j] = np.nanmedian(
            D1[
                max(0, i - b) : min(D1.shape[0], i + b),
                max(0, i - b) : min(D1.shape[1], i + b),
            ]
        )

    # if still nans after this, use 0 instead
    D1 = np.nan_to_num(D1)

    cleaner = nsclean.NSClean(detector.upper(), M)
    D1 = cleaner.clean(D1, buff=True)

    # reinstate the nans and store the result in the all_hdus
    D1[where_nan] = np.nan
    all_hdus[1].data = D1

    # add a comment to the header
    H0["comment"] = "Processed by NSClean Rev. " + nsclean.__version__
    H0["NSCLEANV"] = nsclean.__version__
    all_hdus[0].header = H0
    all_hdus.info()  # info before we write
    Path(args.rate_out).parent.mkdir(parents=True, exist_ok=True)
    all_hdus.writeto(args.rate_out, overwrite=True)
    all_hdus.close()


if __name__ == "__main__":
    main()
