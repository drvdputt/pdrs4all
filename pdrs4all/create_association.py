import glob
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from astropy.io import fits
from argparse import ArgumentParser


def main():
    """Command line tool to create Level 3 associations"""
    ap = ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("files", nargs="+")
    ap.add_argument("-b", "--background_files", nargs="+", default=None)
    args = ap.parse_args()

    fn = f"{args.name}_asn.json"
    if args.background_files is not None:
        other_exptypes = {"background": args.background}
    else:
        other_exptypes = None

    writeasn(args.files, fn, f"Level3_{args.name}", other_exptypes=other_exptypes)
    print("Wrote ", fn)


def sort_files_per_filter(files):
    """Define a convenience function that will split the input _cal.fits
    files into their corresponding filter or channel/band"""
    typ, filename = [], []
    for f in files:
        # using datamodels is too slow, or sometimes even gets stuck in
        # infinite recursion
        header = fits.getheader(f, ext=0)
        filt = header["FILTER"] if "FILTER" in header else None
        chan = header["CHANNEL"] if "CHANNEL" in header else None
        band = header["BAND"] if "BAND" in header else None
        if filt is not None:
            typ.append(filt)
        elif chan is not None and band is not None:
            typ.append(chan + band)
        filename.append(f)
    dic = {}
    for i, j in zip(typ, filename):
        dic.setdefault(i, []).append(j)
    return dic


def glob_and_make_per_filter_dict(directory, pattern):
    """Utility that automatically returns None (no files requested), or {}
    (no files found).
    """
    if directory is None:
        return None
    else:
        files = list(glob.glob(directory + "/" + pattern))
        if len(files) > 0:
            return sort_files_per_filter(files)
        else:
            print(f"{pattern} files not found in {directory}")
            return {}


def writeasn(files, asnfile, prodname, other_exptypes=None):
    """
    Make ASN from file list, add backgrounds/imprints, write to json.

    Parameters
    ----------

    files: list of str
        The science exposures. These will get exptype 'science'

    asnfile: str
        Name for output json file containing the ASN

    prodname: str
        Product name. Must contain "Level2" or "Level3".

    other_exptypes: dictionary {'exptype': ['file1', 'file2'], ...}.
        Can be used to specify backgrounds, imprints, etc. These will
        get exptypes as specified in the keys of the given dictionary.

    """
    print("writing asn with files", files)

    if isinstance(files, str):
        raise ValueError("files should be an iterable of string, not a string!")

    if "Level2" in prodname:
        asn = afl.asn_from_list(files, rule=DMSLevel2bBase, product_name=prodname)

    elif "Level3" in prodname:
        asn = afl.asn_from_list(files, rule=DMS_Level3_Base, product_name=prodname)

    # Replaced old code with that from the miri flight notebook
    # (https://github.com/STScI-MIRI/MRS-ExampleNB/blob/main/Flight_Notebook1/MRS_FlightNB1.ipynb)

    # Add the same background or imprint files to every product in the
    # association (every pointing with same filter gets same background)
    if other_exptypes is not None:
        for product in asn["products"]:
            for exptype, extra_files in other_exptypes.items():
                for f in extra_files:
                    product["members"].append({"expname": f, "exptype": exptype})

    _, serialized = asn.dump()
    with open(asnfile, "w") as outfile:
        outfile.write(serialized)


def create_asn(
    directory,
    pattern,
    level,
    backdir=None,
    impdir=None,
    spectroscopy=False,
    per_pointing=False,
    output_dir="./",
):
    """Creates associations for running stage 2 or 3.

    This function focuses on FINDING and SORTING the right files.
    Grouping the discovered files into associations happens in
    create_asn_per_filter.

    Parameters
    ----------
    directory, pattern : str, str
        The science files are found in the given directory, using the given
        glob pattern. Extra files and backgrounds will be used in the right
        way if the directories where they (and only they) reside are given.

    level : 2 or 3
        The number of association files depends on the level and the
        contents. For level 2, one asn is made per exposure. For level
        3, one is made per filter.

    backdir : str
        Finds and adds background files from the given background dir.
        It is assumed that all files of the right type (cal or x1d) in
        this directory are the right backgrounds. For stage 2 imaging
        with, this function will look _cal files; for stage 3
        spectroscopy it will look for _x1d files. The background are
        automatically sorted per filter. All backgrounds with matching
        filter, are added to each association created.

    impdir : str
        Finds and adds nirspec imprint files from the given imprint dir.
        For stage 2 nirspec only. This function will look for _rate
        files in the given directory. Each ASN created will consist of
        one science exposure and one imprint. The science and imprint
        pairs are matched based on alphabetical sorting.

    spectroscopy : bool
        Easy way to tell if we're doing spectroscopy, without having to
        open a file and check the instrument setting. If True, look for
        x1d background files. If False, look for cal files instead.

    output_dir : str
        ASN files will be written to this directory.

    """
    sortfiles_dict = glob_and_make_per_filter_dict(directory, pattern)
    bkg_dict = None
    imp_dict = None

    if level == 2:
        # imaging or spectroscopy
        bkg_pattern = "stage1/*rate.fits"
    elif level == 3 and spectroscopy:
        # only for spectroscopy (master background)
        bkg_pattern = "stage2/*x1d.fits"

    # print(backdir)
    # print(bkg_pattern)

    asnlist = []
    if level == 2:
        if backdir is not None:
            bkg_dict = glob_and_make_per_filter_dict(backdir, bkg_pattern)
            print("Found backgrounds", bkg_dict)
        if impdir is not None:
            imp_dict = glob_and_make_per_filter_dict(impdir, "stage1/*rate.fits")

        asnlist = create_asn_per_filter(
            sortfiles_dict,
            2,
            bkg_sortfiles_dict=bkg_dict,
            imp_sortfiles_dict=imp_dict,
            output_dir=output_dir,
        )

    if level == 3:
        if backdir is not None:
            bkg_dict = glob_and_make_per_filter_dict(backdir, bkg_pattern)

        asnlist = create_asn_per_filter(
            sortfiles_dict,
            3,
            bkg_sortfiles_dict=bkg_dict,
            output_dir=output_dir,
            per_pointing=per_pointing,
        )

    if len(asnlist) == 0:
        raise RuntimeError("No ASN's created! Something must be wrong.")

    return asnlist


def create_asn_per_filter(
    sortfiles_dict,
    level,
    bkg_sortfiles_dict=None,
    imp_sortfiles_dict=None,
    per_pointing=False,
    output_dir="./",
):
    asnlist = []
    for (filter_key, filter_sci) in sortfiles_dict.items():
        if len(filter_sci) == 0:
            continue

        if bkg_sortfiles_dict is not None:
            filter_backgrounds = bkg_sortfiles_dict.get(filter_key, [])
        else:
            filter_backgrounds = []

        if imp_sortfiles_dict is not None:
            # we assume that alphabetical sorting is sufficient to match the
            # imprints to the right science exposures
            filter_imprints = imp_sortfiles_dict.get(filter_key, [])
            if len(filter_imprints) != len(filter_sci):
                raise RuntimeError(
                    "Number of exposures and imprints is not the same!"
                    f"\n{len(filter_imprints)} imprints and {len(filter_sci)} science for {filter_key}"
                )
        else:
            filter_imprints = []

        sorted_sci = sorted(filter_sci)
        sorted_bkg = sorted(filter_backgrounds)  # not needed, but prettier
        sorted_imp = sorted(filter_imprints)  # needed to match sci to imp!

        if level == 2:
            # one ASN per science exposure
            for i, sci in enumerate(sorted_sci):
                other_exptypes = {}
                if len(filter_backgrounds) > 0:
                    other_exptypes["background"] = sorted_bkg
                if len(filter_imprints) > 0:
                    other_exptypes["imprint"] = [sorted_imp[i]]

                asn_filename = f"{output_dir}/l{level}asn-{filter_key}-{i}.json"
                asn_prodname = f"Level{level}_{filter_key}_{i}"
                asnlist.append(asn_filename)
                writeasn([sci], asn_filename, asn_prodname, other_exptypes)

        elif level == 3:
            # one ASN per filter, with all science exposures
            other_exptypes = {}
            if len(filter_backgrounds) > 0:
                other_exptypes["background"] = filter_backgrounds

            asn_filename = f"{output_dir}/l{level}asn-{filter_key}.json"
            asn_prodname = f"Level{level}_{filter_key}"

            if per_pointing:
                point_groups = group_per_pointing(sorted_sci)
                for point_name, point_files in point_groups.items():
                    point_asn_filename = asn_filename.replace(
                        ".json", f"_{point_name}.json"
                    )
                    point_asn_prodname = asn_prodname + f"_{point_name}"
                    asnlist.append(point_asn_filename)
                    writeasn(
                        point_groups[point_name],
                        point_asn_filename,
                        point_asn_prodname,
                        other_exptypes,
                    )
            else:
                asnlist.append(asn_filename)
                writeasn(sorted_sci, asn_filename, asn_prodname, other_exptypes)

    return asnlist


def group_per_pointing(cal_files):
    """
    Groups per act_id.

    Every filter/pointing combination has a different act_id, common
    between the 4 dither positions.
    """
    pointings = {}
    for f in cal_files:
        header = fits.getheader(f, ext=0)
        act_id = header["ACT_ID"]
        pointings.setdefault(act_id, []).append(f)
    return pointings


if __name__ == "__main__":
    main()
