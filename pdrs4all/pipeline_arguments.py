from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def parse_args():
    """Argument parser for the command line interface of `pypher`"""
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "source_dir", type=str, help="Science folder with /stage0 observations"
    )

    parser.add_argument(
        "-b",
        "--back_dir",
        type=str,
        default=None,
        help="Background folder with /stage0 observations",
    )

    parser.add_argument(
        "-i",
        "--imp_dir",
        type=str,
        default=None,
        help="Imprint folder with /stage0 imprints for observations",
    )

    parser.add_argument(
        "-s",
        "--steps",
        type=str,
        default="123",
        help="Steps of the pipeline (1,2,3,12,23,13,123)",
    )

    parser.add_argument(
        "--mosaic",
        action="store_true",
        help="Merge all pointings in one lvl3 ASN for stage 3.",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default=None,
        help="Output folder. Default is the parent of stage0/.",
    )

    parser.add_argument(
        "--intermediate_dir",
        type=str,
        default=None,
        help="Alternate directory from which intermediate results (stage 1 or 2 products) should be loaded. Defaults to output dir.",
    )

    parser.add_argument(
        "-j",
        "--max_cores",
        type=int,
        default=1,
        help="Number of processes to speed up stage 1 and stage 2. Watch out for memory usage.",
    )

    parser.add_argument(
        "--residual_fringe",
        action="store_true",
        help="Enable residual fringe correction (stage 2 spectroscopy only)",
    )

    parser.add_argument(
        "--custom_options",
        type=str,
        help="""JSON file to define extra pipeline options. Should
        reflect the kwargs given to the call() function of a pipeline.""",
    )

    return parser.parse_args()
