import jwst
from jwst import datamodels
from crds.config import get_crds_env_context
from pdrs4all.pipeline_arguments import parse_args
from pdrs4all import create_association
from pdrs4all.pipeline_settings import (
    pipeline_class_and_options_dict,
    apply_custom_options,
)
from pathlib import Path

print("JWST pipeline version", jwst.__version__)
print("CRDS context version", get_crds_env_context())


def create_strun_jobs(args):
    """Parse default options, custom options, and set up strun
    commands.

    Output
    ------
    Writes job file containing one strun command per input file.

    (Per association file in the case of stages 2 and 3)

    """

    stage_numbers = [i for i in (1, 2, 3) if str(i) in args.steps]
    args.custom_options = args.custom_options
    obsfile = sorted(Path(args.source_dir).glob("*.fits"))[0]
    out_dir = (
        str(Path(args.source_dir) / "..")
        if args.output_dir is None
        else args.output_dir
    )
    intermediate_dir = (
        out_dir if args.intermediate_dir is None else args.intermediate_dir
    )
    per_pointing = not args.mosaic

    # NRC_IMAGE, MIR_MRS, MIR_IMAGE
    with datamodels.open(obsfile) as im:
        instru = im.meta.exposure.type

    is_spectroscopy = instru == "MIR_MRS" or instru == "NRS_IFU"

    # run requested stages
    for stage in stage_numbers:
        # set up output path and make subdir
        output_path = Path(out_dir) / f"stage{stage}/"
        output_path.mkdir(parents=True, exist_ok=True)

        # choose input files
        if stage == 1:
            inputs = sorted(
                [str(p) for p in Path(args.source_dir).glob("*_uncal.fits")]
            )
        elif stage == 2:
            inputs = create_association.create_asn(
                intermediate_dir,
                "stage1/*_rate.fits",
                level=2,
                backdir=args.back_dir,
                impdir=args.imp_dir,
                output_dir=str(Path(intermediate_dir) / "stage1/"),
            )
        elif stage == 3:
            inputs = create_association.create_asn(
                intermediate_dir,
                "stage2/*_cal.fits",
                level=3,
                backdir=args.back_dir,
                spectroscopy=is_spectroscopy,
                per_pointing=per_pointing,
                output_dir=str(Path(intermediate_dir) / "stage2/"),
            )

        # get pipeline type and options
        pipeline, options = pipeline_class_and_options_dict(
            stage, instru, str(output_path)
        )
        # user-defined options, just for this run (not in default pipeline_settings.py)
        if args.custom_options is not None:
            print("Applying custom settings from", args.custom_options)
            apply_custom_options(options, args.custom_options)
            print("Final options are: ", options)

        # convert options dict to command line options string
        strun_command_base = f"strun {pipeline.class_alias}"
        strun_options = ""
        for opt_k, opt_v in options.items():
            if opt_k == "steps":
                for step, step_pars in opt_v.items():
                    for par, par_v in step_pars.items():
                        strun_options += f" --steps.{step}.{par}={par_v}"
            else:
                strun_options += f" --{opt_k}={opt_v}"

        job_fname = f"strun_{pipeline.class_alias}_jobs.sh"
        with open(job_fname, "w") as f:
            for input_fn in inputs:
                strun_command = f"{strun_command_base} {input_fn} {strun_options}"
                f.write(strun_command + "\n")

        print("wrote job file to ", job_fname)


def main():
    args = parse_args()
    create_strun_jobs(args)


if __name__ == "__main__":
    main()
