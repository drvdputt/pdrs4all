import jwst
from jwst import datamodels
from crds.config import get_crds_env_context
from pdr_reduction.pipeline_arguments import parse_args
from pdr_reduction import create_association
from pdr_reduction.parallel_tools import run_stage_many
from pdr_reduction.pipeline_settings import (
    pipeline_class_and_options_dict,
    apply_custom_options,
)
from pathlib import Path

print("JWST pipeline version", jwst.__version__)
print("CRDS context version", get_crds_env_context())


class InstrumentsPipelines:
    def __init__(self, args):
        """
        args: all the things read in by parser.py

        """
        self.stage_numbers = [i for i in (1, 2, 3) if str(i) in args.steps]
        self.obs_dir = args.source_dir
        self.imp_dir = args.imp_dir
        self.back_dir = args.back_dir
        self.custom_options = args.custom_options
        obsfile = sorted(Path(self.obs_dir).glob("*.fits"))[0]
        self.out_dir = (
            str(Path(self.obs_dir) / "..")
            if args.output_dir is None
            else args.output_dir
        )
        self.intermediate_dir = (
            self.out_dir if args.intermediate_dir is None else args.intermediate_dir
        )
        self.per_pointing = not args.mosaic
        self.residual_fringe = args.residual_fringe
        # NRC_IMAGE, MIR_MRS, MIR_IMAGE
        with datamodels.open(obsfile) as im:
            self.instru = im.meta.exposure.type

        self.is_spectroscopy = self.instru == "MIR_MRS" or self.instru == "NRS_IFU"

    def run_pipeline(self, max_cores=1):
        # run requested stages
        for stage in self.stage_numbers:
            # set up output path and make subdir
            output_path = Path(self.out_dir) / f"stage{stage}/"
            output_path.mkdir(parents=True, exist_ok=True)

            # choose input files
            if stage == 1:
                inputs = sorted(
                    [str(p) for p in Path(self.obs_dir).glob("*_uncal.fits")]
                )
            elif stage == 2:
                inputs = create_association.create_asn(
                    self.intermediate_dir,
                    "stage1/*_rate.fits",
                    level=2,
                    backdir=self.back_dir,
                    impdir=self.imp_dir,
                    output_dir=str(Path(self.intermediate_dir) / "stage1/"),
                )
            elif stage == 3:
                inputs = create_association.create_asn(
                    self.intermediate_dir,
                    "stage2/*_cal.fits",
                    level=3,
                    backdir=self.back_dir,
                    spectroscopy=self.is_spectroscopy,
                    per_pointing=self.per_pointing,
                    output_dir=str(Path(self.intermediate_dir) / "stage2/"),
                )

            # get pipeline, options, and run
            pipeline, options = pipeline_class_and_options_dict(
                stage, self.instru, str(output_path)
            )
            # user-defined options, just for this run (not in default pipeline_settings.py)
            if self.custom_options is not None:
                print("Applying custom settings from", self.custom_options)
                apply_custom_options(options, self.custom_options)
                print("Final options are: ", options)

            run_stage_many(max_cores, inputs, pipeline, options)


def main():
    args = parse_args()
    init = InstrumentsPipelines(args)
    init.run_pipeline(args.max_cores)


if __name__ == "__main__":
    main()
