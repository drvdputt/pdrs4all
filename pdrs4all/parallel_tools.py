# this has problems returning the JWST datamodels (i.e. the output of
# the steps), because they're not pickleable. But the import below
# works! It's a third-party package that uses dill for pickling
from multiprocess import get_context

def run_stage_many(max_cores, obsfiles, stage_class, stage_options_dict=None):
    """
    Run the same stage on many files.

    The call() function will be used.

    Parameters
    ----------

    obsfiles: list of str
        Fits files used as input (strings)

    stage_options_dict: dict
        Dictionary of options passed to the call() function for the
        stage. After unpacking using **, it should look like this
        example from the documentation.

        result = Detector1Pipeline.call('jw00017001001_01101_00001_nrca1_uncal.fits',
                                        steps={'jump': {'threshold': 12.0, 'save_results':True}})

    """
    strun_command_base = f"strun {stage_class.class_alias}"
    options = ""
    for opt_k, opt_v in stage_options_dict.items():
        if opt_k == "steps":
            for step, step_pars in opt_v.items():
                for par, par_v in step_pars.items():
                    options += f" --steps.{step}.{par}={par_v}"
        else:
            options += f" --{opt_k}={opt_v}"

    job_fname = f"strun_{stage_class.class_alias}_jobs.sh"
    with open(job_fname, "w") as f:
        for input_fn in obsfiles:
            strun_command = f"{strun_command_base} {input_fn} {options}"
            f.write(strun_command + "\n")

    print("wrote job file to ", job_fname)
    print("to execute jobs in parallel, use:")
    print(f"parallel -j {max_cores} ::: {job_fname}")


def run_function_many(func, args, max_cores):
    """Run a function that takes a single argument in parallel.

    CAUTION: do not use with jwst pipeline or step objects, since I
    found that there are some problems if the same instance is reused.
    E.g., some of the output files get the wrong name or are
    overwritten. The workaround is run_stage_many().

    """
    if max_cores > 1:

        print("running in parallel on ", args)
        with get_context("spawn").Pool(max_cores) as p:
            _ = p.map(func, args)
            p.close()
            p.join()
    else:
        print("running normal for loop over", args)
        for x in args:
            func(x)
