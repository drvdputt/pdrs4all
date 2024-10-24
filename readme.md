# PDRs4All

This repository is a collection of data reduction tools for
[PDRs4All](https://pdrs4all.org). It is a specialized version of
[pdr_reduction](https://github.com/drvdputt/jwst-pdr-reduction), with some code
that is specific for the Orion Bar observations. Some of the general
improvements made here should be ported to pdr_reduction.

In addition to running the pipeline with our preferred settins, there are also
tools to create resolution-matched data cubes, merged cubes, stitched cubes,
extracted template spectra.

## Installation

1. Clone this repository
2. Install the python package in your environment by running `pip install -e .`
   in the root directory of this repository. Alternatively, use `poetry
   install`, and then `poetry shell` to create and activate a new environment.
3. Install NSClean manually, see [Paper on
   arxiv](https://arxiv.org/abs/2306.03250), and [download
   page](https://webb.nasa.gov/content/forScientists/publications.html).
   Download and `nsclean_1.9.tar.gz`, then `cd` into `nsclean_1.9/` and run `pip
   install .` in the environment.

## Workflow Part 1: Running the pipeline

PDRs4All consists of NIRCam imaging, MIRI imaging, NIRSpec IFU spectroscopy, and
MIRI IFU spectroscopy. For each of these observations, there is a shell script
provided that performs the pipeline. These shell scripts are generally
applicable (i.e. they should work out of the box for observations other than
PDRs4All), if the following preparation steps are followed.

### Preparing the data

First, sort the `_uncal` according to
1. object
2. instrument
3. exposure type: science, background, and (nirspec only) science imprint,
   background imprint.

Then, copy the appropriate shell script from the `shell_scripts` directory in
this repository to your working directory (one level above the science and
background directories). For example, your directory structure could look like
this.

- `object_1/`
  + `nirspec/`
    - `nirspec.bash`
    - `science/ (*_uncal.fits files here)`
    - `science_imprint/`
    - `background/`
    - `background_imprint/`
  + `mirifu/`
    - `mirifu.bash`
    - `science/`
    - `background/`
  + `nircam/`
    - `nircam.bash`
    - `science/`
    - `background/`
- `object_2/`
  ...


Note: The historical reason for this was to work around some issues with the
default association files. We coded a simplified association generator which
works by simply globbing all the files in the given directories, so it assumes
that the files are sorted correctly, and does not actually apply any association
rules. This manual approach has the advantage that it is fully predictable what
products will processed and which will be used as a background etc.

### Run the shell script

To run the pipeline, one can simply use
```bash
cd object_1/mirifu/
<activate python environment where pdrs4all is installed>
bash mirifu.bash
```

While the shell script is running it will create output directories. Your mirifu
directory from the above example will then look like this

`mirifu/`
 - `mirifu.bash`
 - `science/` <- uncal input files
 - `background/` <- uncal input files
 - `latestpmap/`
   + `science/`
     - `stage1/` <- rate files + stage 2 ASN files are here
     - `stage2/` <- cal files + stage 3 ASN files are here
     - `stage3/` <- cubes are here
   + `background/`
     - `stage1/` <- used in stage 2 for background subtraction
   + `log/`
     - bkg\_1\_cpu1.log <- saved pipeline output, separate for each process (see warning below)
     - bkg\_1.joblog <- overview of the strun jobs that ran. Useful to inspect exit values (if any job failed).
     - sci\_?.log and sci\_?.joblog <- same for all three stages of the science files
     - sci\_?.joblog
     - versions.txt <- contains pmap and jwst version info

Warning about multiprocessing: set the J variable at the top of each shell
script to 1 for an initial run, to avoid conflicts when the references files are
being downloaded. For re-runs (with different options), multiprocessing can be
used safely as long as the CRDS context has not changed. When J is set > 1, the
scripts will use GNU parallel to run multiple pipelines at the same time. I
recommend having at least 8 GB of RAM per process.

## 1D merged spectrum extraction

We provide a script that performs an aperture extraction on the final cubes,
merges the spectral segments, and collects the results in a plain text table. To
use it, the following is needed:
1. A list of data cubes produced by the pipeline (can be the three NIRSpec
   cubes, the 12 MIRI cubes, or both sets)
2. A single region file in the format as produced by DS9. All regions of
   interest should be in one file, in sky coordinates. Currently, only rectangle
   regions are supported

The command is then for example

```
python pdrs4all/extract_templates.py my_regions.reg nirspec/stage3/*s3d.fits miri/stage3/*s3d.fits --template_names Atomic DF
```

where the number of arguments for the optional `--template_names` should equal
the number of regions in the `.reg` file. The output is a file called
`templates.ecsv`, which can be loaded as an astropy table.

See also the oneliner script in `shell_scripts`.

## Current shell script content

### NIRSpec IFU script overview

1. Stage 1 pipeline for `science, science_imprint, background,
   background_imprint`
2. 1/f noise reduction with NSclean (Rauscher, 2023arXiv230603250R), and custom
   masks. TODO: use flicker noise step of detector1 pipeline instead.
3. Stage 2 for `science, background`. Association files that include the NIRSpec
   leakcal imprints (rate files) are created for this.
4. Stage 3 with master background subtraction. Association files that include
   the x1d background produces are created for this. Cube_build is disabled for
   now, as it takes too much memory. See postprocessing step to build the cubes.

### MIRI IFU script overview

1. Stage 1 pipeline for `science, background`
2. Stage 2 pipeline using association files that contain the background rate
   files, to perform image-to-image background subtraction.
3. Stage 3 pipeline using asssociation files, cubes are built per band (12 cubes
   in total).

### Status of imaging scripts

A basic script running the three stages is provided. Intermediate steps for
astrometric correction, wisp removal, 1/f noise reduction, still have to be
added.

## Installation


## Quick start

1. Sort your data (see above)
2. Copy the appropriate bash script from `shell_scripts/` to your working
   directory
3. Edit the copy of the script. Make sure to check the number of processes and
   the CRDS context (pmap number `N`, `CRDS_PATH`, and `CRDS_SERVER_URL`).
4. Activate the environment in which you installed this package (see
   installation instructions above)
5. Run `bash modified_script.bash`

## Credit

PDRs4All data reduction team.
