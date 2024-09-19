# Where the postprocessing scripts are stored. Needs to be set by user.
ROOT=~/Repositories/pdrs4all

set -e
# point the script to where the pipeline output is
RUNDIR="$1"

JWSTVERSION="$(python -c 'import jwst; print(jwst.__version__)')"
PMAPVERSION="$(python -c 'import crds; print(crds.get_context_name("jwst"))')"

NEWDIR=vX.X_auto_"$PMAPVERSION"_"$JWSTVERSION"
echo mkdir -p "$NEWDIR"
mkdir -p "$NEWDIR"

for d in log crf cubes/default templates
do mkdir -p "$NEWDIR"/"$d"
done

# log files
cp "$RUNDIR"/log/* "$NEWDIR"/log

# crf files
cp "$RUNDIR"/science/stage3/*crf.fits "$NEWDIR"/crf

# default cubes
cp "$RUNDIR"/science/stage3/*s3d.fits "$NEWDIR"/cubes/default

# we are now done copying stuff, so cd into the data release directory
cd "$NEWDIR"

# refined cubes: create crf_asn.json, then run cube_build
mkdir -p cubes/ch1wcs
python3 "$ROOT"/pdrs4all/create_association.py crf crf/*crf.fits
strun cube_build crf_asn.json --output_dir cubes/ch1wcs --cube_pa=250.42338204969806 --scalexy 0.13 --ra_center 83.83535169747073 --dec_center -5.419729392828107  --nspax_x 213 --nspax_y 41 --output_type='band'

# templates
python3 "$ROOT"/pdrs4all/extract_templates.py "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default/*s3d.fits --template_names "HII" "Atomic" "DF3" "DF2" "DF1"
mv templates.ecsv templates/templates_nostitch.ecsv

python3 "$ROOT"/pdrs4all/extract_templates.py "$ROOT"/regions/aper_T_DF_extraction.reg cubes/default/*s3d.fits --template_names "HII" "Atomic" "DF3" "DF2" "DF1" --apply_offsets --reference_segment 2
mv templates.ecsv templates/templates_addstitch.ecsv
