# first create 3 different association files like: create_association ch1wcs
# latestpmap/science/stage3/*crf.fits. I recommend to adjust the file pattern so that you have
# one association file per nirspec filter. I feel it's a bit more memory efficient that way. Can
# build the cubes one by bone.
mkdir -p ch1wcs
strun cube_build $1 --output_dir ch1wcs --cube_pa=250.42338204969806 --scalexy 0.13 --ra_center 83.83535169747073 --dec_center -5.419729392828107  --nspax_x 213 --nspax_y 41 --output_type='band'
