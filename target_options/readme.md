# What is this

This directory contains `.json` files, which can be given to the
`--custom_options` argument of the `pipeline` script.

These options could be stored here for a variety of purposes, but the most
common reason is that certain targets need specific customizations. 

# orion_cube_pa

The first file added here is a good example: `orion_cube_pa.json` simply sets
the position angle for a `cube_build` for the Orion Bar IFU data. The angle was
chosen to minimize empty space, by aligning the cube coordinate with the
orientation of the 9x1 IFU mosaic.
