"""Correct NIRSpec WCS and calibration based on NIRCam data.

This implements the worflow by Ryan to improve the NIRSpec data. By
comparing wavelength-integrated NIRSpec cubes to reprojected NIRCam
images, the WCS and flux calibration are improved. At the end, the
improved products are sitched. Here are main steps, with some
interdependency.

1. Create comparable maps

   a. NIRCam images x NIRSpec WCS -> NIRCam maps on same grid (ncrpj)
   b. NIRSpec cubes x NIRCam filter curves -> synthetic photometry (synphot)

2. WCS correction

   a. Measure centroid of proplyd in synphot -> RA, Dec offset
   b. Apply to cubes and synphot -> wcscorr cubes

3. Calibration factor

   a. ncrpj x synphot -> factors
   b. wcscorr cubes x factors -> calcubes (See Peeters 2024 to see which ones to use)

4. Stitch end results
"""

