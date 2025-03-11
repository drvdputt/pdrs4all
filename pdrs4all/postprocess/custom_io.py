from astropy.nddata import CCDData
from astropy import units as u
from jwst import datamodels

def write_i2d(fits_fn, array, wcs, unit=u.MJy / u.sr):
    """Write numpy array and to FITS file with WCS

    Parameters
    ----------

    fits_fn : path
        FITS file name for output

    array : 2D array-like
        Flux data

    wcs : WCS
        WCS to be added
    """
    ccd = CCDData(array.T, wcs=wcs, unit=unit)
    ccd.write(fits_fn, overwrite=True)

def write_s3d_with_new_crval(output_fn, original_fn, crval):
    """Load cube, edit WCS, and write again.

    Parameters
    ----------

    output_fn : path
        Output file, should end in s3d.fits

    original_fn : path
        Original cube of which the CRVAL will be edited

    crval : pair of float
        Usually new_wcs.wcs.crval
    
    """
    cube_dm = datamodels.open(original_fn)
    cube_dm.meta.wcsinfo.crval1 = crval[0]
    cube_dm.meta.wcsinfo.crval2 = crval[1]
    cube_dm.write(output_fn)
