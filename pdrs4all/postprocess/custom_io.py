from astropy.nddata import CCDData
from astropy import units as u

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
