# Interpolates a .fits image. Used here in order to sample an IRAC PSF to the pixel scale of the Tphot registration catalogue.

import numpy as np
import scipy
from astropy.io import fits
from scipy import ndimage
origimg = fits.open('ker_i_ch1.fits')
origdata = origimg[0].data
origarray = np.array(origdata)
zoomarray = scipy.ndimage.interpolation.zoom(origarray,0.9)
zoomdata = fits.PrimaryHDU(zoomarray)
zoomdata.writeto('ker_i_ch1_zoom0.9_.fits')
origimg.close()

# after this, do in IRAF:
# imcopy ker_i_ch1_zoom0.9_.fits[6:104,6:104] ker_i_ch1_zoom0.9.fits
# imdel ker_i_ch1_zoom0.9_.fits
