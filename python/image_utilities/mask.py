# Simple script to change pixel values in custom regions

import numpy as np
from astropy.io import fits

image = fits.open('gaussswarpWFI2033_F814W_sci.fits')
weight = fits.open('gaussswarpWFI2033_F814W_sci_wht.fits')
datai = image[0].data
dataw = weight[0].data
datai[dataw == 0] = 0
image.writeto('gaussswarpWFI2033_F814W_sci.fits',clobber=True)

