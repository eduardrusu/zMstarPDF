# Create a simple blank image

import numpy as np
from astropy.io import fits

image = fits.open('ch1_4amin_nolens.fits')
data = image[0].data
data[data != 0] = 0
image.writeto('flg.fits',clobber=True)
