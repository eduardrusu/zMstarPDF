# Simple script to change pixel values in custom regions

import numpy as np
from astropy.io import fits

image = fits.open("ch2_4amin_nolens.fits")
imagem = fits.open("msk.fits")
image[0].data[(imagem[0].data == 1) & (image[0].data < 0)] = 0
image.writeto("ch2_4amin_nolens.fits",clobber=True)
