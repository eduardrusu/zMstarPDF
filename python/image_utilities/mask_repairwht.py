# Simple script to change pixel values in custom regions

import numpy as np
from astropy.io import fits

image = fits.open("ch2_4amin_rms.fits")
imagem = fits.open("ch2_4amin_nolens.fits")
image[0].data[imagem[0].data == 0] = 10000000
image.writeto("ch2_4amin_nolens_rms.fits",clobber=True)
