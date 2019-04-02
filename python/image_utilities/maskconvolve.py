# Simple script to change pixel values in custom regions extended by convolution

import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft

image = fits.open('/Users/cerusu/OneDrive - Subaru Telescope/PG1115/tmp_obj.fits')
data = image[0].data
data[data == 0] = -1000
kernel = Gaussian2DKernel(10)
conv =  convolve_fft(data, kernel, allow_huge=True)
conv[conv < -100] = 0
conv[conv != 0] = 1
image[0].data = image[0].data * conv
image.writeto('/Users/cerusu/OneDrive - Subaru Telescope/PG1115/tmp_objexpand.fits',clobber=True)
