# Replace the wings of a PSF with an analytical model, given the PSF, the model and the mask

import numpy as np
from astropy.io import fits

sky = -256.5 # sky level in the psf
psf = fits.open('i_psf.fits')
mask = fits.open('msknegative.fits')
model = fits.open('i_model.fits')
data_psf = psf[0].data - sky
data_mask = mask[0].data
data_model = model[0].data
out = data_psf
out[data_mask == 1] = data_model[data_mask == 1]
model[0].data = out
model.writeto('i_psfhybrid.fits',clobber=True)

