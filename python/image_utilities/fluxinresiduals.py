# Given a mask marking the object footprint and a residual image, adds up the flux in the residuals
# Run as python /Users/cerusu/GITHUB/zMstarPDF/python/image_utilities/fluxinresiduals.py /Users/cerusu/Dropbox/Glikman/analysis/msk_visit1A.fits /Users/cerusu/Dropbox/Glikman/analysis/network/125_1sers4/out_lens1pertpointextendedoriginalpsfmskarc_iter30pymc_1imagonly_subtract.fits

import numpy as np
from astropy.io import fits
import sys

mask = fits.open(str(sys.argv[1]))
mask = mask[0].data
resid = fits.open(str(sys.argv[2]))
resid = resid[0].data
print np.sum(resid[mask==1])
