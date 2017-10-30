# Replace

import numpy as np
from astropy.io import fits

orig = fits.open('patch.fits')
replace = fits.open('LR.20141121.11919_2_bin.fits')
orig_ = orig[4].data
replace_ = replace[0].data
orig_[0:2058,0:369] = replace_
orig[1].data = orig_
orig.writeto('patch.fits',clobber=True)

