# Simple script that combines Sextractor OBJECTS and BACKGROUND_RMS images into a sigma map. NGAIN is NCOMBINExGAIN
# run as python sexnoisemap.py obj.fits rms.fits out.fits NGAIN

import numpy as np
import sys
from astropy.io import fits

ngain = float(str(sys.argv[4]))
im = fits.open(str(sys.argv[1]))
im[0].data = im[0].data / ngain

imr = fits.open(str(sys.argv[2]))
imr[0].data = imr[0].data ** 2

out = np.sqrt(im[0].data + imr[0].data)
out[out < 0] = np.median(out)

im[0].data = out
im.writeto(str(sys.argv[3]),clobber=True)
