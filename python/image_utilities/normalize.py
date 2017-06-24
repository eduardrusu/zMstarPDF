# Simple code to normalize an image

import numpy as np
from astropy.io import fits
import sys

image = fits.open(str(sys.argv[1]))
data = image[0].data
sum = np.sum(image[0].data)
data = image[0].data/sum
imagex = image
imagex[0].data = data
imagex.writeto(str(sys.argv[2]),clobber=True)

