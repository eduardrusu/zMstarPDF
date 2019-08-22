# Create a simple sigma image given the original image and assuming just Poisson noise. This is written for Chih-Fan's iteration code, but can be used in general
# Since the original image is typically sky subtracted, the code uses the sky variance in a selected empty region in the image, and Poisson for everything above 2 sigma of that.
# Also defines a box outside of which the noise is set very large

import sys
import os
import numpy as np
from astropy.io import fits

emptyloc = [30,200] # center of box with no sources in the image, to compute statistics
emptysize = 20 # semilength of box
outsidebox_x = [130,230] # min and max coord outside of wchich noise is considered infinite
outsidebox_y = [80,180]

filein = str(sys.argv[1])
image = fits.open(filein)
data = image[0].data.T # so I work in the natural axis

std = np.std(data[emptyloc[0]-emptysize:emptyloc[0]+emptysize,emptyloc[1]-emptysize:emptyloc[1]+emptysize])
data = np.abs(data) + std ** 2
data = np.sqrt(data)

box = data[outsidebox_x[0]:outsidebox_x[1],outsidebox_y[0]:outsidebox_y[1]]
data[outsidebox_x[0]:outsidebox_x[1],outsidebox_y[0]:outsidebox_y[1]] = -np.abs(box)
data[data > 0] = 100 * np.max(np.abs(box))
data[outsidebox_x[0]:outsidebox_x[1],outsidebox_y[0]:outsidebox_y[1]] = np.abs(box)

image[0].data = data.T.astype("float32") # Hostlens expects this
image.writeto(filein[:-5]+"_sigma.fits",clobber=True)
