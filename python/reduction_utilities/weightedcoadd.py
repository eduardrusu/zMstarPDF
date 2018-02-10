# Weighted coaddition of science frames and weighted frames. A final science frame and variance frame is created. Currently only works for variance frames as weighted frames. Uses equations 22 and 23 from the Swarp manual. If regions of the weight files are NaN, they are replaced with infinity.
# Run as: python weightedcoadd sci1.fits sci2.fits ... var1.fits var2.fits ... outsci.fits outvar.fits

import numpy as np
import sys
import math
import numpy as np
from numpy import inf
from astropy.io import fits

nr = (len(sys.argv)-3)/2

for i in range(nr):
    image = fits.open(sys.argv[1+i])
    weight = fits.open(sys.argv[1+i+nr])
    datai = image[0].data
    dataw = weight[0].data
    dataw[np.isnan(dataw)]=np.inf
    if i == 0:
        finalweight = 1.0/dataw # for variance
        finaldatanominator = 1.0/dataw * datai
    else:
        finalweight = finalweight + 1.0/dataw
        finaldatanominator = finaldatanominator + 1.0/dataw * datai

image[0].data = 1.0 * finaldatanominator/finalweight # to preserve the header info
weight[0].data = 1.0/finalweight
image.writeto(sys.argv[-2],clobber=True)
weight.writeto(sys.argv[-1],clobber=True)

