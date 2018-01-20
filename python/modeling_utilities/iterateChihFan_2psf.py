# I first edit gconfig to contain the correct input names. After that I modify g_nolens_out_file.input to use gpsfcorrect.fits as psf, disalow analytical PSF parameters (except for sky) and run the code as follows:
# run as python iterateChihFan.py i_lens_out_file.input config psfcorrect.fits

import sys
import os
import numpy as np
from astropy.io import fits

iterations = 3
pixscale = 0.256

for i in range(iterations):
    with open(sys.argv[1], 'r') as fileinput:
        hostlens = fileinput.readlines()

    x1 = float(hostlens[54-1].split()[1]) # faint
    y1 = float(hostlens[55-1].split()[1])
    f1 = float(hostlens[61-1].split()[1])
    x2 = float(hostlens[43-1].split()[1]) # bright
    y2 = float(hostlens[44-1].split()[1])
    f2 = float(hostlens[50-1].split()[1])

    image = fits.open(str(sys.argv[1])[:-10]+"psf.fits")
    data = image[0].data * f2
    imagex = image
    imagex[0].data = data
    imagex.writeto(str(sys.argv[3]),clobber=True)

    with open(sys.argv[2], 'r') as fileconfig:
        config = fileconfig.readlines()

    config[6] = "x1_in_arcsec " + "%.6f" % (x1 * pixscale) + "\n"
    config[7] = "y1_in_arcsec " + "%.6f" % (y1 * pixscale) + "\n"
    config[8] = "x2_in_arcsec " + "%.6f" % (x2 * pixscale) + "\n"
    config[9] = "y2_in_arcsec " + "%.6f" % (y2 * pixscale) + "\n"
    config[10] = "intensityof1(weak) " + "%.6e" % f1 + "\n"
    config[11] = "intensityof2(strong) " + "%.6e" % f2 + "\n"

    with open(sys.argv[2], 'w') as fileconfig:
        fileconfig.writelines(config)

    fileinput.close()
    fileconfig.close()

    os.system("python PSF_correction_2psf.py %s %s" %(str(sys.argv[2]),str(sys.argv[3])))
    os.system("hostlens %s" % str(sys.argv[1]))

for i in range(1):
    os.system("hostlens %s" % str(sys.argv[1])) # because hostlens needs a few executions in order to fully converge
