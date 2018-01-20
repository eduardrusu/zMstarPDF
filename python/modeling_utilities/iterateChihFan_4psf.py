# I first copy r_nolens_out_psf.fits into rnolenspsfcorrect.fits. Then I edit config_4psfrnolens to contain the correct input names, and this code to read the correct lines. After that I modify r_nolens_out_file.input to use rnolenspsfcorrect.fits as psf, read sigmar.fits, disalow analytical PSF parameters (except for sky) and run the code as follows:
# run as python iterateChihFan_4psf.py r_nolens_out_file.input config_4psfrnolens rnolenspsfcorrect.fits

import sys
import os 
import numpy as np
from astropy.io import fits

iterations = 3
pixscale = 0.256

for i in range(iterations):
    with open(sys.argv[1], 'r') as fileinput:
        hostlens = fileinput.readlines()

    x1 = float(hostlens[43-1].split()[1]) # brightest
    y1 = float(hostlens[44-1].split()[1])
    f1 = float(hostlens[50-1].split()[1])
    x2 = float(hostlens[54-1].split()[1]) # 2nd brightest
    y2 = float(hostlens[55-1].split()[1])
    f2 = float(hostlens[61-1].split()[1])
    x3 = float(hostlens[65-1].split()[1]) # 3rd brightest
    y3 = float(hostlens[66-1].split()[1])
    f3 = float(hostlens[72-1].split()[1])
    x4 = float(hostlens[76-1].split()[1]) # 3rd brightest
    y4 = float(hostlens[77-1].split()[1])
    f4 = float(hostlens[83-1].split()[1])

    image = fits.open(str(sys.argv[1])[:-10]+"psf.fits")
    data = image[0].data * f2
    imagex = image
    imagex[0].data = data
    imagex.writeto(str(sys.argv[1])[:-10]+"psf.fits",clobber=True)

    with open(sys.argv[2], 'r') as fileconfig:
        config = fileconfig.readlines()

    config[6] = "x1_in_arcsec " + "%.6f" % (x1 * pixscale) + "\n"
    config[7] = "y1_in_arcsec " + "%.6f" % (y1 * pixscale) + "\n"
    config[8] = "intensityof1(weak) " + "%.6e" % f1 + "\n"
    config[9] = "x2_in_arcsec " + "%.6f" % (x2 * pixscale) + "\n"
    config[10] = "y2_in_arcsec " + "%.6f" % (y2 * pixscale) + "\n"
    config[11] = "intensityof2(strong) " + "%.6e" % f2 + "\n"
    config[12] = "x3_in_arcsec " + "%.6f" % (x3 * pixscale) + "\n"
    config[13] = "y3_in_arcsec " + "%.6f" % (y3 * pixscale) + "\n"
    config[14] = "intensityof3(weak) " + "%.6e" % f3 + "\n"
    config[15] = "x4_in_arcsec " + "%.6f" % (x4 * pixscale) + "\n"
    config[16] = "y4_in_arcsec " + "%.6f" % (y4 * pixscale) + "\n"
    config[17] = "intensityof4(weak) " + "%.6e" % f4 + "\n"

    with open(sys.argv[2], 'w') as fileconfig:
        fileconfig.writelines(config)

    fileinput.close()
    fileconfig.close()

    os.system("python PSF_correction_4psf.py %s %s" %(str(sys.argv[2]),str(sys.argv[3])))
    os.system("hostlens %s" % str(sys.argv[1]))

for i in range(3):
    os.system("hostlens %s" % str(sys.argv[1])) # because hostlens needs a few executions in order to fully converge
