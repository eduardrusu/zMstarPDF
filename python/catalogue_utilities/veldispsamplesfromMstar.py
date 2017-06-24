# Calculates the pdf for the velocity dispersion of each galaxy. Check Jabran Zahid et al 2016 for the applicability of the formulas. The input is a file that contains 100 (z,Mstar) pairs for each galaxy, on each line, and is produced by mstarsampling_nobeta_forveldisperr.py

import numpy as np
import scipy
from scipy import stats
import sys
import os
from os import system
from astropy import units as u
from astropy.coordinates import SkyCoord

####################
# Jabran Zahid et al 2017 parameters:
sb = 10**(2.07)
Mb = 10**(10.26)
a1 = 0.403
a2 = 0.293
#veldisperr = lambda x: 0.001*(0.1 * x + 14) # [km/s]; eyeballed from Fig. 9(a)
veldisperr = lambda x: 0.1 * x + 14 # [km/s]; eyeballed from Fig. 9(a)

####################
lens = SkyCoord('20:33:42.16 -47:23:44.20', frame='fk5', unit=(u.hourangle, u.deg)) # center of the lensing galaxy

# read from file
ra = 2
dec = 3
imag = 4
classify = 7
z0 = 8
mstar0 = 9

file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_forveldisperr.cat"
data = np.loadtxt(file)
data = data[(data[:,imag]<=23) & (data[:,classify]>=0)] # only i <= 23 galaxies
veldisp = np.zeros([np.shape(data)[0],100])

for i in range(np.shape(data)[0]):
    mstar = data[i][mstar0:208:2]
    #mstar = mstar - np.log10(0.55) # not converting to Salpeter IMF, because me and Zahid et al 2016 calculated masses with Chabrier
    coord = SkyCoord(ra=data[i][ra]*u.degree, dec=data[i][dec]*u.degree, frame='fk5')
    sep = coord.separation(lens).arcsec
    z = data[i][z0:207:2] # in the current code, it is not used for anything
    veldisp[i][10**mstar < Mb] = np.abs(np.random.normal(sb * ((10**mstar[10**mstar < Mb] / Mb)**a1),veldisperr((10**mstar[10**mstar < Mb] / Mb)**a1)))
    veldisp[i][10**mstar >= Mb] = np.abs(np.random.normal(sb * ((10**mstar[10**mstar >= Mb] / Mb)**a2),veldisperr(10**mstar[10**mstar >= Mb] / Mb)**a2))

np.savetxt(file[:-4] + "_veldisppdf.cat",veldisp,fmt='%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f')
