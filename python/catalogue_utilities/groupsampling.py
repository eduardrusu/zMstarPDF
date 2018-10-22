# Run this code in order to sample possible galaxies part of an incomplete spectroscopic group

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
#import scipy
#from scipy import stats
#import sys
#import os
#from os import system
#import time
zgroup = 0.66
file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
filephotozpool = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/grouppool.cat"
x = np.loadtxt(file,usecols=[2,3,4,8,28,29,30,40,97])
ra = 0
dec = 1
i = 2
id = 3
z = 4
zinf = 5
zsup = 6
spec = 7
cls = 8
center_lens = SkyCoord('20:33:42.080 -47:23:43.00', frame='fk5', unit=(u.hourangle, u.deg))
coord = SkyCoord(ra=x[:,ra]*u.degree, dec=x[:,dec]*u.degree, frame='fk5')
sep = coord.separation(center_lens).arcsec
print x[:,id][(x[:,spec] <= zgroup + 0.01) & (x[:,spec] >= zgroup - 0.01)]
