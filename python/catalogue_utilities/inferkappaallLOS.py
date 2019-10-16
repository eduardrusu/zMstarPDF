# Run as python inferkappaallLOS.py
# Computes kappa for all LOS in the MS, producing a histogram in the same way that the other inferkappa codes do.

import sys
import os
from os import system
import scipy as sp
from scipy.stats import norm
import numpy as np
import time
import fitsio

start_time=time.time()

root = '/lfs08/rusucs/0408/MSwghtratios/'
output = '/lfs08/rusucs/0408/MSkapparesults/kappahistallLOS_plane30.cat'

bin_stat = 2000
min_kappa = -0.20
max_kappa = 1

print "Reading..."
for j in range(8):
    for i in range(8):
        file = "%snobeta30measuredmedinject_griz_lens_0408_GGL_los_8_%s_%s_22.5_45_5arcsecinner_gap_-1.0_-1.0.fits" % (root,str(j),str(i))
        f = fitsio.FITS(file)
        print f # I need to print it, or f.hdu_list will not read
        ext = len(f.hdu_list)
        #if (i!=5) or (j!=5): # because one input field is missing
            #kappa_ = np.loadtxt("%snobeta35measuredmedinject_ugriz_WFI2033_GGL_los_8_%s_%s_45_5arcsecinnermsk.cat" % (root,str(j),str(i)), usecols=[1], unpack=True)
        if (i == 0) and (j == 0): kappa = fitsio.read(file, columns=[1]).astype(float)
        else: kappa = np.append(kappa,fitsio.read(file, columns=[1]).astype(float))
    print j

kappa_hist = np.histogram(kappa, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float)
kappa_hist = kappa_hist/np.sum(kappa_hist)
head = "LOS: "+ str(len(kappa))
np.savetxt(output,kappa_hist,fmt='%s',delimiter='\t',newline='\n',header=head)

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
