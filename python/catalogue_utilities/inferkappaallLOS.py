# Run as python inferkappaallLOS.py
# Computes kappa for all LOS in the MS, producing a histogram in the same way that the other inferkappa codes do.

import sys
import os
from os import system
import scipy as sp
from scipy.stats import norm
import numpy as np
import time

start_time=time.time()

root = '/mfst01a/rusucs/WFI2033/MSwghtratios/'
output = '/mfst01a/rusucs/WFI2033/MSwghtratios/kappahistallLOS_plane35.cat'

bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

print "Reading..."
for j in range(8):
    for i in range(8):
        if (i!=5) or (j!=5): # because one input field is missing
            kappa_ = np.loadtxt("%snobeta35measuredmedinject_ugriz_WFI2033_GGL_los_8_%s_%s_45_5arcsecinnermsk.cat" % (root,str(j),str(i)), usecols=[1], unpack=True)
            if (i == 0) and (j == 0): kappa = kappa_
            else: kappa = np.append(kappa,kappa_)
    print j

kappa_hist = np.histogram(kappa, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float)/np.sum(kappa_hist)
head = "LOS: "+ str(np.sum(kappa_hist))
np.savetxt(output,kappa_hist,fmt='%s',delimiter='\t',newline='\n',header=head)

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
