# Creates the input catalogue for Fig. 4 in Rusu et al. 2017

import numpy as np
import scipy
from scipy import stats
import sys
import os
from os import system
import time

x,y,i,z,mstar,mhalo = np.loadtxt("HE0435IRACbpz_nobeta_i24.cat",usecols=(0,1,4,7,8,9),unpack=True)

wht_1 = np.ones(len(x))
wht_z = 1.69 * z - z * z
wht_M = 10 ** mstar
wht_M2 = (10 ** mstar) ** 2
wht_M3 = (10 ** mstar) ** 3
wht_1r = 1/(np.sqrt((-600+x)**2 + (-600+y)**2) * 0.2)
wht_zr = wht_z * wht_1r
wht_Mr = wht_M * wht_1r
wht_M2r = wht_M2 * wht_1r
wht_M3r = wht_M3 * wht_1r
wht_flexion = wht_M * (wht_1r ** 3)
wht_tidal = wht_M * (wht_1r ** 2)
wht_SIS = np.sqrt(wht_M) * wht_1r
wht_SIShalo = np.sqrt(10 ** mhalo) * wht_1r
data = np.c_[x,y,i,wht_1,wht_z,wht_M,wht_M2,wht_M3,wht_1r,wht_zr,wht_Mr,wht_M2r,wht_M3r,wht_flexion,wht_tidal,wht_SIS,wht_SIShalo]

np.savetxt('HE0435forweightedcircles.cat',data,fmt='%.2f \t %.2f \t %.2f \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e \t %.4e',header = " x \t y \t i \t wht_1 \t wht_z \t wht_M \t wht_M2 \t wht_M3 \t wht_1r \t wht_zr \t wht_Mr \t wht_M2r \t wht_M3r \t wht_flexion \t wht_tidal \t wht_SIS \t wht_SIShalo")
