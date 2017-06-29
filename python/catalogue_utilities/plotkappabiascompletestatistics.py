# The code computes medians and 184h, 84th percentiles from all the unbiased kappa files
# Run as python plotkappabiascompletecompute.py WFI2033 5 23 45

import scipy as sp
from scipy.stats import norm
import numpy as np
import sys
import os
from os import system

lens = str(sys.argv[1])
inner = str(sys.argv[2])
mag = str(sys.argv[3])
rad = str(sys.argv[4])

root = "/Users/eduardrusu/Dropbox/Davis_work/code/GOODCODE/WFI2033kappa/"

def compute(kappa_counts_):
    if np.median(kappa_counts_[4]) > 2:
        median = np.median(kappa_counts_[1] - kappa_counts_[0])
        std = (np.percentile(kappa_counts_[1] - kappa_counts_[0],84) - np.percentile(kappa_counts_[1] - kappa_counts_[0],16))/2.0
    else: median = float('nan'); std = float('nan');
    print np.median((kappa_counts_[3] - kappa_counts_[2])/2.0),np.median(kappa_counts_[4])
    return median, std

kappa_counts1 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts1_gal = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts1_galgamma = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts2 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_z_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts2_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_z_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts3 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts3_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts4 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts4_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass2_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts5 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts5_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass3_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts6 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts7 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_zoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts7_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_zoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts8 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_massoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts8_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_massoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts9 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts9_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass2overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts10 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts10_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass3overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts11 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts11_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass2rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts12 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts12_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass3rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts13 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts13_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass2overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts14 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts14_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_mass3overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts15 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_flexion_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts15_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_flexion_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts16 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_tidal_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts16_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_tidal_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts17 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIS_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts17_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_SIS_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts18 = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIShalo_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts18_ = np.loadtxt("%skappasimbias_%s_%sinnermask_nobeta_gal_oneoverr_SIShalo_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)

median1,stddev1 = compute(kappa_counts1)
median1_gal,stddev1_gal = compute(kappa_counts1_gal)
median1_galgamma,stddev1_galgamma = compute(kappa_counts1_galgamma)
median2,stddev2 = compute(kappa_counts2)
median2_,stddev2_ = compute(kappa_counts2_)
median3,stddev3 = compute(kappa_counts3)
median3_,stddev3_ = compute(kappa_counts3_)
median4,stddev4 = compute(kappa_counts4)
median4_,stddev4_ = compute(kappa_counts4_)
median5,stddev5 = compute(kappa_counts5)
median5_,stddev5_ = compute(kappa_counts5_)
median6,stddev6 = compute(kappa_counts6)
median7,stddev7 = compute(kappa_counts7)
median7_,stddev7_ = compute(kappa_counts7_)
median8,stddev8 = compute(kappa_counts8)
median8_,stddev8_ = compute(kappa_counts8_)
median9,stddev9 = compute(kappa_counts9)
median9_,stddev9_ = compute(kappa_counts9_)
median10,stddev10 = compute(kappa_counts10)
median10_,stddev10_ = compute(kappa_counts10_)
median11,stddev11 = compute(kappa_counts11)
median11_,stddev11_ = compute(kappa_counts11_)
median12,stddev12 = compute(kappa_counts12)
median12_,stddev12_ = compute(kappa_counts12_)
median13,stddev13 = compute(kappa_counts13)
median13_,stddev13_ = compute(kappa_counts13_)
median14,stddev14 = compute(kappa_counts14)
median14_,stddev14_ = compute(kappa_counts14_)
median15,stddev15 = compute(kappa_counts15)
median15_,stddev15_ = compute(kappa_counts15_)
median16,stddev16 = compute(kappa_counts16)
median16_,stddev16_ = compute(kappa_counts16_)
median17,stddev17 = compute(kappa_counts17)
median17_,stddev17_ = compute(kappa_counts17_)
median18,stddev18 = compute(kappa_counts18)
median18_,stddev18_ = compute(kappa_counts18_)

head = "median_1+1/r+ median_1+1/r+gamma+ std_1+1/r+ std_1+1/r+gamma+ "
np.savetxt('%skappacomputebias_%s_%s_%s_%s.lst' % (root,lens,inner,mag,rad),np.c_[np.array([median1_gal,median2_,median3_,median4_,median5_,median6,median7_,median8_,median9_,median10_,median11_,median12_,median13_,median14_,median15_,median16_,median17_,median18_]),np.array([median1_galgamma,median2,median3,median4,median5,median1,median7,median8,median9,median10,median11,median12,median13,median14,median15,median16,median17,median18]),np.array([stddev1_gal,stddev2_,stddev3_,stddev4_,stddev5_,stddev6,stddev7_,stddev8_,stddev9_,stddev10_,stddev11_,stddev12_,stddev13_,stddev14_,stddev15_,stddev16_,stddev17_,stddev18_]),np.array([stddev1_galgamma,stddev2,stddev3,stddev4,stddev5,stddev1,stddev7,stddev8,stddev9,stddev10,stddev11,stddev12,stddev13,stddev14,stddev15,stddev16,stddev17,stddev18])],fmt='%s',delimiter='\t',newline='\n',header=head)
