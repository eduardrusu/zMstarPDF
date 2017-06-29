# The code computes medians and 184h, 84th percentiles from all the unbiased kappa files
# Run as python plotkappacompletestatistics.py WFI2033 5 23 45

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

min_kappa = -0.10
max_kappa = 1
min_kappa_plot = -0.1
max_kappa_plot = 0.30
bin_stat = 2000
halfwidth = (max_kappa - min_kappa) / (bin_stat * 2.0)

root = "/Users/eduardrusu/Dropbox/Davis_work/code/GOODCODE/WFI2033kappa/"

def statistics(kappa_all_,bin_stat_,min_kappa_,max_kappa_):
    a, kappa_values = np.histogram([0], bins = bin_stat_, range=(min_kappa_,max_kappa_)) # create an empty histogram of the correct shape

    sum = np.sum(kappa_all_)
    #meanX = np.sum(kappa_counts * (kappa_values[:-1] + halfwidth)) / sum
    #meanX2 = np.sum(kappa_counts * (kappa_values[:-1] + halfwidth) ** 2) / sum
    #std = np.sqrt(meanX2 - meanX**2)
    
    med = 0
    i = 0
    ok = False
    while (med <= sum/2.0) and (ok == False):
        med = med + kappa_all_[i]
        if med > sum/2.0:
            median = kappa_values[i] + halfwidth
            ok = True
        if med == sum/2.0:
            median = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1
    
    std = 0
    ok = False
    i = 0
    while (std <= sum * 0.16) and (ok == False):
        std = std + kappa_all_[i]
        if std > sum * 0.16:
            std1_ = kappa_values[i] + halfwidth
            ok = True
        if med == sum*0.16:
            std1_ = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1
    
    std = 0
    ok = False
    i = 0
    while (std <= sum*0.84) and (ok == False):
        std = std + kappa_all_[i]
        if std > sum*0.84:
            std1 = kappa_values[i] + halfwidth
            ok = True
        if med == sum*0.84:
            std1 = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1
    
    stddev = (std1 - std1_) / 2
    
    return median,stddev,kappa_values

kappa_counts1 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts1_gal = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_%s_%s_meds_increments2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts1_galgamma = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_%s_%s_meds_increments2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts2 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_z_%s_%s_meds_increments2_2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts2_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_z_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts3 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass_%s_%s_meds_increments2_2_2_6.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts3_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts4 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2_%s_%s_meds_increments2_2_2_8.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts4_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass2_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts5 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3_%s_%s_meds_increments2_2_2_8.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts5_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass3_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts6 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_%s_%s_meds_increments2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts7 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_zoverr_%s_%s_meds_increments2_2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts7_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_zoverr_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts8 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_massoverr_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts8_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_massoverr_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts9 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overr_%s_%s_meds_increments2_2_2_8.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts9_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass2overr_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts10 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overr_%s_%s_meds_increments2_2_2_8.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts10_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass3overrrms_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts11 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2rms_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts11_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass2rms_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts12 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3rms_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts12_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass3rms_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
# NEED TO CORRECT THE NEXT LINE AFTER THE FILE IS FOUND
kappa_counts13 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overrrms_%s_%s_meds_increments2_2_2_8.cat" % (root,lens,inner,mag,str(45)), unpack=True)
kappa_counts13_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass2overrrms_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts14 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overrrms_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts14_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_mass3overrrms_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts15 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_flexion_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts15_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_flexion_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts16 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_tidal_%s_%s_meds_increments2_2_2_4.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts16_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_tidal_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts17 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIS_%s_%s_meds_increments2_2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts17_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_SIS_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts18 = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIShalo_%s_%s_meds_increments2_2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)
kappa_counts18_ = np.loadtxt("%skappahist_%s_%sinnermask_nobeta_gal_oneoverr_SIShalo_%s_%s_meds_increments2_2_2.cat" % (root,lens,inner,mag,rad), unpack=True)

median1,stddev1,kappa_values = statistics(kappa_counts1,bin_stat,min_kappa,max_kappa)
median1_gal,stddev1_gal,kappa_values = statistics(kappa_counts1_gal,bin_stat,min_kappa,max_kappa)
median1_galgamma,stddev1_galgamma,kappa_values = statistics(kappa_counts1_galgamma,bin_stat,min_kappa,max_kappa)
median2,stddev2,kappa_values = statistics(kappa_counts2,bin_stat,min_kappa,max_kappa)
median2_,stddev2_,kappa_values = statistics(kappa_counts2_,bin_stat,min_kappa,max_kappa)
median3,stddev3,kappa_values = statistics(kappa_counts3,bin_stat,min_kappa,max_kappa)
median3_,stddev3_,kappa_values = statistics(kappa_counts3_,bin_stat,min_kappa,max_kappa)
median4,stddev4,kappa_values = statistics(kappa_counts4,bin_stat,min_kappa,max_kappa)
median4_,stddev4_,kappa_values = statistics(kappa_counts4_,bin_stat,min_kappa,max_kappa)
median5,stddev5,kappa_values = statistics(kappa_counts5,bin_stat,min_kappa,max_kappa)
median5_,stddev5_,kappa_values = statistics(kappa_counts5_,bin_stat,min_kappa,max_kappa)
median6,stddev6,kappa_values = statistics(kappa_counts6,bin_stat,min_kappa,max_kappa)
median7,stddev7,kappa_values = statistics(kappa_counts7,bin_stat,min_kappa,max_kappa)
median7_,stddev7_,kappa_values = statistics(kappa_counts7_,bin_stat,min_kappa,max_kappa)
median8,stddev8,kappa_values = statistics(kappa_counts8,bin_stat,min_kappa,max_kappa)
median8_,stddev8_,kappa_values = statistics(kappa_counts8_,bin_stat,min_kappa,max_kappa)
median9,stddev9,kappa_values = statistics(kappa_counts9,bin_stat,min_kappa,max_kappa)
median9_,stddev9_,kappa_values = statistics(kappa_counts9_,bin_stat,min_kappa,max_kappa)
median10,stddev10,kappa_values = statistics(kappa_counts10,bin_stat,min_kappa,max_kappa)
median10_,stddev10_,kappa_values = statistics(kappa_counts10_,bin_stat,min_kappa,max_kappa)
median11,stddev11,kappa_values = statistics(kappa_counts11,bin_stat,min_kappa,max_kappa)
median11_,stddev11_,kappa_values = statistics(kappa_counts11_,bin_stat,min_kappa,max_kappa)
median12,stddev12,kappa_values = statistics(kappa_counts12,bin_stat,min_kappa,max_kappa)
median12_,stddev12_,kappa_values = statistics(kappa_counts12_,bin_stat,min_kappa,max_kappa)
median13,stddev13,kappa_values = statistics(kappa_counts13,bin_stat,min_kappa,max_kappa)
median13_,stddev13_,kappa_values = statistics(kappa_counts13_,bin_stat,min_kappa,max_kappa)
median14,stddev14,kappa_values = statistics(kappa_counts14,bin_stat,min_kappa,max_kappa)
median14_,stddev14_,kappa_values = statistics(kappa_counts14_,bin_stat,min_kappa,max_kappa)
median15,stddev15,kappa_values = statistics(kappa_counts15,bin_stat,min_kappa,max_kappa)
median15_,stddev15_,kappa_values = statistics(kappa_counts15_,bin_stat,min_kappa,max_kappa)
median16,stddev16,kappa_values = statistics(kappa_counts16,bin_stat,min_kappa,max_kappa)
median16_,stddev16_,kappa_values = statistics(kappa_counts16_,bin_stat,min_kappa,max_kappa)
median17,stddev17,kappa_values = statistics(kappa_counts17,bin_stat,min_kappa,max_kappa)
median17_,stddev17_,kappa_values = statistics(kappa_counts17_,bin_stat,min_kappa,max_kappa)
median18,stddev18,kappa_values = statistics(kappa_counts18,bin_stat,min_kappa,max_kappa)
median18_,stddev18_,kappa_values = statistics(kappa_counts18_,bin_stat,min_kappa,max_kappa)

head = "median_1+1/r+ median_1+1/r+gamma+ std_1+1/r+ std_1+1/r+gamma+ "
np.savetxt('%skappastatistics_%s_%s_%s_%s.lst' % (root,lens,inner,mag,rad),np.c_[np.array([median1_gal,median2_,median3_,median4_,median5_,median6,median7_,median8_,median9_,median10_,median11_,median12_,median13_,median14_,median15_,median16_,median17_,median18_]),np.array([median1_galgamma,median2,median3,median4,median5,median1,median7,median8,median9,median10,median11,median12,median13,median14,median15,median16,median17,median18]),np.array([stddev1_gal,stddev2_,stddev3_,stddev4_,stddev5_,stddev6,stddev7_,stddev8_,stddev9_,stddev10_,stddev11_,stddev12_,stddev13_,stddev14_,stddev15_,stddev16_,stddev17_,stddev18_]),np.array([stddev1_galgamma,stddev2,stddev3,stddev4,stddev5,stddev1,stddev7,stddev8,stddev9,stddev10,stddev11,stddev12,stddev13,stddev14,stddev15,stddev16,stddev17,stddev18])],fmt='%s',delimiter='\t',newline='\n',header=head)
