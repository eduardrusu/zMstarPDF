# The code uses the output statistics produced by plotkappacompletestatistics.py/plotkappabiascompletestatistics.py in order to plot bars. Run without arguments. Make sure the uncomment the appropriate ax.set_ylim, ylabel and savefig lines

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np
import sys
import os
from os import system

root = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappa/"
data = np.genfromtxt('%smedstd.dat' % root,dtype=['S1000','f8','f8','f8','f8'])

kappastat = np.array([])
for i in range(np.shape(data)[0]):
    if i == 0:
        kappastat = np.array([data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]])
    else:
        x = np.array([data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]])
        kappastat = np.c_[kappastat,x]


kappastat_45_disjointgaloneoverr = np.c_[   kappastat.T[kappastat[0]=='fiducial_45_gal_22.5_med_increments4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_z_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_22.5_med_increments4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_zoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_massoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2overr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3overr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2rms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3rms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2overrrms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3overrrms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_flexion_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_tidal_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_SIS_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_SIShalo_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

kappastat_120_disjointgaloneoverr = np.c_[  kappastat.T[kappastat[0]=='fiducial_120_gal_22.5_med_increments4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_z_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass2_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass3_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_22.5_med_increments4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_zoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_massoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass2overr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass3overr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass2rms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass3rms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass2overrrms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_mass3overrrms_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_flexion_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_tidal_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_SIS_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_120_SIShalo_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

kappastat_45_disjointgalgammaoneoverr = np.c_[   kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_22.5_med_increments4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_z_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass2_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass3_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_22.5_med_increments4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_zoverr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_massoverr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass2overr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass3overr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass2rms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass3rms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass2overrrms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_mass3overrrms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_flexion_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_tidal_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_SIS_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_45_oneoverr_45_SIShalo_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float)] # SIShalo

kappastat_120_disjointgalgammaoneoverr = np.c_[  kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_22.5_med_increments4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_z_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass2_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass3_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_22.5_med_increments4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_zoverr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_massoverr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass2overr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass3overr_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass2rms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass3rms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass2overrrms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_mass3overrrms_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_flexion_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_tidal_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_SIS_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_120_SIShalo_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float)] # SIShalo

kappastat_45_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr = np.c_[   kappastat.T[kappastat[0]=='fiducial_45_gal_120_gal_120_oneoverr_22.5_med_increments1_1_1_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_z_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass2_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass3_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_zoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_massoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass2overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass3overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass2rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass3rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass2overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_mass3overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_flexion_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_tidal_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_oneoverr_120_SIShalo_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

kappastat_120_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr = np.c_[  kappastat.T[kappastat[0]=='fiducial_120_gal_45_gal_45_zoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_z_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass2_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass3_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_oneoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_massoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass2overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass3overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass2rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass3rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass2overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_mass3overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_flexion_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_tidal_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_zoverr_120_SIShalo_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

kappastat_45_conjointgaloneoverrgalgamma = np.c_[   kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_z_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass2_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass3_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_oneoverr_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_zoverr_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_massoverr_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass2overr_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass3overr_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass2rms_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass3rms_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass2overrrms_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_mass3overrrms_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_flexion_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_tidal_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_120_SIShalo_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float)] # SIShalo

kappastat_120_conjointgaloneoverrgalgamma = np.c_[  kappastat.T[kappastat[0]=='fiducial_45_gal_120_gal_120_gamma_22.5_med_increments4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_z_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_zoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_massoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2overr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3overr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2rms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3rms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass2overrrms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_mass3overrrms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_flexion_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_tidal_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_SIS_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_oneoverr_45_SIShalo_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk_shearwithoutprior.cat'][0][1:].astype(float)] # SIShalo

N = 18
ind = 2.5 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

fig, axs = plt.subplots(nrows=4, ncols=2, sharex=True, figsize=(10, 10))
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0,wspace=0)


col1 = (kappastat_45_disjointgaloneoverr[0])
rects1 = axs[0,0].bar(ind + width, col1, width, color='gray')
col2 = (kappastat_120_disjointgaloneoverr[0])
rects2 = axs[0,0].bar(ind + 2*width, col2, width, color='k')

#ax.set_ylim([0.00,0.05])
axs[0,0].set_ylim([-0.01,0.08])
axs[0,0].set_ylabel('$\kappa^{med}_\mathrm{ext}$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[0,0].set_xticks(ind + 2*width)
axs[0,0].set_xticklabels((), fontsize=10, rotation='vertical')
axs[0,0].set_yticks(np.arange(-0.01, 0.09, 0.01))

col3 = (kappastat_45_disjointgaloneoverr[1])
rects3 = axs[1,0].bar(ind + width, col3, width, color='gray')
col4 = (kappastat_120_disjointgaloneoverr[1])
rects4 = axs[1,0].bar(ind + 2*width, col4, width, color='k')

axs[1,0].set_ylim([0,0.08])
axs[1,0].set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[1,0].set_xticks(ind + width)
axs[1,0].set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
axs[1,0].legend((rects1[0], rects2[0]), ('$P(\kappa_\mathrm{ext}|\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{1/r},\zeta^{120\'\'}_{q})$', '$P(\kappa_\mathrm{ext}|\zeta^{120\'\'}_{1},\zeta^{120\'\'}_{1/r},\zeta^{120\'\'}_{q})$'), bbox_to_anchor=(1.0, 1.0), fontsize=12)
axs[1,0].set_yticks(np.arange(0.00, 0.08, 0.01))


col5 = (kappastat_45_disjointgalgammaoneoverr[0])
rects5 = axs[0,1].bar(ind + width, col5, width, color='gray')
col6 = (kappastat_120_disjointgalgammaoneoverr[0])
rects6 = axs[0,1].bar(ind + 2*width, col6, width, color='k')

#ax.set_ylim([0.00,0.05])
axs[0,1].set_ylim([-0.01,0.08])
#axs[0,1].set_ylabel('$\kappa^{med}_\mathrm{ext}$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[0,1].set_xticks(ind + 2*width)
axs[0,1].set_xticklabels((), fontsize=10, rotation='vertical')
axs[0,1].set_yticklabels((), fontsize=10, rotation='vertical')


col7 = (kappastat_45_disjointgalgammaoneoverr[1])
rects7 = axs[1,1].bar(ind + width, col7, width, color='gray')
col8 = (kappastat_120_disjointgalgammaoneoverr[1])
rects8 = axs[1,1].bar(ind + 2*width, col8, width, color='k')

axs[1,1].set_ylim([0,0.08])
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[1,1].set_xticks(ind + width)
axs[1,1].set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
axs[1,1].legend((rects1[0], rects2[0]), ('$P(\kappa_\mathrm{ext}|\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{1/r},\gamma,\zeta^{45\'\'}_{q})$', '$P(\kappa_\mathrm{ext}|\zeta^{120\'\'}_{1},\zeta^{120\'\'}_{1/r},\gamma,\zeta^{120\'\'}_{q})$'), bbox_to_anchor=(1, 1), loc=1, fontsize=12)
axs[1,1].set_yticklabels((), fontsize=10, rotation='vertical')

col9 = (kappastat_45_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr[0])
rects9 = axs[2,0].bar(ind + width, col9, width, color='gray')
col10 = (kappastat_120_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr[0])
rects10 = axs[2,0].bar(ind + 2*width, col10, width, color='k')


#ax.set_ylim([0.00,0.05])
axs[2,0].set_ylim([-0.01,0.08])
axs[2,0].set_ylabel('$\kappa^{med}_\mathrm{ext}$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[2,0].set_xticks(ind + 2*width)
axs[2,0].set_xticklabels((), fontsize=10, rotation='vertical')
axs[2,0].set_yticks(np.arange(-0.01, 0.08, 0.01))

col11 = (kappastat_45_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr[1])
rects11 = axs[3,0].bar(ind + width, col11, width, color='gray')
col12 = (kappastat_120_conjointgaloneoverrgaloneoverr_and_galzoverrgalzoverr[1])
rects12 = axs[3,0].bar(ind + 2*width, col12, width, color='k')

axs[3,0].set_ylim([0,0.08])
axs[3,0].set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[3,0].set_xticks(ind + width)
axs[3,0].set_xticklabels(('$1-(1/r,z/r)$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
axs[3,0].legend((rects1[0], rects2[0]), ('$P(\kappa_\mathrm{ext}|\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{1/r},\zeta^{120\'\'}_{1},\zeta^{120\'\'}_{1/r},\zeta^{120\'\'}_{q})$', '$P(\kappa_\mathrm{ext}|\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{z/r}$,$\zeta^{120\'\'}_{1},\zeta^{120\'\'}_{z/r},\zeta^{120\'\'}_{q})$'), bbox_to_anchor=(1.0, 1.0), loc=1, fontsize=12)
axs[3,0].set_yticks(np.arange(0.00, 0.08, 0.01))

col13 = (kappastat_45_conjointgaloneoverrgalgamma[0])
rects13 = axs[2,1].bar(ind + width, col13, width, color='gray')
col14 = (kappastat_120_conjointgaloneoverrgalgamma[0])
rects14 = axs[2,1].bar(ind + 2*width, col14, width, color='k')

#ax.set_ylim([0.00,0.05])
axs[2,1].set_ylim([-0.01,0.08])
#axs[0,1].set_ylabel('$\kappa^{med}_\mathrm{ext}$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[2,1].set_xticks(ind + 2*width)
axs[2,1].set_xticklabels((), fontsize=10, rotation='vertical')
axs[2,1].set_yticklabels((), fontsize=10, rotation='vertical')


col15 = (kappastat_45_conjointgaloneoverrgalgamma[1])
rects15 = axs[3,1].bar(ind + width, col15, width, color='gray')
col16 = (kappastat_120_conjointgaloneoverrgalgamma[1])
rects16 = axs[3,1].bar(ind + 2*width, col16, width, color='k')

axs[3,1].set_ylim([0,0.08])
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
axs[3,1].set_xticks(ind + width)
axs[3,1].set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
axs[3,1].legend((rects1[0], rects2[0]), ('$P(\kappa_\mathrm{ext}|\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{1/r},\zeta^{120\'\'}_{1},\gamma,\zeta^{120\'\'}_{q})$', '$P(\kappa_\mathrm{ext}|\zeta^{120\'\'}_{1},\gamma,\zeta^{45\'\'}_{1},\zeta^{45\'\'}_{1/r},\zeta^{45\'\'}_{q})$'), bbox_to_anchor=(1, 1), loc=1, fontsize=12)
axs[3,1].set_yticklabels((), fontsize=10, rotation='vertical')


#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
#plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
plt.subplots_adjust(bottom=0.15,top=0.95)
plt.savefig('%skappashistbar.png' % root, dpi=250)

plt.clf()
