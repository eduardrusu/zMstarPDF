# Run as python inferkappa_unbiasedwithshear.py WFI2033 5 45 23 meds gal gamma oneoverr mass
# Description of arguments: inferkappa_unbiasedwithshear.py lens radius maglim innermask sum/meds gal list_of_weight_constraints
# weight1 should always be "gal", in order to use the galaxy counts when correcting the bias due to different LOS
# the code is written such that, if shear is used as overdensity, it should be the second weight used (unless only one weight is used);

import sys
import os
from os import system
import scipy as sp
from scipy.stats import norm
import numpy as np
import time

start_time=time.time()

lens = str(sys.argv[1])
innermask = str(sys.argv[2])
radius = str(sys.argv[3])
mag = str(sys.argv[4])
mode = str(sys.argv[5])
conjoined = len(sys.argv) - 6 # total number of arguments including code name, minus the number of ones that are not weights

if conjoined == 1:
    weightin1 = str(sys.argv[6])
if conjoined == 2:
    weightin1 = str(sys.argv[6])
    weightin2 = str(sys.argv[7])
if conjoined == 3:
    weightin1 = str(sys.argv[6])
    weightin2 = str(sys.argv[7])
    weightin3 = str(sys.argv[8])
if conjoined == 4:
    weightin1 = str(sys.argv[6])
    weightin2 = str(sys.argv[7])
    weightin3 = str(sys.argv[8])
    weightin4 = str(sys.argv[9])

print "conjoined:", conjoined
root = "/mfst01a/rusucs/WFI2033/MSwghtratios/"
rootout = "/mfst01a/rusucs/WFI2033/MSkapparesults/"
weightsfile = np.loadtxt(root+'weightedcounts_%s_%s_%sarcsec.lst' % (lens,mode,str(innermask)),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

increment1 = 2 # refers to the E interval from Greene et al. 2014
increment2 = 2
increment3 = 2
increment4 = 4

if lens == "WFI2033":
    constr_gamma = 0.16
    constrwidth_gamma_inf = 0.14
    constrwidth_gamma_sup = 0.18
    filters = "ugrizJHK"

if mag == "23" and radius == "45":
    measured_index = 0 # specifies the column index in weightsfile
    measured_index_inf = 1
    measured_index_sup = 2
if mag == "23" and radius == "120":
    measured_index = 3
    measured_index_inf = 4
    measured_index_sup = 5

if mag == "23":
    if weightin1 == "gal": weight1_index = 4
    if weightin1 == "z": weight1_index = 5
    if weightin1 == "mass": weight1_index = 6
    if weightin1 == "mass2": weight1_index = 7
    if weightin1 == "mass3": weight1_index = 8
    if weightin1 == "oneoverr": weight1_index = 9
    if weightin1 == "zoverr": weight1_index = 10
    if weightin1 == "massoverr": weight1_index = 11
    if weightin1 == "mass2overr": weight1_index = 12
    if weightin1 == "mass3overr": weight1_index = 13
    if weightin1 == "mass2rms": weight1_index = 14
    if weightin1 == "mass3rms": weight1_index = 15
    if weightin1 == "mass2overrrms": weight1_index = 16
    if weightin1 == "mass3overrrms": weight1_index = 17
    if weightin1 == "flexion": weight1_index = 18
    if weightin1 == "tidal": weight1_index = 19
    if weightin1 == "SIS": weight1_index = 20
    if weightin1 == "SIShalo": weight1_index = 21
if conjoined >= 2:
    if mag == "23":
        if weightin2 == "gal": weight2_index = 4
        if weightin2 == "z": weight2_index = 5
        if weightin2 == "mass": weight2_index = 6
        if weightin2 == "mass2": weight2_index = 7
        if weightin2 == "mass3": weight2_index = 8
        if weightin2 == "oneoverr": weight2_index = 9
        if weightin2 == "zoverr": weight2_index = 10
        if weightin2 == "massoverr": weight2_index = 11
        if weightin2 == "mass2overr": weight2_index = 12
        if weightin2 == "mass3overr": weight2_index = 13
        if weightin2 == "mass2rms": weight2_index = 14
        if weightin2 == "mass3rms": weight2_index = 15
        if weightin2 == "mass2overrrms": weight2_index = 16
        if weightin2 == "mass3overrrms": weight2_index = 17
        if weightin2 == "flexion": weight2_index = 18
        if weightin2 == "tidal": weight2_index = 19
        if weightin2 == "SIS": weight2_index = 20
        if weightin2 == "SIShalo": weight2_index = 21
    if conjoined >= 3:
        if mag == "23":
            if weightin3 == "gal": weight3_index = 4
            if weightin3 == "z": weight3_index = 5
            if weightin3 == "mass": weight3_index = 6
            if weightin3 == "mass2": weight3_index = 7
            if weightin3 == "mass3": weight3_index = 8
            if weightin3 == "oneoverr": weight3_index = 9
            if weightin3 == "zoverr": weight3_index = 10
            if weightin3 == "massoverr": weight3_index = 11
            if weightin3 == "mass2overr": weight3_index = 12
            if weightin3 == "mass3overr": weight3_index = 13
            if weightin3 == "mass2rms": weight3_index = 14
            if weightin3 == "mass3rms": weight3_index = 15
            if weightin3 == "mass2overrrms": weight3_index = 16
            if weightin3 == "mass3overrrms": weight3_index = 17
            if weightin3 == "flexion": weight3_index = 18
            if weightin3 == "tidal": weight3_index = 19
            if weightin3 == "SIS": weight3_index = 20
            if weightin3 == "SIShalo": weight3_index = 21
        if conjoined == 4:
            if mag == "23":
                if weightin4 == "gal": weight4_index = 4
                if weightin4 == "z": weight4_index = 5
                if weightin4 == "mass": weight4_index = 6
                if weightin4 == "mass2": weight4_index = 7
                if weightin4 == "mass3": weight4_index = 8
                if weightin4 == "oneoverr": weight4_index = 9
                if weightin4 == "zoverr": weight4_index = 10
                if weightin4 == "massoverr": weight4_index = 11
                if weightin4 == "mass2overr": weight4_index = 12
                if weightin4 == "mass3overr": weight4_index = 13
                if weightin4 == "mass2rms": weight4_index = 14
                if weightin4 == "mass3rms": weight4_index = 15
                if weightin4 == "mass2overrrms": weight4_index = 16
                if weightin4 == "mass3overrrms": weight4_index = 17
                if weightin4 == "flexion": weight4_index = 18
                if weightin4 == "tidal": weight4_index = 19
                if weightin4 == "SIS": weight4_index = 20
                if weightin4 == "SIShalo": weight4_index = 21

constr_gal_meds = weightsfile[measured_index][0]
constrwidth_gal_meds_inf = weightsfile[measured_index_inf][0]
constrwidth_gal_meds_sup = weightsfile[measured_index_sup][0]

constr_z_meds = weightsfile[measured_index][1]
constrwidth_z_meds_inf = weightsfile[measured_index_inf][1]
constrwidth_z_meds_sup = weightsfile[measured_index_sup][1]

constr_mass_meds = weightsfile[measured_index][2]
constrwidth_mass_meds_inf = weightsfile[measured_index_inf][2]
constrwidth_mass_meds_sup = weightsfile[measured_index_sup][2]

constr_mass2_meds = weightsfile[measured_index][3]
constrwidth_mass2_meds_inf = weightsfile[measured_index_inf][3]
constrwidth_mass2_meds_sup = weightsfile[measured_index_sup][3]

constr_mass3_meds = weightsfile[measured_index][4]
constrwidth_mass3_meds_inf = weightsfile[measured_index_inf][4]
constrwidth_mass3_meds_sup = weightsfile[measured_index_sup][4]

constr_oneoverr_meds = weightsfile[measured_index][5]
constrwidth_oneoverr_meds_inf = weightsfile[measured_index_inf][5]
constrwidth_oneoverr_meds_sup = weightsfile[measured_index_sup][5]

constr_zoverr_meds = weightsfile[measured_index][6]
constrwidth_zoverr_meds_inf = weightsfile[measured_index_inf][6]
constrwidth_zoverr_meds_sup = weightsfile[measured_index_sup][6]

constr_massoverr_meds = weightsfile[measured_index][7]
constrwidth_massoverr_meds_inf = weightsfile[measured_index_inf][7]
constrwidth_massoverr_meds_sup = weightsfile[measured_index_sup][7]

constr_mass2overr_meds = weightsfile[measured_index][8]
constrwidth_mass2overr_meds_inf = weightsfile[measured_index_inf][8]
constrwidth_mass2overr_meds_sup = weightsfile[measured_index_sup][8]

constr_mass3overr_meds = weightsfile[measured_index][9]
constrwidth_mass3overr_meds_inf = weightsfile[measured_index_inf][9]
constrwidth_mass3overr_meds_sup = weightsfile[measured_index_sup][9]

constr_mass2rms_meds = weightsfile[measured_index][10]
constrwidth_mass2rms_meds_inf = weightsfile[measured_index_inf][10]
constrwidth_mass2rms_meds_sup = weightsfile[measured_index_sup][10]

constr_mass3rms_meds = weightsfile[measured_index][11]
constrwidth_mass3rms_meds_inf = weightsfile[measured_index_inf][11]
constrwidth_mass3rms_meds_sup = weightsfile[measured_index_sup][11]

constr_mass2overrrms_meds = weightsfile[measured_index][12]
constrwidth_mass2overrrms_meds_inf = weightsfile[measured_index_inf][12]
constrwidth_mass2overrrms_meds_sup = weightsfile[measured_index_sup][12]

constr_mass3overrrms_meds = weightsfile[measured_index][13]
constrwidth_mass3overrrms_meds_inf = weightsfile[measured_index_inf][13]
constrwidth_mass3overrrms_meds_sup = weightsfile[measured_index_sup][13]

constr_flexion_meds = weightsfile[measured_index][14]
constrwidth_flexion_meds_inf = weightsfile[measured_index_inf][14]
constrwidth_flexion_meds_sup = weightsfile[measured_index_sup][14]

constr_tidal_meds = weightsfile[measured_index][15]
constrwidth_tidal_meds_inf = weightsfile[measured_index_inf][15]
constrwidth_tidal_meds_sup = weightsfile[measured_index_sup][15]

constr_SIS_meds = weightsfile[measured_index][16]
constrwidth_SIS_meds_inf = weightsfile[measured_index_inf][16]
constrwidth_SIS_meds_sup = weightsfile[measured_index_sup][16]

constr_SIShalo_meds = weightsfile[measured_index][17]
constrwidth_SIShalo_meds_inf = weightsfile[measured_index_inf][17]
constrwidth_SIShalo_meds_sup = weightsfile[measured_index_sup][17]

if conjoined == 4:
        if weightin4 == "gal": constr_weight4 = constr_gal_meds; constrwidth_weight4_inf = constrwidth_gal_meds_inf; constrwidth_weight4_sup = constrwidth_gal_meds_sup
        if weightin4 == "z": constr_weight4 = constr_z_meds; constrwidth_weight4_inf = constrwidth_z_meds_inf; constrwidth_weight4_sup = constrwidth_z_meds_sup
        if weightin4 == "mass": constr_weight4 = constr_mass_meds; constrwidth_weight4_inf = constrwidth_mass_meds_inf; constrwidth_weight4_sup = constrwidth_mass_meds_sup
        if weightin4 == "mass2": constr_weight4 = constr_mass2_meds; constrwidth_weight4_inf = constrwidth_mass2_meds_inf; constrwidth_weight4_sup = constrwidth_mass2_meds_sup
        if weightin4 == "mass3": constr_weight4 = constr_mass3_meds; constrwidth_weight4_inf = constrwidth_mass3_meds_inf; constrwidth_weight4_sup = constrwidth_mass3_meds_sup
        if weightin4 == "oneoverr": constr_weight4 = constr_oneoverr_meds; constrwidth_weight4_inf = constrwidth_oneoverr_meds_inf; constrwidth_weight4_sup = constrwidth_oneoverr_meds_sup
        if weightin4 == "zoverr": constr_weight4 = constr_zoverr_meds; constrwidth_weight4_inf = constrwidth_zoverr_meds_inf; constrwidth_weight4_sup = constrwidth_zoverr_meds_sup
        if weightin4 == "massoverr": constr_weight4 = constr_massoverr_meds; constrwidth_weight4_inf = constrwidth_massoverr_meds_inf; constrwidth_weight4_sup = constrwidth_massoverr_meds_sup
        if weightin4 == "mass2overr": constr_weight4 = constr_mass2overr_meds; constrwidth_weight4_inf = constrwidth_mass2overr_meds_inf; constrwidth_weight4_sup = constrwidth_mass2overr_meds_sup
        if weightin4 == "mass3overr": constr_weight4 = constr_mass3overr_meds; constrwidth_weight4_inf = constrwidth_mass3overr_meds_inf; constrwidth_weight4_sup = constrwidth_mass3overr_meds_sup
        if weightin4 == "mass2rms": constr_weight4 = constr_mass2rms_meds; constrwidth_weight4_inf = constrwidth_mass2rms_meds_inf; constrwidth_weight4_sup = constrwidth_mass2rms_meds_sup
        if weightin4 == "mass3rms": constr_weight4 = constr_mass3rms_meds; constrwidth_weight4_inf = constrwidth_mass3rms_meds_inf; constrwidth_weight4_sup = constrwidth_mass3rms_meds_sup
        if weightin4 == "mass2overrrms": constr_weight4 = constr_mass2overrrms_meds; constrwidth_weight4_inf = constrwidth_mass2overrrms_meds_inf; constrwidth_weight4_sup = constrwidth_mass2overrrms_meds_sup
        if weightin4 == "mass3overrrms": constr_weight4 = constr_mass3overrrms_meds; constrwidth_weight4_inf = constrwidth_mass3overrrms_meds_inf; constrwidth_weight4_sup = constrwidth_mass3overrrms_meds_sup
        if weightin4 == "flexion": constr_weight4 = constr_flexion_meds; constrwidth_weight4_inf = constrwidth_flexion_meds_inf; constrwidth_weight4_sup = constrwidth_flexion_meds_sup
        if weightin4 == "tidal": constr_weight4 = constr_tidal_meds; constrwidth_weight4_inf = constrwidth_tidal_meds_inf; constrwidth_weight4_sup = constrwidth_tidal_meds_sup
        if weightin4 == "SIS": constr_weight4 = constr_SIS_meds; constrwidth_weight4_inf = constrwidth_SIS_meds_inf; constrwidth_weight4_sup = constrwidth_SIS_meds_sup
        if weightin4 == "SIShalo": constr_weight4 = constr_SIShalo_meds; constrwidth_weight4_inf = constrwidth_SIShalo_meds_inf; constrwidth_weight4_sup = constrwidth_SIShalo_meds_sup
        if weightin4 == "gamma": constr_weight4 = constr_gamma; constrwidth_weight4_inf = constrwidth_gamma_inf; constrwidth_weight4_sup = constrwidth_gamma_sup

if (conjoined == 3) | (conjoined == 4):
        if weightin3 == "gal": constr_weight3 = constr_gal_meds; constrwidth_weight3_inf = constrwidth_gal_meds_inf; constrwidth_weight3_sup = constrwidth_gal_meds_sup
        if weightin3 == "z": constr_weight3 = constr_z_meds; constrwidth_weight3_inf = constrwidth_z_meds_inf; constrwidth_weight3_sup = constrwidth_z_meds_sup
        if weightin3 == "mass": constr_weight3 = constr_mass_meds; constrwidth_weight3_inf = constrwidth_mass_meds_inf; constrwidth_weight3_sup = constrwidth_mass_meds_sup
        if weightin3 == "mass2": constr_weight3 = constr_mass2_meds; constrwidth_weight3_inf = constrwidth_mass2_meds_inf; constrwidth_weight3_sup = constrwidth_mass2_meds_sup
        if weightin3 == "mass3": constr_weight3 = constr_mass3_meds; constrwidth_weight3_inf = constrwidth_mass3_meds_inf; constrwidth_weight3_sup = constrwidth_mass3_meds_sup
        if weightin3 == "oneoverr": constr_weight3 = constr_oneoverr_meds; constrwidth_weight3_inf = constrwidth_oneoverr_meds_inf; constrwidth_weight3_sup = constrwidth_oneoverr_meds_sup
        if weightin3 == "zoverr": constr_weight3 = constr_zoverr_meds; constrwidth_weight3_inf = constrwidth_zoverr_meds_inf; constrwidth_weight3_sup = constrwidth_zoverr_meds_sup
        if weightin3 == "massoverr": constr_weight3 = constr_massoverr_meds; constrwidth_weight3_inf = constrwidth_massoverr_meds_inf; constrwidth_weight3_sup = constrwidth_massoverr_meds_sup
        if weightin3 == "mass2overr": constr_weight3 = constr_mass2overr_meds; constrwidth_weight3_inf = constrwidth_mass2overr_meds_inf; constrwidth_weight3_sup = constrwidth_mass2overr_meds_sup
        if weightin3 == "mass3overr": constr_weight3 = constr_mass3overr_meds; constrwidth_weight3_inf = constrwidth_mass3overr_meds_inf; constrwidth_weight3_sup = constrwidth_mass3overr_meds_sup
        if weightin3 == "mass2rms": constr_weight3 = constr_mass2rms_meds; constrwidth_weight3_inf = constrwidth_mass2rms_meds_inf; constrwidth_weight3_sup = constrwidth_mass2rms_meds_sup
        if weightin3 == "mass3rms": constr_weight3 = constr_mass3rms_meds; constrwidth_weight3_inf = constrwidth_mass3rms_meds_inf; constrwidth_weight3_sup = constrwidth_mass3rms_meds_sup
        if weightin3 == "mass2overrrms": constr_weight3 = constr_mass2overrrms_meds; constrwidth_weight3_inf = constrwidth_mass2overrrms_meds_inf; constrwidth_weight3_sup = constrwidth_mass2overrrms_meds_sup
        if weightin3 == "mass3overrrms": constr_weight3 = constr_mass3overrrms_meds; constrwidth_weight3_inf = constrwidth_mass3overrrms_meds_inf; constrwidth_weight3_sup = constrwidth_mass3overrrms_meds_sup
        if weightin3 == "flexion": constr_weight3 = constr_flexion_meds; constrwidth_weight3_inf = constrwidth_flexion_meds_inf; constrwidth_weight3_sup = constrwidth_flexion_meds_sup
        if weightin3 == "tidal": constr_weight3 = constr_tidal_meds; constrwidth_weight3_inf = constrwidth_tidal_meds_inf; constrwidth_weight3_sup = constrwidth_tidal_meds_sup
        if weightin3 == "SIS": constr_weight3 = constr_SIS_meds; constrwidth_weight3_inf = constrwidth_SIS_meds_inf; constrwidth_weight3_sup = constrwidth_SIS_meds_sup
        if weightin3 == "SIShalo": constr_weight3 = constr_SIShalo_meds; constrwidth_weight3_inf = constrwidth_SIShalo_meds_inf; constrwidth_weight3_sup = constrwidth_SIShalo_meds_sup
        if weightin3 == "gamma": constr_weight3 = constr_gamma; constrwidth_weight3_inf = constrwidth_gamma_inf; constrwidth_weight3_sup = constrwidth_gamma_sup

if (conjoined == 2) | (conjoined == 3) | (conjoined == 4):
        if weightin2 == "gal": constr_weight2 = constr_gal_meds; constrwidth_weight2_inf = constrwidth_gal_meds_inf; constrwidth_weight2_sup = constrwidth_gal_meds_sup
        if weightin2 == "z": constr_weight2 = constr_z_meds; constrwidth_weight2_inf = constrwidth_z_meds_inf; constrwidth_weight2_sup = constrwidth_z_meds_sup
        if weightin2 == "mass": constr_weight2 = constr_mass_meds; constrwidth_weight2_inf = constrwidth_mass_meds_inf; constrwidth_weight2_sup = constrwidth_mass_meds_sup
        if weightin2 == "mass2": constr_weight2 = constr_mass2_meds; constrwidth_weight2_inf = constrwidth_mass2_meds_inf; constrwidth_weight2_sup = constrwidth_mass2_meds_sup
        if weightin2 == "mass3": constr_weight2 = constr_mass3_meds; constrwidth_weight2_inf = constrwidth_mass3_meds_inf; constrwidth_weight2_sup = constrwidth_mass3_meds_sup
        if weightin2 == "oneoverr": constr_weight2 = constr_oneoverr_meds; constrwidth_weight2_inf = constrwidth_oneoverr_meds_inf; constrwidth_weight2_sup = constrwidth_oneoverr_meds_sup
        if weightin2 == "zoverr": constr_weight2 = constr_zoverr_meds; constrwidth_weight2_inf = constrwidth_zoverr_meds_inf; constrwidth_weight2_sup = constrwidth_zoverr_meds_sup
        if weightin2 == "massoverr": constr_weight2 = constr_massoverr_meds; constrwidth_weight2_inf = constrwidth_massoverr_meds_inf; constrwidth_weight2_sup = constrwidth_massoverr_meds_sup
        if weightin2 == "mass2overr": constr_weight2 = constr_mass2overr_meds; constrwidth_weight2_inf = constrwidth_mass2overr_meds_inf; constrwidth_weight2_sup = constrwidth_mass2overr_meds_sup
        if weightin2 == "mass3overr": constr_weight2 = constr_mass3overr_meds; constrwidth_weight2_inf = constrwidth_mass3overr_meds_inf; constrwidth_weight2_sup = constrwidth_mass3overr_meds_sup
        if weightin2 == "mass2rms": constr_weight2 = constr_mass2rms_meds; constrwidth_weight2_inf = constrwidth_mass2rms_meds_inf; constrwidth_weight2_sup = constrwidth_mass2rms_meds_sup
        if weightin2 == "mass3rms": constr_weight2 = constr_mass3rms_meds; constrwidth_weight2_inf = constrwidth_mass3rms_meds_inf; constrwidth_weight2_sup = constrwidth_mass3rms_meds_sup
        if weightin2 == "mass2overrrms": constr_weight2 = constr_mass2overrrms_meds; constrwidth_weight2_inf = constrwidth_mass2overrrms_meds_inf; constrwidth_weight2_sup = constrwidth_mass2overrrms_meds_sup
        if weightin2 == "mass3overrrms": constr_weight2 = constr_mass3overrrms_meds; constrwidth_weight2_inf = constrwidth_mass3overrrms_meds_inf; constrwidth_weight2_sup = constrwidth_mass3overrrms_meds_sup
        if weightin2 == "flexion": constr_weight2 = constr_flexion_meds; constrwidth_weight2_inf = constrwidth_flexion_meds_inf; constrwidth_weight2_sup = constrwidth_flexion_meds_sup
        if weightin2 == "tidal": constr_weight2 = constr_tidal_meds; constrwidth_weight2_inf = constrwidth_tidal_meds_inf; constrwidth_weight2_sup = constrwidth_tidal_meds_sup
        if weightin2 == "SIS": constr_weight2 = constr_SIS_meds; constrwidth_weight2_inf = constrwidth_SIS_meds_inf; constrwidth_weight2_sup = constrwidth_SIS_meds_sup
        if weightin2 == "SIShalo": constr_weight2 = constr_SIShalo_meds; constrwidth_weight2_inf = constrwidth_SIShalo_meds_inf; constrwidth_weight2_sup = constrwidth_SIShalo_meds_sup
        if weightin2 == "gamma": constr_weight2 = constr_gamma; constrwidth_weight2_inf = constrwidth_gamma_inf; constrwidth_weight2_sup = constrwidth_gamma_sup

if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4):
        if weightin1 == "gal": constr_weight1 = constr_gal_meds; constrwidth_weight1_inf = constrwidth_gal_meds_inf; constrwidth_weight1_sup = constrwidth_gal_meds_sup
        if weightin1 == "z": constr_weight1 = constr_z_meds; constrwidth_weight1_inf = constrwidth_z_meds_inf; constrwidth_weight1_sup = constrwidth_z_meds_sup
        if weightin1 == "mass": constr_weight1 = constr_mass_meds; constrwidth_weight1_inf = constrwidth_mass_meds_inf; constrwidth_weight1_sup = constrwidth_mass_meds_sup
        if weightin1 == "mass2": constr_weight1 = constr_mass2_meds; constrwidth_weight1_inf = constrwidth_mass2_meds_inf; constrwidth_weight1_sup = constrwidth_mass2_meds_sup
        if weightin1 == "mass3": constr_weight1 = constr_mass3_meds; constrwidth_weight1_inf = constrwidth_mass3_meds_inf; constrwidth_weight1_sup = constrwidth_mass3_meds_sup
        if weightin1 == "oneoverr": constr_weight1 = constr_oneoverr_meds; constrwidth_weight1_inf = constrwidth_oneoverr_meds_inf; constrwidth_weight1_sup = constrwidth_oneoverr_meds_sup
        if weightin1 == "zoverr": constr_weight1 = constr_zoverr_meds; constrwidth_weight1_inf = constrwidth_zoverr_meds_inf; constrwidth_weight1_sup = constrwidth_zoverr_meds_sup
        if weightin1 == "massoverr": constr_weight1 = constr_massoverr_meds; constrwidth_weight1_inf = constrwidth_massoverr_meds_inf; constrwidth_weight1_sup = constrwidth_massoverr_meds_sup
        if weightin1 == "mass2overr": constr_weight1 = constr_mass2overr_meds; constrwidth_weight1_inf = constrwidth_mass2overr_meds_inf; constrwidth_weight1_sup = constrwidth_mass2overr_meds_sup
        if weightin1 == "mass3overr": constr_weight1 = constr_mass3overr_meds; constrwidth_weight1_inf = constrwidth_mass3overr_meds_inf; constrwidth_weight1_sup = constrwidth_mass3overr_meds_sup
        if weightin1 == "mass2rms": constr_weight1 = constr_mass2rms_meds; constrwidth_weight1_inf = constrwidth_mass2rms_meds_inf; constrwidth_weight1_sup = constrwidth_mass2rms_meds_sup
        if weightin1 == "mass3rms": constr_weight1 = constr_mass3rms_meds; constrwidth_weight1_inf = constrwidth_mass3rms_meds_inf; constrwidth_weight1_sup = constrwidth_mass3rms_meds_sup
        if weightin1 == "mass2overrrms": constr_weight1 = constr_mass2overrrms_meds; constrwidth_weight1_inf = constrwidth_mass2overrrms_meds_inf; constrwidth_weight1_sup = constrwidth_mass2overrrms_meds_sup
        if weightin1 == "mass3overrrms": constr_weight1 = constr_mass3overrrms_meds; constrwidth_weight1_inf = constrwidth_mass3overrrms_meds_inf; constrwidth_weight1_sup = constrwidth_mass3overrrms_meds_sup
        if weightin1 == "flexion": constr_weight1 = constr_flexion_meds; constrwidth_weight1_inf = constrwidth_flexion_meds_inf; constrwidth_weight1_sup = constrwidth_flexion_meds_sup
        if weightin1 == "tidal": constr_weight1 = constr_tidal_meds; constrwidth_weight1_inf = constrwidth_tidal_meds_inf; constrwidth_weight1_sup = constrwidth_tidal_meds_sup
        if weightin1 == "SIS": constr_weight1 = constr_SIS_meds; constrwidth_weight1_inf = constrwidth_SIS_meds_inf; constrwidth_weight1_sup = constrwidth_SIS_meds_sup
        if weightin1 == "SIShalo": constr_weight1 = constr_SIShalo_meds; constrwidth_weight1_inf = constrwidth_SIShalo_meds_inf; constrwidth_weight1_sup = constrwidth_SIShalo_meds_sup
        if weightin1 == "gamma": constr_weight1 = constr_gamma; constrwidth_weight1_inf = constrwidth_gamma_inf; constrwidth_weight1_sup = constrwidth_gamma_sup

print "Reading..."

#if mode == "sum":
    #str1 = ""
if mode == "meds":
    str1 = "med"

if conjoined == 4:
    output = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,weightin3,weightin4,mag,radius,mode,increment1,increment2,increment3,increment4)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,weightin3,weightin4,mag,radius,mode,increment1,increment2,increment3,increment4)
if conjoined == 3:
    output = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_%s_increments%s_%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,weightin3,mag,radius,mode,increment1,increment2,increment3)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,weightin3,mag,radius,mode,increment1,increment2,increment3)
if conjoined == 2:
    output = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_increments%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,mag,radius,mode,increment1,increment2)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_%s_LOS_increments%s_%s.cat' % (rootout,lens,innermask,weightin1,weightin2,mag,radius,mode,increment1,increment2)
if conjoined == 1:
    output = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_increments%s.cat' % (rootout,lens,innermask,weightin1,mag,radius,mode,increment1)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta_%s_%s_%s_%s_LOS_increments%s.cat' % (rootout,lens,innermask,weightin1,mag,radius,mode,increment1)

if conjoined == 1:
    ''' Here I only read the columns of interest, without kappa, for ugriz, in order to find the medians of their values over the whole MS.'''
    med1 = np.zeros(8)
    for j in range(8):
      for i in range(8):
        if weightin1 != "gamma":
            weight1_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=[weight1_index], unpack=True)
            if i == 0:
                weight1 = weight1_
            else:
                weight1 = np.append(weight1,weight1_)
        else:
            weight1_1_,weight1_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=[2,3], unpack=True)
            if i == 0:
                weight1_1 = weight1_1_
                weight1_2 = weight1_2_
            else:
                weight1_1 = np.append(weight1_1,weight1_1_)
                weight1_2 = np.append(weight1_2,weight1_2_)
      if weightin1 != "gamma":
        med1[j] = np.median(weight1)
      else:
        med1[j] = np.median(np.sqrt(weight1_1**2 + weight1_2**2))
      print j
    med_weight1 = np.mean(med1) # throughout the code I use med_weight1 when computing intervals, following Green et al. For this, weight1 should always refer to simple galaxy number counts
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))]) # absolute number, e.g. of galaxies within the lower width
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    if weightin1 == "gamma":
        constr_weight1 = constr_weight1 / med_weight1 # for gamma, measured shear divided by the median value of shear in MS; this turns it into an overdensity, like the other weights
    
    ''' Here I read ugrizJHK, converting weighted counts into overdensities, and recording the kappa values only for overdensities satisfying the constraint. I consider the full range of the constraint.'''
    for j in range(8):
      for i in range(8):
        if weightin1 != "gamma":
            kappa_, weight1_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index), unpack=True)
            weight1_ = weight1_ / med_weight1
        else:
            kappa_, gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,2,3), unpack=True)
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = gamma / med_weight1_
        weight = np.copy(weight1_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        if (i == 0) and (j == 0):
            kappa = kappa_
            weight1 = weight1_
        else:
            kappa = np.append(kappa,kappa_)
            weight1 = np.append(weight1,weight1_)
      print j

if conjoined == 2:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=(weight1_index,weight2_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
        else:
            weight1_,weight2_1_,weight2_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=[weight1_index,1,2], unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2_1 = weight2_1_
                weight2_2 = weight2_2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2_1 = np.append(weight2_1,weight2_1_)
                weight2_2 = np.append(weight2_2,weight2_2_)
      if weightin2 != "gamma":
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
      print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    if weightin2 == "gamma":
        constr_weight2 = constr_weight2 / med_weight2

    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,weight2_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
        else:
            kappa_, weight1_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,2,3), unpack=True)
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
        weight = np.copy(weight1_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        del weight
        if (i == 0) and (j == 0):
            kappa = kappa_
            weight1 = weight1_
            weight2 = weight2_
        else:
            kappa = np.append(kappa,kappa_)
            weight1 = np.append(weight1,weight1_)
            weight2 = np.append(weight2,weight2_)
      print j

if conjoined == 3:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=(weight1_index,weight2_index,weight3_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
                weight3 = weight3_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
                weight3 = np.append(weight3,weight3_)
        else:
            weight1_,weight2_1_,weight2_2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=(weight1_index,1,2,weight3_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2_1 = weight2_1_
                weight2_2 = weight2_2_
                weight3 = weight3_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2_1 = np.append(weight2_1,weight2_1_)
                weight2_2 = np.append(weight2_2,weight2_2_)
                weight3 = np.append(weight3,weight3_)
      if weightin2 != "gamma":
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
        med3[j] = np.median(weight3)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
        med3[j] = np.median(weight3)
      print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    E_w3_inf = np.max([1, round(med_weight1 * (constr_weight3 - constrwidth_weight3_inf))])
    E_w3_sup = np.max([1, round(med_weight1 * (-constr_weight3 + constrwidth_weight3_sup))])
    if weightin2 == "gamma":
        constr_weight2 = constr_weight2 / med_weight2

    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,weight2_index,weight3_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
        else:
            kappa_, weight1_,weight3_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,weight3_index,2,3), unpack=True)
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
            weight3_ = weight3_ / med_weight3
        weight = np.copy(weight1_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        del weight
        weight = np.copy(weight3_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        del weight
        if (i == 0) and (j == 0):
            kappa = kappa_
            weight1 = weight1_
            weight2 = weight2_
            weight3 = weight3_
        else:
            kappa = np.append(kappa,kappa_)
            weight1 = np.append(weight1,weight1_)
            weight2 = np.append(weight2,weight2_)
            weight3 = np.append(weight3,weight3_)
      print j

if conjoined == 4:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    med4 = np.zeros(8)
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=(weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
                weight3 = weight3_
                weight4 = weight4_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
                weight3 = np.append(weight3,weight3_)
                weight4 = np.append(weight4,weight4_)
        else:
            weight1_,weight2_1_,weight2_2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,str(j),str(i),radius,innermask), usecols=(weight1_index,1,2,weight3_index,weight4_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2_1 = weight2_1_
                weight2_2 = weight2_2_
                weight3 = weight3_
                weight4 = weight4_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2_1 = np.append(weight2_1,weight2_1_)
                weight2_2 = np.append(weight2_2,weight2_2_)
                weight3 = np.append(weight3,weight3_)
                weight4 = np.append(weight4,weight4_)
      if weightin2 != "gamma":
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
        med3[j] = np.median(weight3)
        med4[j] = np.median(weight4)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
        med3[j] = np.median(weight3)
        med4[j] = np.median(weight4)
      print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    med_weight4 = np.mean(med4)
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    E_w3_inf = np.max([1, round(med_weight1 * (constr_weight3 - constrwidth_weight3_inf))])
    E_w3_sup = np.max([1, round(med_weight1 * (-constr_weight3 + constrwidth_weight3_sup))])
    E_w4_inf = np.max([1, round(med_weight1 * (constr_weight4 - constrwidth_weight4_inf))])
    E_w4_sup = np.max([1, round(med_weight1 * (-constr_weight4 + constrwidth_weight4_sup))])
    if weightin2 == "gamma":
        constr_weight2 = constr_weight2 / med_weight2

    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        else:
            kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_WFI2033_GGL_los_8_%s_%s_%s_%sarcsecinnermsk.cat" % (root,str1,filters,str(j),str(i),radius,innermask), usecols=(1,weight1_index,weight3_index,weight4_index,2,3), unpack=True)
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        weight = np.copy(weight1_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        del weight
        weight = np.copy(weight3_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        del weight
        if (i == 0) and (j == 0):
            kappa = kappa_
            weight1 = weight1_
            weight2 = weight2_
            weight3 = weight3_
            weight4 = weight4_
        else:
            kappa = np.append(kappa,kappa_)
            weight1 = np.append(weight1,weight1_)
            weight2 = np.append(weight2,weight2_)
            weight3 = np.append(weight3,weight3_)
            weight4 = np.append(weight4,weight4_)
      print j

print(" Read in %s seconds" % (time.time() - start_time))

gauss = sp.stats.norm(0, 1)
start1 = time.time()
LOS = 0

if conjoined == 4:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1): # use as specific value
        for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2):
            for E3 in np.arange(-limsigma * E_w3_inf, limsigma * E_w3_sup + 1, increment3):
                for E4 in np.arange(-limsigma * E_w4_inf, limsigma * E_w4_sup + 1, increment4):
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") ", "E3 = ", E3, "in (", -limsigma * E_w3_inf, ",", limsigma * E_w3_sup, ") ", "E4 = ", E4, "in (", -limsigma * E_w4_inf, ",", limsigma * E_w4_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight1 >= round(constr_weight4 * med_weight1) + E4 - increment4/2.0) & (weight4 * med_weight1 < round(constr_weight4 * med_weight1) + E4 + increment4/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf)
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        if E2 < 0: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_inf)
                        else: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_sup)
                        if E3 < 0: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_inf)
                        else: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_sup)
                        if E4 < 0: gauss_factorE4 = gauss.pdf(float(E4)/E_w4_inf)
                        else: gauss_factorE4 = gauss.pdf(float(E4)/E_w4_sup)
                        kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 * gauss_factorE3 * gauss_factorE4 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained
                    LOS = LOS + data.size

if conjoined == 3:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1):
        for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2):
            for E3 in np.arange(-limsigma * E_w3_inf, limsigma * E_w3_sup + 1, increment3):
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") ", "E3 = ", E3, "in (", -limsigma * E_w3_inf, ",", limsigma * E_w3_sup, ") "#, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf)
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        if E2 < 0: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_inf)
                        else: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_sup)
                        if E3 < 0: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_inf)
                        else: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_sup)
                        kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 * gauss_factorE3 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained
                    LOS = LOS + data.size

if conjoined == 2:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1):
        for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2):
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf)
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        if E2 < 0: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_inf)
                        else: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_sup)
                        kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained
                    LOS = LOS + data.size

if conjoined == 1:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1):
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf) # for asymmetric limits, implement a gaussian on each side
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained # I tested that this addition works correctly
                    LOS = LOS + data.size

np.savetxt(output,unbiased_kappa_constrained,fmt='%s',delimiter='\t',newline='\n')
np.savetxt(outputLOS,np.array([LOS]),fmt='%s',delimiter='\t',newline='\n')
print(" time for computing kappa %s seconds" % (time.time() - start1))

if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4):
    print "increment1 = ", increment1
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4):
    print "increment2 = ", increment2
if (conjoined == 3) | (conjoined == 4):
    print "increment3 = ", increment3
if conjoined == 4:
    print "increment4 = ", increment4

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
