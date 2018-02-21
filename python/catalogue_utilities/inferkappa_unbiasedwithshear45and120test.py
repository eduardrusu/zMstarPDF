# CE Rusu Feb 19 2018
# Run as python /lfs08/rusucs/code/inferkappa_unbiasedwithshear45and120.py WFI2033 -1.0 -1.0 yes fiducial 5 23 meds 120_gal 120_gamma 120_oneoverr 45_gal 45_oneoverr
# do not use more than 5 constraints in total, of which maximum 4 can refer to the same radius; do not mix order of 45_ and 120_; e.g.: 45_gal 45_oneoverr 120_gamma correct, but 45_gal 120_gamma 45_oneoverr incorrect
# the code currently works for maglim 23 (WFI2033)
# Description of arguments: inferkappa_unbiasedwithshear.py lens radius maglim innermask sum/meds gal list_of_weight_constraints
# for each redius, weight1 should always be "gal", in order to use the galaxy counts when correcting the bias due to different LOS
# the code is written such that, if shear is used as overdensity, it should be the second weight used in either radius

import sys
import os
from os import system
import scipy as sp
from scipy.stats import norm
import numpy as np
import time

start_time=time.time()

lens = str(sys.argv[1])
zinf = str(sys.argv[2])
zsup = str(sys.argv[3])
handpicked = str(sys.argv[4])
other = str(sys.argv[5]) # refers to an optional suffix for the shear constraint
innermask = str(sys.argv[6])
mag = str(sys.argv[7])
mode = str(sys.argv[8])
conjoined = len(sys.argv) - 9 # total number of arguments including code name, minus the number of ones that are not weights

if handpicked == 'yes': handpickedstr = '_handpicked'
else: handpickedstr = ''

if conjoined == 1:
    weightin1 = str(sys.argv[9])
if conjoined == 2:
    weightin1 = str(sys.argv[9])
    weightin2 = str(sys.argv[10])
if conjoined == 3:
    weightin1 = str(sys.argv[9])
    weightin2 = str(sys.argv[10])
    weightin3 = str(sys.argv[11])
if conjoined == 4:
    weightin1 = str(sys.argv[9])
    weightin2 = str(sys.argv[10])
    weightin3 = str(sys.argv[11])
    weightin4 = str(sys.argv[12])
if conjoined == 5:
    weightin1 = str(sys.argv[9])
    weightin2 = str(sys.argv[10])
    weightin3 = str(sys.argv[11])
    weightin4 = str(sys.argv[12])
    weightin5 = str(sys.argv[13])

print "conjoined:", conjoined
#root = "/lfs08/rusucs/%s/MSwghtratios/" % lens
root = "/Volumes/LaCieSubaru/MSweights/"
#rootout = "/lfs08/rusucs/%s/MSkapparesults/" % lens
rootout = "/Volumes/LaCieSubaru/kapparesults/"
#weightsfile = np.loadtxt(root+'weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(lens,mode,innermask,handpickedstr,zinf,zsup),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
weightsfile = np.loadtxt('/Users/cerusu/Dropbox/Davis_work/code/WFI2033/weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(lens,mode,innermask,handpickedstr,zinf,zsup),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

increment1 = 2 # refers to the E interval from Greene et al. 2014
increment2 = 2
increment3 = 2
increment4 = 2
increment5 = 2

# define the shear constraints
if lens == "WFI2033":
    if other == 'fiducial' and handpicked == 'no' and float(zsup) < 0 and innermask == '5':
        constr_gamma = 0.154
        constrwidth_gamma_inf = 0.139
        constrwidth_gamma_sup = 0.169
    if other == 'chameleon' and handpicked == 'no' and float(zsup) < 0 and innermask == '5':
        constr_gamma = 0.193
        constrwidth_gamma_inf = 0.178
        constrwidth_gamma_sup = 0.208
    if other == 'fiducial' and (handpicked == 'yes' or innermask == '15') and float(zsup) < 0:
        constr_gamma = 0.09
        constrwidth_gamma_inf = 0.075
        constrwidth_gamma_sup = 0.105
    if other == 'fiducial' and (handpicked == 'yes' or innermask == '15') and float(zsup) < 0:
        constr_gamma = 0.09
        constrwidth_gamma_inf = 0.075
        constrwidth_gamma_sup = 0.105
    filters = "ugrizJHK"
    print 'shear: ',constr_gamma

# declare which weights to read
if mag == "23": # the only option currently working
    measured_index45 = 0 # specifies the column index in weightsfile
    measured_index_inf45 = 1
    measured_index_sup45 = 2
    measured_index120 = 3
    measured_index_inf120 = 4
    measured_index_sup120 = 5

def declareweight(weightin): # column index in the files created by kappamed_insertstarsnobeta.py
    if weightin.split('_')[1] == "gal": weight_index = 4
    if weightin.split('_')[1] == "z": weight_index = 5
    if weightin.split('_')[1] == "mass": weight_index = 6
    if weightin.split('_')[1] == "mass2": weight_index = 7
    if weightin.split('_')[1] == "mass3": weight_index = 8
    if weightin.split('_')[1] == "oneoverr": weight_index = 9
    if weightin.split('_')[1] == "zoverr": weight_index = 10
    if weightin.split('_')[1] == "massoverr": weight_index = 11
    if weightin.split('_')[1] == "mass2overr": weight_index = 12
    if weightin.split('_')[1] == "mass3overr": weight_index = 13
    if weightin.split('_')[1] == "mass2rms": weight_index = 14
    if weightin.split('_')[1] == "mass3rms": weight_index = 15
    if weightin.split('_')[1] == "mass2overrrms": weight_index = 16
    if weightin.split('_')[1] == "mass3overrrms": weight_index = 17
    if weightin.split('_')[1] == "flexion": weight_index = 18
    if weightin.split('_')[1] == "tidal": weight_index = 19
    if weightin.split('_')[1] == "SIS": weight_index = 20
    if weightin.split('_')[1] == "SIShalo": weight_index = 21
    if weightin.split('_')[1] == "gamma": weight_index = None
    return weight_index

if mag == "23":
    weight1_index = declareweight(weightin1)
if conjoined >= 2:
    if mag == "23":
        weight2_index = declareweight(weightin2)
    if conjoined >= 3:
        if mag == "23":
            weight3_index = declareweight(weightin3)
        if conjoined == 4:
            if mag == "23":
                weight4_index = declareweight(weightin4)
            if conjoined == 5:
                if mag == "23":
                    weight5_index = declareweight(weightin5)

# read weight constraints
constr_gal_meds45 = weightsfile[measured_index45][0]
constrwidth_gal_meds_inf45 = weightsfile[measured_index_inf45][0]
constrwidth_gal_meds_sup45 = weightsfile[measured_index_sup45][0]

constr_z_meds45 = weightsfile[measured_index45][1]
constrwidth_z_meds_inf45 = weightsfile[measured_index_inf45][1]
constrwidth_z_meds_sup45 = weightsfile[measured_index_sup45][1]

constr_mass_meds45 = weightsfile[measured_index45][2]
constrwidth_mass_meds_inf45 = weightsfile[measured_index_inf45][2]
constrwidth_mass_meds_sup45 = weightsfile[measured_index_sup45][2]

constr_mass2_meds45 = weightsfile[measured_index45][3]
constrwidth_mass2_meds_inf45 = weightsfile[measured_index_inf45][3]
constrwidth_mass2_meds_sup45 = weightsfile[measured_index_sup45][3]

constr_mass3_meds45 = weightsfile[measured_index45][4]
constrwidth_mass3_meds_inf45 = weightsfile[measured_index_inf45][4]
constrwidth_mass3_meds_sup45 = weightsfile[measured_index_sup45][4]

constr_oneoverr_meds45 = weightsfile[measured_index45][5]
constrwidth_oneoverr_meds_inf45 = weightsfile[measured_index_inf45][5]
constrwidth_oneoverr_meds_sup45 = weightsfile[measured_index_sup45][5]

constr_zoverr_meds45 = weightsfile[measured_index45][6]
constrwidth_zoverr_meds_inf45 = weightsfile[measured_index_inf45][6]
constrwidth_zoverr_meds_sup45 = weightsfile[measured_index_sup45][6]

constr_massoverr_meds45 = weightsfile[measured_index45][7]
constrwidth_massoverr_meds_inf45 = weightsfile[measured_index_inf45][7]
constrwidth_massoverr_meds_sup45 = weightsfile[measured_index_sup45][7]

constr_mass2overr_meds45 = weightsfile[measured_index45][8]
constrwidth_mass2overr_meds_inf45 = weightsfile[measured_index_inf45][8]
constrwidth_mass2overr_meds_sup45 = weightsfile[measured_index_sup45][8]

constr_mass3overr_meds45 = weightsfile[measured_index45][9]
constrwidth_mass3overr_meds_inf45 = weightsfile[measured_index_inf45][9]
constrwidth_mass3overr_meds_sup45 = weightsfile[measured_index_sup45][9]

constr_mass2rms_meds45 = weightsfile[measured_index45][10]
constrwidth_mass2rms_meds_inf45 = weightsfile[measured_index_inf45][10]
constrwidth_mass2rms_meds_sup45 = weightsfile[measured_index_sup45][10]

constr_mass3rms_meds45 = weightsfile[measured_index45][11]
constrwidth_mass3rms_meds_inf45 = weightsfile[measured_index_inf45][11]
constrwidth_mass3rms_meds_sup45 = weightsfile[measured_index_sup45][11]

constr_mass2overrrms_meds45 = weightsfile[measured_index45][12]
constrwidth_mass2overrrms_meds_inf45 = weightsfile[measured_index_inf45][12]
constrwidth_mass2overrrms_meds_sup45 = weightsfile[measured_index_sup45][12]

constr_mass3overrrms_meds45 = weightsfile[measured_index45][13]
constrwidth_mass3overrrms_meds_inf45 = weightsfile[measured_index_inf45][13]
constrwidth_mass3overrrms_meds_sup45 = weightsfile[measured_index_sup45][13]

constr_flexion_meds45 = weightsfile[measured_index45][14]
constrwidth_flexion_meds_inf45 = weightsfile[measured_index_inf45][14]
constrwidth_flexion_meds_sup45 = weightsfile[measured_index_sup45][14]

constr_tidal_meds45 = weightsfile[measured_index45][15]
constrwidth_tidal_meds_inf45 = weightsfile[measured_index_inf45][15]
constrwidth_tidal_meds_sup45 = weightsfile[measured_index_sup45][15]

constr_SIS_meds45 = weightsfile[measured_index45][16]
constrwidth_SIS_meds_inf45 = weightsfile[measured_index_inf45][16]
constrwidth_SIS_meds_sup45 = weightsfile[measured_index_sup45][16]

constr_SIShalo_meds45 = weightsfile[measured_index45][17]
constrwidth_SIShalo_meds_inf45 = weightsfile[measured_index_inf45][17]
constrwidth_SIShalo_meds_sup45 = weightsfile[measured_index_sup45][17]


constr_gal_meds120 = weightsfile[measured_index120][0]
constrwidth_gal_meds_inf120 = weightsfile[measured_index_inf120][0]
constrwidth_gal_meds_sup120 = weightsfile[measured_index_sup120][0]

constr_z_meds120 = weightsfile[measured_index120][1]
constrwidth_z_meds_inf120 = weightsfile[measured_index_inf120][1]
constrwidth_z_meds_sup120 = weightsfile[measured_index_sup120][1]

constr_mass_meds120 = weightsfile[measured_index120][2]
constrwidth_mass_meds_inf120 = weightsfile[measured_index_inf120][2]
constrwidth_mass_meds_sup120 = weightsfile[measured_index_sup120][2]

constr_mass2_meds120 = weightsfile[measured_index120][3]
constrwidth_mass2_meds_inf120 = weightsfile[measured_index_inf120][3]
constrwidth_mass2_meds_sup120 = weightsfile[measured_index_sup120][3]

constr_mass3_meds120 = weightsfile[measured_index120][4]
constrwidth_mass3_meds_inf120 = weightsfile[measured_index_inf120][4]
constrwidth_mass3_meds_sup120 = weightsfile[measured_index_sup120][4]

constr_oneoverr_meds120 = weightsfile[measured_index120][5]
constrwidth_oneoverr_meds_inf120 = weightsfile[measured_index_inf120][5]
constrwidth_oneoverr_meds_sup120 = weightsfile[measured_index_sup120][5]

constr_zoverr_meds120 = weightsfile[measured_index120][6]
constrwidth_zoverr_meds_inf120 = weightsfile[measured_index_inf120][6]
constrwidth_zoverr_meds_sup120 = weightsfile[measured_index_sup120][6]

constr_massoverr_meds120 = weightsfile[measured_index120][7]
constrwidth_massoverr_meds_inf120 = weightsfile[measured_index_inf120][7]
constrwidth_massoverr_meds_sup120 = weightsfile[measured_index_sup120][7]

constr_mass2overr_meds120 = weightsfile[measured_index120][8]
constrwidth_mass2overr_meds_inf120 = weightsfile[measured_index_inf120][8]
constrwidth_mass2overr_meds_sup120 = weightsfile[measured_index_sup120][8]

constr_mass3overr_meds120 = weightsfile[measured_index120][9]
constrwidth_mass3overr_meds_inf120 = weightsfile[measured_index_inf120][9]
constrwidth_mass3overr_meds_sup120 = weightsfile[measured_index_sup120][9]

constr_mass2rms_meds120 = weightsfile[measured_index120][10]
constrwidth_mass2rms_meds_inf120 = weightsfile[measured_index_inf120][10]
constrwidth_mass2rms_meds_sup120 = weightsfile[measured_index_sup120][10]

constr_mass3rms_meds120 = weightsfile[measured_index120][11]
constrwidth_mass3rms_meds_inf120 = weightsfile[measured_index_inf120][11]
constrwidth_mass3rms_meds_sup120 = weightsfile[measured_index_sup120][11]

constr_mass2overrrms_meds120 = weightsfile[measured_index120][12]
constrwidth_mass2overrrms_meds_inf120 = weightsfile[measured_index_inf120][12]
constrwidth_mass2overrrms_meds_sup120 = weightsfile[measured_index_sup120][12]

constr_mass3overrrms_meds120 = weightsfile[measured_index120][13]
constrwidth_mass3overrrms_meds_inf120 = weightsfile[measured_index_inf120][13]
constrwidth_mass3overrrms_meds_sup120 = weightsfile[measured_index_sup120][13]

constr_flexion_meds120 = weightsfile[measured_index120][14]
constrwidth_flexion_meds_inf120 = weightsfile[measured_index_inf120][14]
constrwidth_flexion_meds_sup120 = weightsfile[measured_index_sup120][14]

constr_tidal_meds120 = weightsfile[measured_index120][15]
constrwidth_tidal_meds_inf120 = weightsfile[measured_index_inf120][15]
constrwidth_tidal_meds_sup120 = weightsfile[measured_index_sup120][15]

constr_SIS_meds120 = weightsfile[measured_index120][16]
constrwidth_SIS_meds_inf120 = weightsfile[measured_index_inf120][16]
constrwidth_SIS_meds_sup120 = weightsfile[measured_index_sup120][16]

constr_SIShalo_meds120 = weightsfile[measured_index120][17]
constrwidth_SIShalo_meds_inf120 = weightsfile[measured_index_inf120][17]
constrwidth_SIShalo_meds_sup120 = weightsfile[measured_index_sup120][17]

def declareweight(weightin):
    if weightin.split('_')[0] == "45":
        if weightin.split('_')[1] == "gal": constr_weight = constr_gal_meds45; constrwidth_weight_inf = constrwidth_gal_meds_inf45; constrwidth_weight_sup = constrwidth_gal_meds_sup45
        if weightin.split('_')[1] == "z": constr_weight = constr_z_meds45; constrwidth_weight_inf = constrwidth_z_meds_inf45; constrwidth_weight_sup = constrwidth_z_meds_sup45
        if weightin.split('_')[1] == "mass": constr_weight = constr_mass_meds45; constrwidth_weight_inf = constrwidth_mass_meds_inf45; constrwidth_weight_sup = constrwidth_mass_meds_sup45
        if weightin.split('_')[1] == "mass2": constr_weight = constr_mass2_meds45; constrwidth_weight_inf = constrwidth_mass2_meds_inf45; constrwidth_weight_sup = constrwidth_mass2_meds_sup45
        if weightin.split('_')[1] == "mass3": constr_weight = constr_mass3_meds45; constrwidth_weight_inf = constrwidth_mass3_meds_inf45; constrwidth_weight_sup = constrwidth_mass3_meds_sup45
        if weightin.split('_')[1] == "oneoverr": constr_weight = constr_oneoverr_meds45; constrwidth_weight_inf = constrwidth_oneoverr_meds_inf45; constrwidth_weight_sup = constrwidth_oneoverr_meds_sup45
        if weightin.split('_')[1] == "zoverr": constr_weight = constr_zoverr_meds45; constrwidth_weight_inf = constrwidth_zoverr_meds_inf45; constrwidth_weight_sup = constrwidth_zoverr_meds_sup45
        if weightin.split('_')[1] == "massoverr": constr_weight = constr_massoverr_meds45; constrwidth_weight_inf = constrwidth_massoverr_meds_inf45; constrwidth_weight_sup = constrwidth_massoverr_meds_sup45
        if weightin.split('_')[1] == "mass2overr": constr_weight = constr_mass2overr_meds45; constrwidth_weight_inf = constrwidth_mass2overr_meds_inf45; constrwidth_weight_sup = constrwidth_mass2overr_meds_sup45
        if weightin.split('_')[1] == "mass3overr": constr_weight = constr_mass3overr_meds45; constrwidth_weight_inf = constrwidth_mass3overr_meds_inf45; constrwidth_weight_sup = constrwidth_mass3overr_meds_sup45
        if weightin.split('_')[1] == "mass2rms": constr_weight = constr_mass2rms_meds45; constrwidth_weight_inf = constrwidth_mass2rms_meds_inf45; constrwidth_weight_sup = constrwidth_mass2rms_meds_sup45
        if weightin.split('_')[1] == "mass3rms": constr_weight = constr_mass3rms_meds45; constrwidth_weight_inf = constrwidth_mass3rms_meds_inf45; constrwidth_weight_sup = constrwidth_mass3rms_meds_sup45
        if weightin.split('_')[1] == "mass2overrrms": constr_weight = constr_mass2overrrms_meds45; constrwidth_weight_inf = constrwidth_mass2overrrms_meds_inf45; constrwidth_weight_sup = constrwidth_mass2overrrms_meds_sup45
        if weightin.split('_')[1] == "mass3overrrms": constr_weight = constr_mass3overrrms_meds45; constrwidth_weight_inf = constrwidth_mass3overrrms_meds_inf45; constrwidth_weight_sup = constrwidth_mass3overrrms_meds_sup45
        if weightin.split('_')[1] == "flexion": constr_weight = constr_flexion_meds45; constrwidth_weight_inf = constrwidth_flexion_meds_inf45; constrwidth_weight_sup = constrwidth_flexion_meds_sup45
        if weightin.split('_')[1] == "tidal": constr_weight = constr_tidal_meds45; constrwidth_weight_inf = constrwidth_tidal_meds_inf45; constrwidth_weight_sup = constrwidth_tidal_meds_sup45
        if weightin.split('_')[1] == "SIS": constr_weight = constr_SIS_meds45; constrwidth_weight_inf = constrwidth_SIS_meds_inf45; constrwidth_weight_sup = constrwidth_SIS_meds_sup45
        if weightin.split('_')[1] == "SIShalo": constr_weight = constr_SIShalo_meds45; constrwidth_weight_inf = constrwidth_SIShalo_meds_inf45; constrwidth_weight_sup = constrwidth_SIShalo_meds_sup45
    if weightin.split('_')[0] == "120":
        if weightin.split('_')[1] == "gal": constr_weight = constr_gal_meds120; constrwidth_weight_inf = constrwidth_gal_meds_inf120; constrwidth_weight_sup = constrwidth_gal_meds_sup120
        if weightin.split('_')[1] == "z": constr_weight = constr_z_meds120; constrwidth_weight_inf = constrwidth_z_meds_inf120; constrwidth_weight_sup = constrwidth_z_meds_sup120
        if weightin.split('_')[1] == "mass": constr_weight = constr_mass_meds120; constrwidth_weight_inf = constrwidth_mass_meds_inf120; constrwidth_weight_sup = constrwidth_mass_meds_sup120
        if weightin.split('_')[1] == "mass2": constr_weight = constr_mass2_meds120; constrwidth_weight_inf = constrwidth_mass2_meds_inf120; constrwidth_weight_sup = constrwidth_mass2_meds_sup120
        if weightin.split('_')[1] == "mass3": constr_weight = constr_mass3_meds120; constrwidth_weight_inf = constrwidth_mass3_meds_inf120; constrwidth_weight_sup = constrwidth_mass3_meds_sup120
        if weightin.split('_')[1] == "oneoverr": constr_weight = constr_oneoverr_meds120; constrwidth_weight_inf = constrwidth_oneoverr_meds_inf120; constrwidth_weight_sup = constrwidth_oneoverr_meds_sup120
        if weightin.split('_')[1] == "zoverr": constr_weight = constr_zoverr_meds120; constrwidth_weight_inf = constrwidth_zoverr_meds_inf120; constrwidth_weight_sup = constrwidth_zoverr_meds_sup120
        if weightin.split('_')[1] == "massoverr": constr_weight = constr_massoverr_meds120; constrwidth_weight_inf = constrwidth_massoverr_meds_inf120; constrwidth_weight_sup = constrwidth_massoverr_meds_sup120
        if weightin.split('_')[1] == "mass2overr": constr_weight = constr_mass2overr_meds120; constrwidth_weight_inf = constrwidth_mass2overr_meds_inf120; constrwidth_weight_sup = constrwidth_mass2overr_meds_sup120
        if weightin.split('_')[1] == "mass3overr": constr_weight = constr_mass3overr_meds120; constrwidth_weight_inf = constrwidth_mass3overr_meds_inf120; constrwidth_weight_sup = constrwidth_mass3overr_meds_sup120
        if weightin.split('_')[1] == "mass2rms": constr_weight = constr_mass2rms_meds120; constrwidth_weight_inf = constrwidth_mass2rms_meds_inf120; constrwidth_weight_sup = constrwidth_mass2rms_meds_sup120
        if weightin.split('_')[1] == "mass3rms": constr_weight = constr_mass3rms_meds120; constrwidth_weight_inf = constrwidth_mass3rms_meds_inf120; constrwidth_weight_sup = constrwidth_mass3rms_meds_sup120
        if weightin.split('_')[1] == "mass2overrrms": constr_weight = constr_mass2overrrms_meds120; constrwidth_weight_inf = constrwidth_mass2overrrms_meds_inf120; constrwidth_weight_sup = constrwidth_mass2overrrms_meds_sup120
        if weightin.split('_')[1] == "mass3overrrms": constr_weight = constr_mass3overrrms_meds120; constrwidth_weight_inf = constrwidth_mass3overrrms_meds_inf120; constrwidth_weight_sup = constrwidth_mass3overrrms_meds_sup120
        if weightin.split('_')[1] == "flexion": constr_weight = constr_flexion_meds120; constrwidth_weight_inf = constrwidth_flexion_meds_inf120; constrwidth_weight_sup = constrwidth_flexion_meds_sup120
        if weightin.split('_')[1] == "tidal": constr_weight = constr_tidal_meds120; constrwidth_weight_inf = constrwidth_tidal_meds_inf120; constrwidth_weight_sup = constrwidth_tidal_meds_sup120
        if weightin.split('_')[1] == "SIS": constr_weight = constr_SIS_meds120; constrwidth_weight_inf = constrwidth_SIS_meds_inf120; constrwidth_weight_sup = constrwidth_SIS_meds_sup120
        if weightin.split('_')[1] == "SIShalo": constr_weight = constr_SIShalo_meds120; constrwidth_weight_inf = constrwidth_SIShalo_meds_inf120; constrwidth_weight_sup = constrwidth_SIShalo_meds_sup120
    if weightin.split('_')[1] == "gamma": constr_weight = constr_gamma; constrwidth_weight_inf = constrwidth_gamma_inf; constrwidth_weight_sup = constrwidth_gamma_sup
    return constr_weight, constrwidth_weight_inf, constrwidth_weight_sup

if conjoined == 5: constr_weight5, constrwidth_weight5_inf, constrwidth_weight5_sup = declareweight(weightin5)
if (conjoined == 4) | (conjoined == 5): constr_weight4, constrwidth_weight4_inf, constrwidth_weight4_sup = declareweight(weightin4)
if (conjoined == 3) | (conjoined == 4) | (conjoined == 5): constr_weight3, constrwidth_weight3_inf, constrwidth_weight3_sup = declareweight(weightin3)
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4) | (conjoined == 5): constr_weight2, constrwidth_weight2_inf, constrwidth_weight2_sup = declareweight(weightin2)
if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4) | (conjoined == 5): constr_weight1, constrwidth_weight1_inf, constrwidth_weight1_sup = declareweight(weightin1)

print "Reading..."

if mode == "sum": str1 = "sum"
if mode == "meds": str1 = "med"
if conjoined == 5:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4,increment5)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4,increment5)
if conjoined == 4:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4)
if conjoined == 3:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,mode,increment1,increment2,increment3)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,mode,increment1,increment2,increment3)
if conjoined == 2:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_increments%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,mode,increment1,increment2)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_LOS_increments%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,mode,increment1,increment2)
if conjoined == 1:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_increments%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,mode,increment1)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_LOS_increments%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,mode,increment1)
            
def readconjoined1_ugriz(radius,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup):
    ''' Here I only read the columns of interest, without kappa, for ugriz, in order to find the medians of their values over the whole MS.'''
    med1 = np.zeros(8)
    for i in range(8):
      if weightin1.split('_')[1] != "gamma":
          weight1_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index], unpack=True)
          if i == 0:
              weight1 = weight1_
          else:
              weight1 = np.append(weight1,weight1_)
      else:
          weight1_1_,weight1_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=[2,3], unpack=True)
          if i == 0:
              weight1_1 = weight1_1_
              weight1_2 = weight1_2_
          else:
              weight1_1 = np.append(weight1_1,weight1_1_)
              weight1_2 = np.append(weight1_2,weight1_2_)
    if weightin1.split('_')[1] != "gamma":
      med1 = np.median(weight1)
    else:
      med1 = np.median(np.sqrt(weight1_1**2 + weight1_2**2))
    print j
    med_weight1 = np.mean(med1) # throughout the code I use med_weight1 when computing intervals, following Green et al. For this, weight1 should always refer to simple galaxy number counts
    if weightin1.split('_')[1] == "gamma":
        constr_weight1 = constr_weight1 / med_weight1 # for gamma, measured shear divided by the median value of shear in MS; this turns it into an overdensity, like the other weights, so that it is meaningful to multiply by med_weight1
        constrwidth_weight1_inf = constrwidth_weight1_inf / med_weight1
        constrwidth_weight1_sup = constrwidth_weight1_sup / med_weight1
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))]) # absolute number, e.g. of galaxies within the lower width
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    return constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup

def readconjoined2_ugriz(radius,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup):
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,len,str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index), unpack=True)
          if i == 0:
              weight1 = weight1_
              weight2 = weight2_
          else:
              weight1 = np.append(weight1,weight1_)
              weight2 = np.append(weight2,weight2_)
      else:
          weight1_,weight2_1_,weight2_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index,1,2], unpack=True)
          if i == 0:
              weight1 = weight1_
              weight2_1 = weight2_1_
              weight2_2 = weight2_2_
          else:
              weight1 = np.append(weight1,weight1_)
              weight2_1 = np.append(weight2_1,weight2_1_)
              weight2_2 = np.append(weight2_2,weight2_2_)
    if weightin2.split('_')[1] != "gamma":
      med1 = np.median(weight1)
      med2 = np.median(weight2)
    else:
      med1 = np.median(weight1)
      med2 = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
    print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    if weightin2.split('_')[1] == "gamma":
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    return constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup

def readconjoined3_ugriz(radius,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup):
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index), unpack=True)
          if i == 0:
              weight1 = weight1_
              weight2 = weight2_
              weight3 = weight3_
          else:
              weight1 = np.append(weight1,weight1_)
              weight2 = np.append(weight2,weight2_)
              weight3 = np.append(weight3,weight3_)
      else:
          weight1_,weight2_1_,weight2_2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,1,2,weight3_index), unpack=True)
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
    if weightin2.split('_')[1] != "gamma":
      med1 = np.median(weight1)
      med2 = np.median(weight2)
      med3 = np.median(weight3)
    else:
      med1 = np.median(weight1)
      med2 = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
      med3 = np.median(weight3)
    print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    if weightin2.split('_')[1] == "gamma":
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    E_w3_inf = np.max([1, round(med_weight1 * (constr_weight3 - constrwidth_weight3_inf))])
    E_w3_sup = np.max([1, round(med_weight1 * (-constr_weight3 + constrwidth_weight3_sup))])
    return constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup

def readconjoined4_ugriz(radius,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup):
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    med4 = np.zeros(8)
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
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
          weight1_,weight2_1_,weight2_2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,1,2,weight3_index,weight4_index), unpack=True)
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
    if weightin2.split('_')[1] != "gamma":
      med1 = np.median(weight1)
      med2 = np.median(weight2)
      med3 = np.median(weight3)
      med4 = np.median(weight4)
    else:
      med1 = np.median(weight1)
      med2 = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
      med3 = np.median(weight3)
      med4 = np.median(weight4)
    print j
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    med_weight4 = np.mean(med4)
    if weightin2.split('_')[1] == "gamma":
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    E_w3_inf = np.max([1, round(med_weight1 * (constr_weight3 - constrwidth_weight3_inf))])
    E_w3_sup = np.max([1, round(med_weight1 * (-constr_weight3 + constrwidth_weight3_sup))])
    E_w4_inf = np.max([1, round(med_weight1 * (constr_weight4 - constrwidth_weight4_inf))])
    E_w4_sup = np.max([1, round(med_weight1 * (-constr_weight4 + constrwidth_weight4_sup))])
    return constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup


def readconjoined1_ugrizJHK(radius,constr_weight1):
    ''' Here I read ugrizJHK, converting weighted counts into overdensities, and recording the kappa values only for overdensities satisfying the constraint. I consider the full range of the constraint.'''
    for i in range(8):
      if weightin1.split('_')[1] != "gamma":
          id_,kappa_, weight1_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index), unpack=True)
          weight1_ = weight1_ / med_weight1
      else:
          id_,kappa_, gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,2,3), unpack=True)
          gamma1 = gamma1_
          gamma2 = gamma2_
          gamma = gamma1 # just so that the array has the correct shape
          gamma = np.sqrt(gamma1**2 + gamma2**2)
          weight1_ = gamma / med_weight1
      weight = np.copy(weight1_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      #del weight
      if (i == 0) and (j == 0):
          id = id_
          kappa = kappa_
          weight1 = weight1_
          ind1 = np.ones(np.shape(id_)) * j # this is to record the field name, and will be used together with id when matching different apertures
          ind2 = np.ones(np.shape(id_)) * i
      else:
          id = np.append(id,id_)
          ind1 = np.append(ind1,np.ones(np.shape(id_)) * j)
          ind2 = np.append(ind2,np.ones(np.shape(id_)) * i)
          kappa = np.append(kappa,kappa_)
          weight1 = np.append(weight1,weight1_)
    print j
    return id,ind1,ind2,kappa,weight1

def readconjoined2_ugrizJHK(radius,constr_weight1,constr_weight2):
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          id_,kappa_, weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,weight2_index), unpack=True)
          weight1_ = weight1_ / med_weight1
          weight2_ = weight2_ / med_weight2
      else:
          id_,kappa_,weight1_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,2,3), unpack=True)
          gamma1 = gamma1_
          gamma2 = gamma2_
          gamma = gamma1 # just so that the array has the correct shape
          gamma = np.sqrt(gamma1**2 + gamma2**2)
          weight1_ = weight1_ / med_weight1
          weight2_ = gamma / med_weight2
      weight = np.copy(weight1_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      #del weight
      weight = np.copy(weight2_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      #del weight
      if (i == 0) and (j == 0):
          id = id_
          kappa = kappa_
          weight1 = weight1_
          weight2 = weight2_
          ind1 = np.ones(np.shape(id_)) * j
          ind2 = np.ones(np.shape(id_)) * i
      else:
          id = np.append(id,id_)
          ind1 = np.append(ind1,np.ones(np.shape(id_)) * j)
          ind2 = np.append(ind2,np.ones(np.shape(id_)) * i)
          kappa = np.append(kappa,kappa_)
          weight1 = np.append(weight1,weight1_)
          weight2 = np.append(weight2,weight2_)
    print j
    return id,ind1,ind2,kappa,weight1,weight2

def readconjoined3_ugrizJHK(radius,constr_weight1,constr_weight2,constr_weight3):
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          id_,kappa_, weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,weight2_index,weight3_index), unpack=True)
          weight1_ = weight1_ / med_weight1
          weight2_ = weight2_ / med_weight2
          weight3_ = weight3_ / med_weight3
      else:
          id_,kappa_, weight1_,weight3_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,weight3_index,2,3), unpack=True)
          gamma1 = gamma1_
          gamma2 = gamma2_
          gamma = gamma1 # just so that the array has the correct shape
          gamma = np.sqrt(gamma1**2 + gamma2**2)
          weight1_ = weight1_ / med_weight1
          weight2_ = gamma / med_weight2
          weight3_ = weight3_ / med_weight3
      weight = np.copy(weight1_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      #del weight
      weight = np.copy(weight2_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      #del weight
      weight = np.copy(weight3_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      #del weight
      if (i == 0) and (j == 0):
          id = id_
          kappa = kappa_
          weight1 = weight1_
          weight2 = weight2_
          weight3 = weight3_
          ind1 = np.ones(np.shape(id_)) * j
          ind2 = np.ones(np.shape(id_)) * i
      else:
          id = np.append(id,id_)
          ind1 = np.append(ind1,np.ones(np.shape(id_)) * j)
          ind2 = np.append(ind2,np.ones(np.shape(id_)) * i)
          kappa = np.append(kappa,kappa_)
          weight1 = np.append(weight1,weight1_)
          weight2 = np.append(weight2,weight2_)
          weight3 = np.append(weight3,weight3_)
    print j
    return id,ind1,ind2,kappa,weight1,weight2,weight3

def readconjoined4_ugrizJHK(radius,constr_weight1,constr_weight2,constr_weight3,constr_weight4):
    for i in range(8):
      if weightin2.split('_')[1] != "gamma":
          id_,kappa_, weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
          weight1_ = weight1_ / med_weight1
          weight2_ = weight2_ / med_weight2
          weight3_ = weight3_ / med_weight3
          weight4_ = weight4_ / med_weight4
      else:
          id_,kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_0_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(i),mag,radius,innermask,zinf,zsup), usecols=(0,1,weight1_index,weight3_index,weight4_index,2,3), unpack=True)
          gamma1 = gamma1_
          gamma2 = gamma2_
          gamma = gamma1 # just so that the array has the correct shape
          gamma = np.sqrt(gamma1**2 + gamma2**2)
          weight1_ = weight1_ / med_weight1
          weight2_ = gamma / med_weight2
          weight3_ = weight3_ / med_weight3
          weight4_ = weight4_ / med_weight4
      weight = np.copy(weight1_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
      #del weight
      weight = np.copy(weight2_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
      #del weight
      weight = np.copy(weight3_)
      id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
      #del weight
      weight = np.copy(weight4_)
      kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
      weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
      weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
      weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
      weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
      #del weight
      if (i == 0) and (j == 0):
          id = id_
          kappa = kappa_
          weight1 = weight1_
          weight2 = weight2_
          weight3 = weight3_
          weight4 = weight4_
          ind1 = np.ones(np.shape(id_)) * j
          ind2 = np.ones(np.shape(id_)) * i
      else:
          id = np.append(id,id_)
          ind1 = np.append(ind1,np.ones(np.shape(id_)) * j)
          ind2 = np.append(ind2,np.ones(np.shape(id_)) * i)
          kappa = np.append(kappa,kappa_)
          weight1 = np.append(weight1,weight1_)
          weight2 = np.append(weight2,weight2_)
          weight3 = np.append(weight3,weight3_)
          weight4 = np.append(weight4,weight4_)
    print j
    return id,ind1,ind2,kappa,weight1,weight2,weight3,weight4

j = 0
if conjoined == 1:
    constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup) # the final three arguments are necessary in case of 'gamma', because they will be modified by the code
    id,ind1,ind2,kappa,weight1 = readconjoined1_ugrizJHK(weightin1.split('_')[0])
    #del id,ind1,ind2

if conjoined == 2:
    if weightin1.split('_')[0] == weightin2.split('_')[0]:
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        id,ind1,ind2,kappa,weight1,weight2 = readconjoined2_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2)
        #del id,ind1,ind2
    else:
        constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
        id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],constr_weight1)
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight2,E_w2_inf,E_w2_sup = readconjoined1_ugriz(weightin2.split('_')[0],constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_ = readconjoined1_ugrizJHK(weightin2.split('_')[0],constr_weight2)
        for i in range(8):
                id_rad1_ij = id_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                id_rad2_ij = id_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind1_rad1_ij = ind1_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind1_rad2_ij = ind1_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind2_rad1_ij = ind2_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind2_rad2_ij = ind2_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                kappa_rad1_ij = kappa_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                kappa_rad2_ij = kappa_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                weight1_ij = weight1_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                weight2_ij = weight2_[(ind1_rad2 == j) & (ind2_rad2 == i)]

                ind1_rad2_ij = ind1_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                ind2_rad2_ij = ind2_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                kappa_rad2_ij = kappa_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                weight2_ij = weight2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                id_rad2_ij = id_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
            
                ind1_rad1_ij = ind1_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                ind2_rad1_ij = ind2_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                kappa_rad1_ij = kappa_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                weight1_ij = weight1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                id_rad1_ij = id_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
            
                sort1 = np.argsort(id_rad1_ij)
                kappa_rad1_ij = kappa_rad1_ij[sort1]
                weight1_ij = weight1_ij[sort1]
                sort2 = np.argsort(id_rad2_ij)
                weight2_ij = weight2_ij[sort2]
                kappa_rad2_ij = kappa_rad2_ij[sort2]
                
                if kappa_rad1_ij != kappa_rad2_ij: print "error kappa"  # testing sanity
                if (i == 0) and (j == 0):
                    kappa = kappa_rad1_ij
                    weight1 = weight1_ij
                    weight2 = weight2_ij
                else:
                    kappa = np.append(kappa,kappa_rad1_ij)
                    weight1 = np.append(weight1,weight1_ij)
                    weight2 = np.append(weight2,weight2_ij)
        #del weight1_ij
        #del weight2_ij
        #del id_rad1_ij,id_rad1
        #del id_rad2_ij,id_rad2
        #del ind1_rad1_ij,ind1_rad1
        #del ind1_rad2_ij,ind1_rad2
        #del ind2_rad1_ij,ind2_rad1
        #del ind2_rad2_ij,ind2_rad2
        #del kappa_rad1_ij,kappa_rad1
        #del kappa_rad2_ij,kappa_rad2

if conjoined == 3:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
        id,ind1,ind2,kappa,weight1,weight2,weight3 = readconjoined3_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2,constr_weight3)
        #del id,ind1,ind2
    else:
        if weightin1.split('_')[0] == weightin2.split('_')[0]:
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
            id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight3,E_w3_inf,E_w3_sup = readconjoined1_ugriz(weightin3.split('_')[0],constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_ = readconjoined1_ugrizJHK(weightin3.split('_')[0],constr_weight3)
        if weightin2.split('_')[0] == weightin3.split('_')[0]:
            constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
            id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],constr_weight1)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined2_ugriz(weightin2.split('_')[0],constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_ = readconjoined2_ugrizJHK(weightin2.split('_')[0],constr_weight2,constr_weight3)
        for i in range(8):
                id_rad1_ij = id_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                id_rad2_ij = id_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind1_rad1_ij = ind1_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind1_rad2_ij = ind1_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind2_rad1_ij = ind2_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind2_rad2_ij = ind2_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                kappa_rad1_ij = kappa_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                kappa_rad2_ij = kappa_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                weight1_ij = weight1_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if weightin1.split('_')[0] == weightin2.split('_')[0]:
                    weight2_ij = weight2_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if weightin2.split('_')[0] == weightin3.split('_')[0]:
                    weight2_ij = weight2_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                weight3_ij = weight3_[(ind1_rad2 == j) & (ind2_rad2 == i)]
            
                ind1_rad2_ij = ind1_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                ind2_rad2_ij = ind2_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                kappa_rad2_ij = kappa_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if weightin2.split('_')[0] == weightin3.split('_')[0]:
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                weight3_ij = weight3_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                id_rad2_ij = id_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
            
                ind1_rad1_ij = ind1_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                ind2_rad1_ij = ind2_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                kappa_rad1_ij = kappa_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                weight1_ij = weight1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if weightin1.split('_')[0] == weightin2.split('_')[0]:
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                id_rad1_ij = id_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
            
                sort1 = np.argsort(id_rad1_ij)
                kappa_rad1_ij = kappa_rad1_ij[sort1]
                weight1_ij = weight1_ij[sort1]
                sort2 = np.argsort(id_rad2_ij)
                if weightin2.split('_')[0] == weightin3.split('_')[0]:
                    weight2_ij = weight2_ij[sort2]
                if weightin1.split('_')[0] == weightin2.split('_')[0]:
                    weight2_ij = weight2_ij[sort1]
                weight3_ij = weight3_ij[sort2]
                kappa_rad2_ij = kappa_rad2_ij[sort2]

                if kappa_rad1_ij != kappa_rad2_ij: print "error kappa"  # testing sanity
            
                if (i == 0) and (j == 0):
                    kappa = kappa_rad1_ij
                    weight1 = weight1_ij
                    weight2 = weight2_ij
                    weight3 = weight3_ij
                else:
                    kappa = np.append(kappa,kappa_rad1_ij)
                    weight1 = np.append(weight1,weight1_ij)
                    weight2 = np.append(weight2,weight2_ij)
                    weight3 = np.append(weight3,weight3_ij)
        #del weight1_ij,weight1_
        #del weight2_ij,weight2_
        #del weight3_ij,weight3_
        #del id_rad1_ij,id_rad1
        #del id_rad2_ij,id_rad2
        #del ind1_rad1_ij,ind1_rad1
        #del ind1_rad2_ij,ind1_rad2
        #del ind2_rad1_ij,ind2_rad1
        #del ind2_rad2_ij,ind2_rad2
        #del kappa_rad1_ij,kappa_rad1
        #del kappa_rad2_ij,kappa_rad2

if conjoined == 4:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined4_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
        id,ind1,ind2,kappa,weight1,weight2,weight3,weight4 = readconjoined4_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2,constr_weight3,constr_weight4)
        #del id,ind1,ind2
    else:
        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2,constr_weight3)
            constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight4,E_w4_inf,E_w4_sup = readconjoined1_ugriz(weightin4.split('_')[0],constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_ = readconjoined1_ugrizJHK(weightin4.split('_')[0],constr_weight4)
        if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
            constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
            id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],constr_weight1)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,med_weight4,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined3_ugriz(weightin2.split('_')[0],constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_ = readconjoined3_ugrizJHK(weightin2.split('_')[0])
        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
            id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2)
            constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight3,med_weight4,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined2_ugriz(weightin3.split('_')[0],constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_ = readconjoined2_ugrizJHK(weightin3.split('_')[0],constr_weight3,constr_weight4)
        for i in range(8):
                id_rad1_ij = id_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                id_rad2_ij = id_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind1_rad1_ij = ind1_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind1_rad2_ij = ind1_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind2_rad1_ij = ind2_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind2_rad2_ij = ind2_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                kappa_rad1_ij = kappa_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                kappa_rad2_ij = kappa_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                weight1_ij = weight1_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                    weight2_ij = weight2_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight2_ij = weight2_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if ((weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight3_ij = weight3_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                    weight3_ij = weight3_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                weight4_ij = weight4_[(ind1_rad2 == j) & (ind2_rad2 == i)]
            
                ind1_rad2_ij = ind1_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                ind2_rad2_ij = ind2_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                kappa_rad2_ij = kappa_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if ((weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight3_ij = weight3_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                weight4_ij = weight4_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                id_rad2_ij = id_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
            
                ind1_rad1_ij = ind1_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                ind2_rad1_ij = ind2_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                kappa_rad1_ij = kappa_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                weight1_ij = weight1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                    weight3_ij = weight3_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                id_rad1_ij = id_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
            
                sort1 = np.argsort(id_rad1_ij)
                kappa_rad1_ij = kappa_rad1_ij[sort1]
                weight1_ij = weight1_ij[sort1]
                sort2 = np.argsort(id_rad2_ij)
                if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                    weight2_ij = weight2_ij[sort2]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight2_ij = weight2_ij[sort1]
                if ((weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0])):
                    weight3_ij = weight3_ij[sort2]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                    weight3_ij = weight3_ij[sort1]
                weight4_ij = weight4_ij[sort2]
                kappa_rad2_ij = kappa_rad2_ij[sort2]

                if kappa_rad1_ij != kappa_rad2_ij: print "error kappa"  # testing sanity
            
                if (i == 0) and (j == 0):
                    kappa = kappa_rad1_ij
                    weight1 = weight1_ij
                    weight2 = weight2_ij
                    weight3 = weight3_ij
                    weight4 = weight4_ij
                else:
                    kappa = np.append(kappa,kappa_rad1_ij)
                    weight1 = np.append(weight1,weight1_ij)
                    weight2 = np.append(weight2,weight2_ij)
                    weight3 = np.append(weight3,weight3_ij)
                    weight4 = np.append(weight4,weight4_ij)
        #del weight1_ij,weight1_
        #del weight2_ij,weight2_
        #del weight3_ij,weight3_
        #del weight4_ij,weight4_
        #del id_rad1_ij,id_rad1
        #del id_rad2_ij,id_rad2
        #del ind1_rad1_ij,ind1_rad1
        #del ind1_rad2_ij,ind1_rad2
        #del ind2_rad1_ij,ind2_rad1
        #del ind2_rad2_ij,ind2_rad2
        #del kappa_rad1_ij,kappa_rad1
        #del kappa_rad2_ij,kappa_rad2

if conjoined == 5:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        sys.exit('For a single given radius, a maximum of 4 conjoined constraints is allowed (5 given).')
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined4_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
        id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_,weight4_ = readconjoined4_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2,constr_weight3,constr_weight4)
        constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup,med_weight5,E_w5_inf,E_w5_sup = readconjoined1_ugriz(weightin5.split('_')[0],constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight5_ = readconjoined1_ugrizJHK(weightin5.split('_')[0],constr_weight5)
    if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
        id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],constr_weight1)
        constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,med_weight4,med_weight5,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined4_ugriz(weightin2.split('_')[0],constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_,weight5_ = readconjoined4_ugrizJHK(weightin2.split('_')[0],constr_weight2,constr_weight3,constr_weight4,constr_weight5)
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
        id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2,constr_weight3)
        constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup,med_weight4,med_weight5,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined2_ugriz(weightin4.split('_')[0],constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_,weight5_ = readconjoined2_ugrizJHK(weightin4.split('_')[0],constr_weight4,constr_weight5)
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],constr_weight1,constr_weight2)
        constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight3,med_weight4,med_weight5,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined3_ugriz(weightin3.split('_')[0],constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_,weight5_ = readconjoined3_ugrizJHK(weightin3.split('_')[0],constr_weight4,constr_weight5)
    
    for i in range(8):
                id_rad1_ij = id_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                id_rad2_ij = id_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind1_rad1_ij = ind1_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind1_rad2_ij = ind1_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                ind2_rad1_ij = ind2_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                ind2_rad2_ij = ind2_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                kappa_rad1_ij = kappa_rad1[(ind1_rad1 == j) & (ind2_rad1 == i)]
                kappa_rad2_ij = kappa_rad2[(ind1_rad2 == j) & (ind2_rad2 == i)]
                weight1_ij = weight1_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                    weight2_ij = weight2_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight2_ij = weight2_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight4_ij = weight4_[(ind1_rad2 == j) & (ind2_rad2 == i)]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
                    weight4_ij = weight4_[(ind1_rad1 == j) & (ind2_rad1 == i)]
                weight5_ij = weight5_[(ind1_rad2 == j) & (ind2_rad2 == i)]
            
                ind1_rad2_ij = ind1_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                ind2_rad2_ij = ind2_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                kappa_rad2_ij = kappa_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight4_ij = weight4_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                weight5_ij = weight5_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
                id_rad2_ij = id_rad2_ij[np.where(np.in1d(id_rad2_ij, id_rad1_ij))[0]]
            
                ind1_rad1_ij = ind1_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                ind2_rad1_ij = ind2_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                kappa_rad1_ij = kappa_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                weight1_ij = weight1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
                    weight4_ij = weight4_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight2_ij = weight2_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
                id_rad1_ij = id_rad1_ij[np.where(np.in1d(id_rad1_ij, id_rad2_ij))[0]]
            
                sort1 = np.argsort(id_rad1_ij)
                kappa_rad1_ij = kappa_rad1_ij[sort1]
                weight1_ij = weight1_ij[sort1]
                sort2 = np.argsort(id_rad2_ij)
                if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                    weight2_ij = weight2_ij[sort2]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight2_ij = weight2_ij[sort1]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_ij[sort2]
                if ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight3_ij = weight3_ij[sort1]
                if ((weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])) | ((weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0])):
                    weight4_ij = weight4_ij[sort2]
                if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
                    weight4_ij = weight4_ij[sort1]
                weight5_ij = weight5_ij[sort2]
                kappa_rad2_ij = kappa_rad2_ij[sort2]

                if kappa_rad1_ij != kappa_rad2_ij: print "error kappa"  # testing sanity
            
                if (i == 0) and (j == 0):
                    kappa = kappa_rad1_ij
                    weight1 = weight1_ij
                    weight2 = weight2_ij
                    weight3 = weight3_ij
                    weight4 = weight4_ij
                    weight5 = weight5_ij
                else:
                    kappa = np.append(kappa,kappa_rad1_ij)
                    weight1 = np.append(weight1,weight1_ij)
                    weight2 = np.append(weight2,weight2_ij)
                    weight3 = np.append(weight3,weight3_ij)
                    weight4 = np.append(weight4,weight4_ij)
                    weight5 = np.append(weight5,weight5_ij)
    #del weight1_ij,weight1_
    #del weight2_ij,weight2_
    #del weight3_ij,weight3_
    #del weight4_ij,weight4_
    #del weight5_ij,weight5_
    #del id_rad1_ij,id_rad1
    #del id_rad2_ij,id_rad2
    #del ind1_rad1_ij,ind1_rad1
    #del ind1_rad2_ij,ind1_rad2
    #del ind2_rad1_ij,ind2_rad1
    #del ind2_rad2_ij,ind2_rad2
    #del kappa_rad1_ij,kappa_rad1
    #del kappa_rad2_ij,kappa_rad2

print(" Read in %s seconds" % (time.time() - start_time))

gauss = sp.stats.norm(0, 1)
start1 = time.time()
LOS = 0

if conjoined == 5:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1): # use as specific value
        for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2):
            for E3 in np.arange(-limsigma * E_w3_inf, limsigma * E_w3_sup + 1, increment3):
                for E4 in np.arange(-limsigma * E_w4_inf, limsigma * E_w4_sup + 1, increment4):
                    for E5 in np.arange(-limsigma * E_w5_inf, limsigma * E_w5_sup + 1, increment5):
                        print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") ", "E3 = ", E3, "in (", -limsigma * E_w3_inf, ",", limsigma * E_w3_sup, ") ", "E4 = ", E4, "in (", -limsigma * E_w4_inf, ",", limsigma * E_w4_sup, "E5 = ", E5, "in (", -limsigma * E_w5_inf, ",", limsigma * E_w5_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight1 >= round(constr_weight4 * med_weight1) + E4 - increment4/2.0) & (weight4 * med_weight1 < round(constr_weight4 * med_weight1) + E4 + increment4/2.0) & (weight5 * med_weight1 >= round(constr_weight5 * med_weight1) + E5 - increment5/2.0) & (weight5 * med_weight1 < round(constr_weight5 * med_weight1) + E5 + increment5/2.0)] # this is equation 3 in Greene et al.
                        if data.size > 0:
                            if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf)
                            else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                            if E2 < 0: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_inf)
                            else: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_sup)
                            if E3 < 0: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_inf)
                            else: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_sup)
                            if E4 < 0: gauss_factorE4 = gauss.pdf(float(E4)/E_w4_inf)
                            else: gauss_factorE4 = gauss.pdf(float(E4)/E_w4_sup)
                            if E5 < 0: gauss_factorE5 = gauss.pdf(float(E5)/E_w5_inf)
                            else: gauss_factorE5 = gauss.pdf(float(E5)/E_w5_sup)
                            kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 * gauss_factorE3 * gauss_factorE4 * gauss_factorE5 / data.shape[0]
                            if LOS == 0:
                                unbiased_kappa_constrained = kappa_constrained
                            else:
                                unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained
                        LOS = LOS + data.size
            
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
    if weightin1.split('_')[1] == "gamma":
        E_w1_inf = E_w1_inf45
        E_w1_sup = E_w1_inf45
        med_weight1 = med_weight1_45
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
    else:
        E_w1_inf = E_w1_inf45
        E_w1_sup = E_w1_inf45
        E_w2_inf = E_w1_inf120
        E_w2_sup = E_w1_inf120
        med_weight1 = med_weight1_45
        med_weight2 = med_weight1_120
        constr_weight1 = constr_weight1_45
        constr_weight2 = constr_weight1_120
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

np.savetxt(output,unbiased_kappa_constrained,fmt='%s', delimiter='\t',newline='\n')
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
