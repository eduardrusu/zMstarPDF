# CE Rusu Feb 14 2018
# Run as python /lfs08/rusucs/code/inferkappa_unbiasedwithshear.py WFI2033 -1.0 -1.0 yes fiducial 5 23 meds gal gamma oneoverr mass
# the code currently works for maglim 23 (WFI2033)
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

print "conjoined:", conjoined
root = "/lfs08/rusucs/%s/MSwghtratios/" % lens
rootout = "/lfs08/rusucs/%s/MSkapparesults/" % lens
weightsfile = np.loadtxt(root+'weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(lens,type,inner,handpickedstr,zinf,zsup)),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

increment1 = 2 # refers to the E interval from Greene et al. 2014
increment2 = 2
increment3 = 2
increment4 = 2

# define the shear constraints
if lens == "WFI2033":
    if other == 'fiducial' and handpicked == 'no' and int(zsup) < 0 and innermask == '5':
        constr_gamma = 0.154
        constrwidth_gamma_inf = 0.139
        constrwidth_gamma_sup = 0.169
    if other == 'chameleon' and handpicked == 'no' and int(zsup) < 0 and innermask == '5':
        constr_gamma = 0.193
        constrwidth_gamma_inf = 0.178
        constrwidth_gamma_sup = 0.208
    if other == 'fiducial' and (handpicked == 'yes' or innermask == '15') and int(zsup) < 0:
        constr_gamma = 0.09
        constrwidth_gamma_inf = 0.075
        constrwidth_gamma_sup = 0.105
    if other == 'fiducial' and (handpicked == 'yes' or innermask == '15') and int(zsup) < 0:
        constr_gamma = 0.09
        constrwidth_gamma_inf = 0.075
        constrwidth_gamma_sup = 0.105
    filters = "ugrizJHK"
    print 'shear: ',constr_gamma

# declare which weights to read
if mag == "23" and radius == "45":
    measured_index45 = 0 # specifies the column index in weightsfile
    measured_index_inf45 = 1
    measured_index_sup45 = 2
    measured_index120 = 3
    measured_index_inf120 = 4
    measured_index_sup120 = 5

def declareweight(weightin):
    if weightin == "gal": weight_index = 4 # column index in the files created by kappamed_insertstarsnobeta.py
    if weightin == "z": weight_index = 5
    if weightin == "mass": weight_index = 6
    if weightin == "mass2": weight_index = 7
    if weightin == "mass3": weight_index = 8
    if weightin == "oneoverr": weight_index = 9
    if weightin == "zoverr": weight_index = 10
    if weightin == "massoverr": weight_index = 11
    if weightin == "mass2overr": weight_index = 12
    if weightin == "mass3overr": weight_index = 13
    if weightin == "mass2rms": weight_index = 14
    if weightin == "mass3rms": weight_index = 15
    if weightin == "mass2overrrms": weight_index = 16
    if weightin == "mass3overrrms": weight_index = 17
    if weightin == "flexion": weight_index = 18
    if weightin == "tidal": weight_index = 19
    if weightin == "SIS": weight_index = 20
    if weightin == "SIShalo": weight_index = 21
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

def declareweight45(weightin):
    if weightin == "gal": constr_weight = constr_gal_meds45; constrwidth_weight_inf = constrwidth_gal_meds_inf45; constrwidth_weight_sup = constrwidth_gal_meds_sup45
    if weightin == "z": constr_weight = constr_z_meds45; constrwidth_weight_inf = constrwidth_z_meds_inf45; constrwidth_weight_sup = constrwidth_z_meds_sup45
    if weightin == "mass": constr_weight = constr_mass_meds45; constrwidth_weight_inf = constrwidth_mass_meds_inf45; constrwidth_weight_sup = constrwidth_mass_meds_sup45
    if weightin == "mass2": constr_weight = constr_mass2_meds45; constrwidth_weight_inf = constrwidth_mass2_meds_inf45; constrwidth_weight_sup = constrwidth_mass2_meds_sup45
    if weightin == "mass3": constr_weight = constr_mass3_meds45; constrwidth_weight_inf = constrwidth_mass3_meds_inf45; constrwidth_weight_sup = constrwidth_mass3_meds_sup45
    if weightin == "oneoverr": constr_weight = constr_oneoverr_meds45; constrwidth_weight_inf = constrwidth_oneoverr_meds_inf45; constrwidth_weight_sup = constrwidth_oneoverr_meds_sup45
    if weightin == "zoverr": constr_weight = constr_zoverr_meds45; constrwidth_weight_inf = constrwidth_zoverr_meds_inf45; constrwidth_weight_sup = constrwidth_zoverr_meds_sup45
    if weightin == "massoverr": constr_weight = constr_massoverr_meds45; constrwidth_weight_inf = constrwidth_massoverr_meds_inf45; constrwidth_weight_sup = constrwidth_massoverr_meds_sup45
    if weightin == "mass2overr": constr_weight = constr_mass2overr_meds45; constrwidth_weight_inf = constrwidth_mass2overr_meds_inf45; constrwidth_weight_sup = constrwidth_mass2overr_meds_sup45
    if weightin == "mass3overr": constr_weight = constr_mass3overr_meds45; constrwidth_weight_inf = constrwidth_mass3overr_meds_inf45; constrwidth_weight_sup = constrwidth_mass3overr_meds_sup45
    if weightin == "mass2rms": constr_weight = constr_mass2rms_meds45; constrwidth_weight_inf = constrwidth_mass2rms_meds_inf45; constrwidth_weight_sup = constrwidth_mass2rms_meds_sup45
    if weightin == "mass3rms": constr_weight = constr_mass3rms_meds45; constrwidth_weight_inf = constrwidth_mass3rms_meds_inf45; constrwidth_weight_sup = constrwidth_mass3rms_meds_sup45
    if weightin == "mass2overrrms": constr_weight = constr_mass2overrrms_meds45; constrwidth_weight_inf = constrwidth_mass2overrrms_meds_inf45; constrwidth_weight_sup = constrwidth_mass2overrrms_meds_sup45
    if weightin == "mass3overrrms": constr_weight = constr_mass3overrrms_meds45; constrwidth_weight_inf = constrwidth_mass3overrrms_meds_inf45; constrwidth_weight_sup = constrwidth_mass3overrrms_meds_sup45
    if weightin == "flexion": constr_weight = constr_flexion_meds45; constrwidth_weight_inf = constrwidth_flexion_meds_inf45; constrwidth_weight_sup = constrwidth_flexion_meds_sup45
    if weightin == "tidal": constr_weight = constr_tidal_meds45; constrwidth_weight_inf = constrwidth_tidal_meds_inf45; constrwidth_weight_sup = constrwidth_tidal_meds_sup45
    if weightin == "SIS": constr_weight = constr_SIS_meds45; constrwidth_weight_inf = constrwidth_SIS_meds_inf45; constrwidth_weight_sup = constrwidth_SIS_meds_sup45
    if weightin == "SIShalo": constr_weight = constr_SIShalo_meds45; constrwidth_weight_inf = constrwidth_SIShalo_meds_inf45; constrwidth_weight_sup = constrwidth_SIShalo_meds_sup45
    if weightin == "gamma": constr_weight = constr_gamma45; constrwidth_weight_inf = constrwidth_gamma_inf45; constrwidth_weight_sup = constrwidth_gamma_sup45
    return constr_weight, constrwidth_weight_inf, constrwidth_weight_sup

def declareweight120(weightin):
    if weightin == "gal": constr_weight = constr_gal_meds120; constrwidth_weight_inf = constrwidth_gal_meds_inf120; constrwidth_weight_sup = constrwidth_gal_meds_sup120
    if weightin == "z": constr_weight = constr_z_meds120; constrwidth_weight_inf = constrwidth_z_meds_inf120; constrwidth_weight_sup = constrwidth_z_meds_sup120
    if weightin == "mass": constr_weight = constr_mass_meds120; constrwidth_weight_inf = constrwidth_mass_meds_inf120; constrwidth_weight_sup = constrwidth_mass_meds_sup120
    if weightin == "mass2": constr_weight = constr_mass2_meds120; constrwidth_weight_inf = constrwidth_mass2_meds_inf120; constrwidth_weight_sup = constrwidth_mass2_meds_sup120
    if weightin == "mass3": constr_weight = constr_mass3_meds120; constrwidth_weight_inf = constrwidth_mass3_meds_inf120; constrwidth_weight_sup = constrwidth_mass3_meds_sup120
    if weightin == "oneoverr": constr_weight = constr_oneoverr_meds120; constrwidth_weight_inf = constrwidth_oneoverr_meds_inf120; constrwidth_weight_sup = constrwidth_oneoverr_meds_sup120
    if weightin == "zoverr": constr_weight = constr_zoverr_meds120; constrwidth_weight_inf = constrwidth_zoverr_meds_inf120; constrwidth_weight_sup = constrwidth_zoverr_meds_sup120
    if weightin == "massoverr": constr_weight = constr_massoverr_meds120; constrwidth_weight_inf = constrwidth_massoverr_meds_inf120; constrwidth_weight_sup = constrwidth_massoverr_meds_sup120
    if weightin == "mass2overr": constr_weight = constr_mass2overr_meds120; constrwidth_weight_inf = constrwidth_mass2overr_meds_inf120; constrwidth_weight_sup = constrwidth_mass2overr_meds_sup120
    if weightin == "mass3overr": constr_weight = constr_mass3overr_meds120; constrwidth_weight_inf = constrwidth_mass3overr_meds_inf120; constrwidth_weight_sup = constrwidth_mass3overr_meds_sup120
    if weightin == "mass2rms": constr_weight = constr_mass2rms_meds120; constrwidth_weight_inf = constrwidth_mass2rms_meds_inf120; constrwidth_weight_sup = constrwidth_mass2rms_meds_sup120
    if weightin == "mass3rms": constr_weight = constr_mass3rms_meds120; constrwidth_weight_inf = constrwidth_mass3rms_meds_inf120; constrwidth_weight_sup = constrwidth_mass3rms_meds_sup120
    if weightin == "mass2overrrms": constr_weight = constr_mass2overrrms_meds120; constrwidth_weight_inf = constrwidth_mass2overrrms_meds_inf120; constrwidth_weight_sup = constrwidth_mass2overrrms_meds_sup120
    if weightin == "mass3overrrms": constr_weight = constr_mass3overrrms_meds120; constrwidth_weight_inf = constrwidth_mass3overrrms_meds_inf120; constrwidth_weight_sup = constrwidth_mass3overrrms_meds_sup120
    if weightin == "flexion": constr_weight = constr_flexion_meds120; constrwidth_weight_inf = constrwidth_flexion_meds_inf120; constrwidth_weight_sup = constrwidth_flexion_meds_sup120
    if weightin == "tidal": constr_weight = constr_tidal_meds120; constrwidth_weight_inf = constrwidth_tidal_meds_inf120; constrwidth_weight_sup = constrwidth_tidal_meds_sup120
    if weightin == "SIS": constr_weight = constr_SIS_meds120; constrwidth_weight_inf = constrwidth_SIS_meds_inf120; constrwidth_weight_sup = constrwidth_SIS_meds_sup120
    if weightin == "SIShalo": constr_weight = constr_SIShalo_meds120; constrwidth_weight_inf = constrwidth_SIShalo_meds_inf120; constrwidth_weight_sup = constrwidth_SIShalo_meds_sup120
    if weightin == "gamma": constr_weight = constr_gamma120; constrwidth_weight_inf = constrwidth_gamma_inf120; constrwidth_weight_sup = constrwidth_gamma_sup120
    return constr_weight, constrwidth_weight_inf, constrwidth_weight_sup

if conjoined == 4: constr_weight4_45, constrwidth_weight4_inf45, constrwidth_weight4_sup45 = declareweight45(weightin4)
if (conjoined == 3) | (conjoined == 4): constr_weight3_45, constrwidth_weight3_inf45, constrwidth_weight3_sup45 = declareweight45(weightin3)
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight2_45, constrwidth_weight2_inf45, constrwidth_weight2_sup45 = declareweight45(weightin2)
if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight1_45, constrwidth_weight1_inf45, constrwidth_weight1_sup45 = declareweight45(weightin1)

if conjoined == 4: constr_weight4_120, constrwidth_weight4_inf120, constrwidth_weight4_sup120 = declareweight120(weightin4)
if (conjoined == 3) | (conjoined == 4): constr_weight3_120, constrwidth_weight3_inf120, constrwidth_weight3_sup120 = declareweight120(weightin3)
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight2_120, constrwidth_weight2_inf120, constrwidth_weight2_sup120 = declareweight120(weightin2)
if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight1_120, constrwidth_weight1_inf120, constrwidth_weight1_sup120 = declareweight120(weightin1)

print "Reading..."

if mode == "sum": str1 = "sum"
if mode == "meds": str1 = "med"

if conjoined == 4:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_45120_%s_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_45120_%s_LOS_increments%s_%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4)
if conjoined == 3:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_45120_%s_increments%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,mode,increment1,increment2,increment3)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_45120_%s_LOS_increments%s_%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,mode,increment1,increment2,increment3)
if conjoined == 2:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_45120_%s_increments%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,mode,increment1,increment2)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_45120_%s_LOS_increments%s_%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,mode,increment1,increment2)
if conjoined == 1:
    output = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_45120_%s_increments%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,mode,increment1)
    outputLOS = '%skappahist_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_45120_%s_LOS_increments%s.cat' % (rootout,lens,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,mode,increment1)




















if conjoined == 1:
    ''' Here I only read the columns of interest, without kappa, for ugriz, in order to find the medians of their values over the whole MS.'''
    med1 = np.zeros(8)
    for j in range(8):
      for i in range(8):
        if weightin1 != "gamma":
            weight1_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index], unpack=True)
            if i == 0:
                weight1 = weight1_
            else:
                weight1 = np.append(weight1,weight1_)
        else:
            weight1_1_,weight1_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[2,3], unpack=True)
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
            kappa_, weight1_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index), unpack=True)
            weight1_ = weight1_ / med_weight1
        else:
            kappa_, gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,2,3), unpack=True)
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
            weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
        else:
            weight1_,weight2_1_,weight2_2_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index,1,2], unpack=True)
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
            kappa_, weight1_,weight2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
        else:
            kappa_, weight1_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,2,3), unpack=True)
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
            weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index), unpack=True)
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
                weight3 = weight3_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
                weight3 = np.append(weight3,weight3_)
        else:
            weight1_,weight2_1_,weight2_2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,1,2,weight3_index), unpack=True)
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
            kappa_, weight1_,weight2_,weight3_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index,weight3_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
        else:
            kappa_, weight1_,weight3_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight3_index,2,3), unpack=True)
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
            weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
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
            weight1_,weight2_1_,weight2_2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_ugriz_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,1,2,weight3_index,weight4_index), unpack=True)
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
            kappa_, weight1_,weight2_,weight3_,weight4_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index,weight3_index,weight4_index), unpack=True)
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        else:
            kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_ = np.loadtxt("%snobeta35measured%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.cat" % (root,str1,filters,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight3_index,weight4_index,2,3), unpack=True)
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
