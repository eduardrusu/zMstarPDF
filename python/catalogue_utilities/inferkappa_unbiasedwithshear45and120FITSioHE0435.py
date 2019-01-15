# CE Rusu July 21 2018
# NEED MAKE CHANGES WHEN RUNNING ALL BUT J1206 BECAUSE I WILL HAVE DIFFERENT INPUT FILES AND COLUMNS FOR 23 and 24
# Compared to inferkappa_unbiasedwithshear45and120.py, this code takes another argument ('empty' or != 'empty'); in case of 'empty', it only considers (after propoerly computing statistics using all LOS) only LOS without galaxies inside the inner mask. It requires that the input weight files have a final column which shows the number of galaxies inside the inner mask
# Run as python /lfs08/rusucs/code/inferkappa_unbiasedwithshear45and120FITSioHE0435.py HE0435 -1.0 -1.0 nohandpicked fiducial notempty notremovegroups 5 24 measured med 45_gal 45_gamma 45_oneoverr
# if shear is being used, it should be the second weight used in either radius
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
import fitsio # https://github.com/esheldon/fitsio

start_time=time.time()
only8 = False # in this case run only 8/64 MS fields

lens = str(sys.argv[1])
zinf = str(sys.argv[2])
zsup = str(sys.argv[3])
handpicked = str(sys.argv[4])
other = str(sys.argv[5]) # refers to an optional suffix for the shear constraint
empty = str(sys.argv[6])
removegroups = str(sys.argv[7])
innermask = str(sys.argv[8])
mag = str(sys.argv[9])
compmeas = str(sys.argv[10])
mode = str(sys.argv[11])
conjoined = len(sys.argv) - 12 # total number of arguments including code name, minus the number of ones that are not weights

handpickedstr = str(sys.argv[4])

if conjoined == 1:
    weightin1 = str(sys.argv[12])
if conjoined == 2:
    weightin1 = str(sys.argv[12])
    weightin2 = str(sys.argv[13])
if conjoined == 3:
    weightin1 = str(sys.argv[12])
    weightin2 = str(sys.argv[13])
    weightin3 = str(sys.argv[14])
if conjoined == 4:
    weightin1 = str(sys.argv[12])
    weightin2 = str(sys.argv[13])
    weightin3 = str(sys.argv[14])
    weightin4 = str(sys.argv[15])
if conjoined == 5:
    weightin1 = str(sys.argv[12])
    weightin2 = str(sys.argv[13])
    weightin3 = str(sys.argv[14])
    weightin4 = str(sys.argv[15])
    weightin5 = str(sys.argv[16])

print "conjoined:", conjoined
root = "/lfs08/rusucs/%s/MSwghtratios/" % lens
#root = "/mnt/scratch/rusucs/%s/MSwghtratios/" % lens
#root = "/Volumes/LaCieSubaru/MSweights/"
rootcode = "/lfs08/rusucs/code/"
#rootcode = "/Users/cerusu/Dropbox/Davis_work/code/J1206/"
#rootcode = "/mnt/scratch/rusucs/code/"
rootout = "/lfs08/rusucs/%s/MSkapparesults/" % lens
#rootout = "/Volumes/LaCieSubaru/kapparesults/"
#rootout = "/mnt/scratch/rusucs/%s/MSkapparesults/" % lens
weightsfile = np.loadtxt(rootcode+'weightedcounts_%s_%ss_%s_%sinner_%s.cat' %(lens,mode,mag,innermask,handpickedstr),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
if removegroups == 'removegroups': groupsfile = np.loadtxt(rootcode+'8_0_0groups.cat',usecols=[2,3,8],unpack=True)
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

increment1 = 4 # refers to the E interval from Greene et al. 2014
increment2 = 4
increment3 = 4
increment4 = 4
increment5 = 4

# these quantities are only for dealing with galaxy groups
degree = np.pi / 180
L_field = 4.0 * degree
N_pix_per_dim = 4096
L_pix = L_field / N_pix_per_dim
rad = degree / 3600

# define the shear constraints
# OLD VALUES
#if lens == "HE0435":
#    if other == 'fiducial' and innermask == '5':
#        constr_gamma = 0.030
#        constrwidth_gamma_inf = 0.027
#        constrwidth_gamma_sup = 0.033
#    if other == 'composite' and innermask == '5':
#        constr_gamma = 0.004
#        constrwidth_gamma_inf = 0.001
#        constrwidth_gamma_sup = 0.007
#    filters = "ugriJHK"
#    print 'shear: ',constr_gamma

if lens == "HE0435":
    if other == 'powerlawG1' and innermask == '5':
        constr_gamma = 0.047
        constrwidth_gamma_inf = 0.036
        constrwidth_gamma_sup = 0.058
    if other == 'powerlaw5pert' and innermask == '5':
        constr_gamma = 0.052
        constrwidth_gamma_inf = 0.041
        constrwidth_gamma_sup = 0.063
    if other == 'compositeG1' and innermask == '5':
        constr_gamma = 0.029
        constrwidth_gamma_inf = 0.004
        constrwidth_gamma_sup = 0.053
    filters = "ugriJHK"
    print 'shear: ',constr_gamma



# declare which weights to read
measured_index45 = 0 # specifies the column index in weightsfile
measured_index_inf45 = 1
measured_index_sup45 = 2
measured_index120 = 3
measured_index_inf120 = 4
measured_index_sup120 = 5

if mag == "24":
    def declareweight(weightin):
        if weightin.split('_')[1] == "gal": weight_index = 4
        if weightin.split('_')[1] == "z": weight_index = 6
        if weightin.split('_')[1] == "mass": weight_index = 8
        if weightin.split('_')[1] == "mass2": weight_index = 10
        if weightin.split('_')[1] == "mass3": weight_index = 12
        if weightin.split('_')[1] == "oneoverr": weight_index = 14
        if weightin.split('_')[1] == "zoverr": weight_index = 16
        if weightin.split('_')[1] == "massoverr": weight_index = 18
        if weightin.split('_')[1] == "mass2overr": weight_index = 20
        if weightin.split('_')[1] == "mass3overr": weight_index = 22
        if weightin.split('_')[1] == "mass2rms": weight_index = 24
        if weightin.split('_')[1] == "mass3rms": weight_index = 26
        if weightin.split('_')[1] == "mass2overrrms": weight_index = 28
        if weightin.split('_')[1] == "mass3overrrms": weight_index = 30
        if weightin.split('_')[1] == "flexion": weight_index = 32
        if weightin.split('_')[1] == "tidal": weight_index = 34
        if weightin.split('_')[1] == "SIS": weight_index = 36
        if weightin.split('_')[1] == "SIShalo": weight_index = 38
        if weightin.split('_')[1] == "gamma": weight_index = None
        return weight_index
if mag == "23":
    def declareweight(weightin):
        if weightin.split('_')[1] == "gal": weight_index = 5
        if weightin.split('_')[1] == "z": weight_index = 7
        if weightin.split('_')[1] == "mass": weight_index = 9
        if weightin.split('_')[1] == "mass2": weight_index = 11
        if weightin.split('_')[1] == "mass3": weight_index = 13
        if weightin.split('_')[1] == "oneoverr": weight_index = 15
        if weightin.split('_')[1] == "zoverr": weight_index = 17
        if weightin.split('_')[1] == "massoverr": weight_index = 19
        if weightin.split('_')[1] == "mass2overr": weight_index = 21
        if weightin.split('_')[1] == "mass3overr": weight_index = 23
        if weightin.split('_')[1] == "mass2rms": weight_index = 25
        if weightin.split('_')[1] == "mass3rms": weight_index = 27
        if weightin.split('_')[1] == "mass2overrrms": weight_index = 29
        if weightin.split('_')[1] == "mass3overrrms": weight_index = 31
        if weightin.split('_')[1] == "flexion": weight_index = 33
        if weightin.split('_')[1] == "tidal": weight_index = 35
        if weightin.split('_')[1] == "SIS": weight_index = 37
        if weightin.split('_')[1] == "SIShalo": weight_index = 39
        if weightin.split('_')[1] == "gamma": weight_index = None
        return weight_index

weight1_index = declareweight(weightin1)
if conjoined >= 2:
    weight2_index = declareweight(weightin2)
    if conjoined >= 3:
        weight3_index = declareweight(weightin3)
        if conjoined >= 4:
            weight4_index = declareweight(weightin4)
            if conjoined == 5:
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

if empty == 'empty': emptystr = '_emptymsk'
else: emptystr = ''
if removegroups == 'removegroups': groupsstr = '_removegroups'
else: groupsstr = ''
if conjoined == 5:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s_%s%s%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,weightin5,mag,mode,increment1,increment2,increment3,increment4,increment5,emptystr,groupsstr)
if conjoined == 4:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s%s%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,mode,increment1,increment2,increment3,increment4,emptystr,groupsstr)
if conjoined == 3:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s%s%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,mode,increment1,increment2,increment3,emptystr,groupsstr)
if conjoined == 2:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_increments%s_%s%s%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,mode,increment1,increment2,emptystr,groupsstr)
if conjoined == 1:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_increments%s%s%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,mode,increment1,emptystr,groupsstr)

def readfile(file,usecols):
    f = fitsio.FITS(file)
    print f # I need to print it, or f.hdu_list will not read
    ext = len(f.hdu_list)
    for i in range(ext - 1):
        if i == 0:
            data = fitsio.read(file, columns=usecols, ext=i+1)
        else:
            data = np.r_[data,fitsio.read(file, columns=usecols, ext=i+1)]
    # for speed, fitsio always returns columns and rows in order, so for instance in [1,2,3] even when usecols=[2,3,1]
    sort = np.argsort(np.argsort(usecols))
    if len(usecols) == 1:
        return data[data.dtype.names[0]]
    if len(usecols) == 2:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]]
    if len(usecols) == 3:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]]
    if len(usecols) == 4:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]],data[data.dtype.names[sort[3]]]
    if len(usecols) == 5:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]],data[data.dtype.names[sort[3]]],data[data.dtype.names[sort[4]]]
    if len(usecols) == 6:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]],data[data.dtype.names[sort[3]]],data[data.dtype.names[sort[4]]],data[data.dtype.names[sort[5]]]
    if len(usecols) == 7:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]],data[data.dtype.names[sort[3]]],data[data.dtype.names[sort[4]]],data[data.dtype.names[sort[5]]],data[data.dtype.names[sort[6]]]
    if len(usecols) == 8:
        return data[data.dtype.names[sort[0]]],data[data.dtype.names[sort[1]]],data[data.dtype.names[sort[2]]],data[data.dtype.names[sort[3]]],data[data.dtype.names[sort[4]]],data[data.dtype.names[sort[5]]],data[data.dtype.names[sort[6]]],data[data.dtype.names[sort[7]]]

def readconjoined1_ugriz(radius,weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup):
    ''' Here I only read the columns of interest, without kappa, for ugriz, in order to find the medians of their values over the whole MS.'''
    if only8 == True: field = 1
    else: field = 8
    med1 = np.zeros(field)
    filters1 = "ugriz"
    for j in range(field):
      for i in range(8):
        if type(weight1_index) == int:
            weight1_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=[weight1_index])
            if i == 0:
                weight1 = weight1_
            else:
                weight1 = np.append(weight1,weight1_)
        else:
            weight1_1_,weight1_2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=[2,3])
            if i == 0:
                weight1_1 = weight1_1_
                weight1_2 = weight1_2_
            else:
                weight1_1 = np.append(weight1_1,weight1_1_)
                weight1_2 = np.append(weight1_2,weight1_2_)
        #print j,i
      if type(weight1_index) == int:
        med1[j] = np.median(weight1)
      else:
        med1[j] = np.median(np.sqrt(weight1_1**2 + weight1_2**2))
    med_weight1 = np.mean(med1) # throughout the code I use med_weight1 when computing intervals, following Green et al. For this, weight1 should always refer to simple galaxy number counts
    if type(weight1_index) != int:
        constr_weight1 = constr_weight1 / med_weight1 # for gamma, measured shear divided by the median value of shear in MS; this turns it into an overdensity, like the other weights, so that it is meaningful to multiply later by the median number of galaxies
        constrwidth_weight1_inf = constrwidth_weight1_inf / med_weight1
        constrwidth_weight1_sup = constrwidth_weight1_sup / med_weight1
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))]) # absolute number, e.g. of galaxies within the lower width
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    return constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup

def readconjoined2_ugriz(radius,weight1_index,weight2_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup):
    if only8 == True: field = 1
    else: field = 8
    med1 = np.zeros(field)
    med2 = np.zeros(field)
    filters1 = "ugriz"
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            weight1_,weight2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(weight1_index,weight2_index))
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
        else:
            weight1_,weight2_1_,weight2_2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=[weight1_index,2,3])
            if i == 0:
                weight1 = weight1_
                weight2_1 = weight2_1_
                weight2_2 = weight2_2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2_1 = np.append(weight2_1,weight2_1_)
                weight2_2 = np.append(weight2_2,weight2_2_)
        #print j,i
      if type(weight2_index) == int:
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    if type(weight2_index) != int:
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    return constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup

def readconjoined3_ugriz(radius,weight1_index,weight2_index,weight3_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup):
    if only8 == True: field = 1
    else: field = 8
    med1 = np.zeros(field)
    med2 = np.zeros(field)
    med3 = np.zeros(field)
    filters1 = "ugriz"
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            weight1_,weight2_,weight3_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(weight1_index,weight2_index,weight3_index))
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
                weight3 = weight3_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
                weight3 = np.append(weight3,weight3_)
        else:
            weight1_,weight2_1_,weight2_2_,weight3_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(weight1_index,2,3,weight3_index))
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
        #print j,i
      if type(weight2_index) == int:
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
        med3[j] = np.median(weight3)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
        med3[j] = np.median(weight3)
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    if type(weight2_index) != int:
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

def readconjoined4_ugriz(radius,weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup):
    if only8 == True: field = 1
    else: field = 8
    med1 = np.zeros(field)
    med2 = np.zeros(field)
    med3 = np.zeros(field)
    med4 = np.zeros(field)
    filters1 = "ugriz"
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            weight1_,weight2_,weight3_,weight4_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(weight1_index,weight2_index,weight3_index,weight4_index))
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
            weight1_,weight2_1_,weight2_2_,weight3_,weight4_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(weight1_index,2,3,weight3_index,weight4_index))
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
        #print j,i
      if type(weight2_index) == int:
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
        med3[j] = np.median(weight3)
        med4[j] = np.median(weight4)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
        med3[j] = np.median(weight3)
        med4[j] = np.median(weight4)
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    med_weight4 = np.mean(med4)
    if type(weight2_index) != int:
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
    return constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup


def readconjoined1_ugrizJHK(radius,weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup):
    ''' Here I read ugrizJHK, converting weighted counts into overdensities, and recording the kappa values only for overdensities satisfying the constraint. I consider the full range of the constraint.'''
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        if type(weight1_index) == int:
            id_,kappa_, weight1_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index))
            weight1_ = weight1_ / med_weight1
        else:
            id_,kappa_, gamma1_,gamma2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = gamma / med_weight1
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        print np.shape(kappa_)
        del weight
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
        #print j,i
      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          for k in range(len(groupsfile[0])):
              print k,len(id)
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1


def readconjoined2_ugrizJHK(radius,weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            id_,kappa_, weight1_,weight2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,weight2_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
        else:
            id_,kappa_,weight1_,gamma1_,gamma2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight2_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        print np.shape(kappa_)
        del weight
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
        #print j,i
      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          #global x
          #global y
          x = np.array([])
          y = np.array([])
          for k in range(len(groupsfile[0])):
              print k,len(id)
              #x=np.append(x,k+1)
              #y=np.append(y,len(id))
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight2 = weight2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1,weight2

def readconjoined3_ugrizJHK(radius,weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            id_,kappa_, weight1_,weight2_,weight3_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,weight2_index,weight3_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
        else:
            id_,kappa_, weight1_,weight3_,gamma1_,gamma2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,weight3_index,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
            weight3_ = weight3_ / med_weight3
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight2_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight3_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        print np.shape(kappa_)
        del weight
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
        #print j,i
    return id,ind1,ind2,kappa,weight1,weight2,weight3

def readconjoined4_ugrizJHK(radius,weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        if type(weight2_index) == int:
            id_,kappa_, weight1_,weight2_,weight3_,weight4_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,weight2_index,weight3_index,weight4_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        else:
            id_,kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_ = readfile("%snobetaave3435NEWMEASUREDmedinject_%s_%s_GGL_los_8_%s_%s_%s.fits" % (root,filters1,lens,str(j),str(i),radius), usecols=(0,1,weight1_index,weight3_index,weight4_index,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight2_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight3_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        print np.shape(kappa_)
        del weight
        weight = np.copy(weight4_)
        id_ = id_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        print np.shape(kappa_)
        del weight
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
        #print j,i
    return id,ind1,ind2,kappa,weight1,weight2,weight3,weight4



def readconjoined1galinner_ugrizJHK(radius,weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        file = "%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup)
        f = fitsio.FITS(file)
        if 'galinner' not in f[1].get_colnames()[-1]:
            sys.exit('Only empty inner mask lines of sight requested, but the input file lacks the requested final column with the number of galaxies inside the mask.')
        else:
            if type(weight1_index) == int:
                id_,kappa_, weight1_,galinner_ = readfile(file, usecols=(0,1,weight1_index,len(f[1].get_colnames()) - 1))
                weight1_ = weight1_ / med_weight1
            else:
                id_,kappa_, gamma1_,gamma2_,galinner_ = np.loadtxt(file, usecols=(0,1,2,3,len(f[1].get_colnames()) - 1))
                gamma1 = gamma1_
                gamma2 = gamma2_
                gamma = gamma1 # just so that the array has the correct shape
                gamma = np.sqrt(gamma1**2 + gamma2**2)
                weight1_ = gamma / med_weight1
            weight = np.copy(weight1_)
            print np.shape(kappa_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            print np.shape(kappa_)
            del weight
            id_ = id_[galinner_ == 0]
            kappa_ = kappa_[galinner_ == 0]
            weight1_ = weight1_[galinner_ == 0]
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
            #print j,i
      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          for k in range(len(groupsfile[0])):
              print k,len(id)
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1

def readconjoined2galinner_ugrizJHK(radius,weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        file = "%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup)
        f = fitsio.FITS(file)
        if 'galinner' not in f[1].get_colnames()[-1]:
            sys.exit('Only empty inner mask lines of sight requested, but the input file lacks the requested final column with the number of galaxies inside the mask.')
        else:
            if type(weight2_index) == int:
                id_,kappa_, weight1_,weight2_,galinner_ = readfile(file, usecols=(0,1,weight1_index,weight2_index,len(f[1].get_colnames()) - 1))
                weight1_ = weight1_ / med_weight1
                weight2_ = weight2_ / med_weight2
            else:
                id_,kappa_,weight1_,gamma1_,gamma2_,galinner_ = readfile(file, usecols=(0,1,weight1_index,2,3,len(f[1].get_colnames()) - 1))
                gamma1 = gamma1_
                gamma2 = gamma2_
                gamma = gamma1 # just so that the array has the correct shape
                gamma = np.sqrt(gamma1**2 + gamma2**2)
                weight1_ = weight1_ / med_weight1
                weight2_ = gamma / med_weight2
            weight = np.copy(weight1_)
            print np.shape(kappa_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight2_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            print np.shape(kappa_)
            del weight
            id_ = id_[galinner_ == 0]
            kappa_ = kappa_[galinner_ == 0]
            weight1_ = weight1_[galinner_ == 0]
            weight2_ = weight2_[galinner_ == 0]
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
            #print j,i
      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          for k in range(len(groupsfile[0])):
              print k,len(id)
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight2 = weight2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1,weight2

def readconjoined3galinner_ugrizJHK(radius,weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        file = "%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup)
        f = fitsio.FITS(file)
        if 'galinner' not in f[1].get_colnames()[-1]:
            sys.exit('Only empty inner mask lines of sight requested, but the input file lacks the requested final column with the number of galaxies inside the mask.')
        else:
            if type(weight2_index) == int:
                id_,kappa_, weight1_,weight2_,weight3_,galinner_ = readfile(file, usecols=(0,1,weight1_index,weight2_index,weight3_index,len(f[1].get_colnames()) - 1))
                weight1_ = weight1_ / med_weight1
                weight2_ = weight2_ / med_weight2
                weight3_ = weight3_ / med_weight3
            else:
                id_,kappa_, weight1_,weight3_,gamma1_,gamma2_,galinner_ = readfile(file, usecols=(0,1,weight1_index,weight3_index,2,3,len(f[1].get_colnames()) - 1))
                gamma1 = gamma1_
                gamma2 = gamma2_
                gamma = gamma1 # just so that the array has the correct shape
                gamma = np.sqrt(gamma1**2 + gamma2**2)
                weight1_ = weight1_ / med_weight1
                weight2_ = gamma / med_weight2
                weight3_ = weight3_ / med_weight3
            weight = np.copy(weight1_)
            print np.shape(kappa_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight2_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight3_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            print np.shape(kappa_)
            del weight
            id_ = id_[galinner_ == 0]
            kappa_ = kappa_[galinner_ == 0]
            weight1_ = weight1_[galinner_ == 0]
            weight2_ = weight2_[galinner_ == 0]
            weight3_ = weight3_[galinner_ == 0]
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
      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          for k in range(len(groupsfile[0])):
              print k,len(id)
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight2 = weight2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight3 = weight3[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1,weight2,weight3

def readconjoined4galinner_ugrizJHK(radius,weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup):
    filters1 = filters
    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
      for i in range(8):
        file = "%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup)
        f = fitsio.FITS(file)
        if 'galinner' not in f[1].get_colnames()[-1]:
            sys.exit('Only empty inner mask lines of sight requested, but the input file lacks the requested final column with the number of galaxies inside the mask.')
        else:
            if type(weight2_index) == int:
                id_,kappa_, weight1_,weight2_,weight3_,weight4_,galinner_ = readfile(file, usecols=(0,1,weight1_index,weight2_index,weight3_index,weight4_index,len(f[1].get_colnames()) - 1))
                weight1_ = weight1_ / med_weight1
                weight2_ = weight2_ / med_weight2
                weight3_ = weight3_ / med_weight3
                weight4_ = weight4_ / med_weight4
            else:
                id_,kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_,galinner_ = readfile(file, usecols=(0,1,weight1_index,weight3_index,weight4_index,2,3,len(f[1].get_colnames()) - 1))
                gamma1 = gamma1_
                gamma2 = gamma2_
                gamma = gamma1 # just so that the array has the correct shape
                gamma = np.sqrt(gamma1**2 + gamma2**2)
                weight1_ = weight1_ / med_weight1
                weight2_ = gamma / med_weight2
                weight3_ = weight3_ / med_weight3
                weight4_ = weight4_ / med_weight4
            weight = np.copy(weight1_)
            print np.shape(kappa_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight2_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight3_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
            print np.shape(kappa_)
            del weight
            weight = np.copy(weight4_)
            id_ = id_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            galinner_ = galinner_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
            print np.shape(kappa_)
            del weight
            id_ = id_[galinner_ == 0]
            kappa_ = kappa_[galinner_ == 0]
            weight1_ = weight1_[galinner_ == 0]
            weight2_ = weight2_[galinner_ == 0]
            weight3_ = weight3_[galinner_ == 0]
            weight4_ = weight4_[galinner_ == 0]
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
            #print j,i

      if removegroups == 'removegroups':
          idposx = -0.5 * L_field  + (1 + (id / 4096) + 0.5) * L_pix
          idposy = -0.5 * L_field  + (1 + (id % 4096) + 0.5) * L_pix
          idposx = idposx[ind1 == 0]
          idposy = idposy[ind1 == 0]
          print "K: ", np.shape(id)
          for k in range(len(groupsfile[0])):
              print k,len(id)
              id = id[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind1 = ind1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              ind2 = ind2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              kappa = kappa[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight1 = weight1[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight2 = weight2[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight3 = weight3[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              weight4 = weight4[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx_ = idposx[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposy_ = idposy[(idposx - groupsfile[0][k])**2 + (idposy - groupsfile[1][k])**2 > (groupsfile[2][k]*rad)**2]
              idposx = idposx_
              idposy = idposy_
          print "K: ", np.shape(kappa)
    return id,ind1,ind2,kappa,weight1,weight2,weight3,weight4

if conjoined == 1:
    constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup) # the constr argument representing 'gamma' (if 'gamma' is used) needs to be used here, so that it can be updated by the code. In addition, all constr arguments must be used because, as I am using multiple radii, global arguments such as constr_weight1 for one radius will correspond to a different global argument for another radius
    if empty != 'empty': id,ind1,ind2,kappa,weight1 = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
    else: id,ind1,ind2,kappa,weight1 = readconjoined1galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
    del id,ind1,ind2

if conjoined == 2:
    if weightin1.split('_')[0] == weightin2.split('_')[0]:
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        if empty != 'empty': id,ind1,ind2,kappa,weight1,weight2 = readconjoined2_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
        else: id,ind1,ind2,kappa,weight1,weight2 = readconjoined2_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
        del id,ind1,ind2
    else:
        constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
        if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
        else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight2,E_w2_inf,E_w2_sup = readconjoined1_ugriz(weightin2.split('_')[0],weight2_index,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_ = readconjoined1_ugrizJHK(weightin2.split('_')[0],weight2_index,constr_weight2,increment2,med_weight2,E_w2_inf,E_w2_sup)
        else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_ = readconjoined1galinner_ugrizJHK(weightin2.split('_')[0],weight2_index,constr_weight2,increment2,med_weight2,E_w2_inf,E_w2_sup)
        if only8 == True: field = 1
        else: field = 8
        for j in range(field):
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

                diff = kappa_rad1_ij - kappa_rad2_ij
                if diff.any()!=0: print "error kappa"  # testing sanity

                if (i == 0) and (j == 0):
                    kappa = kappa_rad1_ij
                    weight1 = weight1_ij
                    weight2 = weight2_ij
                else:
                    kappa = np.append(kappa,kappa_rad1_ij)
                    weight1 = np.append(weight1,weight1_ij)
                    weight2 = np.append(weight2,weight2_ij)
        del sort1,sort2
        del weight1_ij
        del weight2_ij
        del id_rad1_ij,id_rad1
        del id_rad2_ij,id_rad2
        del ind1_rad1_ij,ind1_rad1
        del ind1_rad2_ij,ind1_rad2
        del ind2_rad1_ij,ind2_rad1
        del ind2_rad2_ij,ind2_rad2
        del kappa_rad1_ij,kappa_rad1
        del kappa_rad2_ij,kappa_rad2

if conjoined == 3:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
        if empty != 'empty': id,ind1,ind2,kappa,weight1,weight2,weight3 = readconjoined3_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        else: id,ind1,ind2,kappa,weight1,weight2,weight3 = readconjoined3galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        del id,ind1,ind2
    else:
        if weightin1.split('_')[0] == weightin2.split('_')[0]:
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
            if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
            else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight3,E_w3_inf,E_w3_sup = readconjoined1_ugriz(weightin3.split('_')[0],weight3_index,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_ = readconjoined1_ugrizJHK(weightin3.split('_')[0],weight3_index,constr_weight3,increment3,med_weight3,E_w3_inf,E_w3_sup)
            else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_ = readconjoined1galinner_ugrizJHK(weightin3.split('_')[0],weight3_index,constr_weight3,increment3,med_weight3,E_w3_inf,E_w3_sup)
        if weightin2.split('_')[0] == weightin3.split('_')[0]:
            constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
            if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
            else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined2_ugriz(weightin2.split('_')[0],weight2_index,weight3_index,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_ = readconjoined2_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,constr_weight2,constr_weight3,increment2,increment3,med_weight2,med_weight3,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
            else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_ = readconjoined2galinner_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,constr_weight2,constr_weight3,increment2,increment3,med_weight2,med_weight3,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        if only8 == True: field = 1
        else: field = 8
        for j in range(field):
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

                diff = kappa_rad1_ij - kappa_rad2_ij
                if diff.any()!=0: print "error kappa"  # testing sanity

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
        del sort1,sort2
        del weight1_ij,weight1_
        del weight2_ij,weight2_
        del weight3_ij,weight3_
        del id_rad1_ij,id_rad1
        del id_rad2_ij,id_rad2
        del ind1_rad1_ij,ind1_rad1
        del ind1_rad2_ij,ind1_rad2
        del ind2_rad1_ij,ind2_rad1
        del ind2_rad2_ij,ind2_rad2
        del kappa_rad1_ij,kappa_rad1
        del kappa_rad2_ij,kappa_rad2

if conjoined == 4:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined4_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
        if empty != 'empty': id,ind1,ind2,kappa,weight1,weight2,weight3,weight4 = readconjoined4_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        else: id,ind1,ind2,kappa,weight1,weight2,weight3,weight4 = readconjoined4galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        del id,ind1,ind2
    else:
        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
            if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
            else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
            constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight4,E_w4_inf,E_w4_sup = readconjoined1_ugriz(weightin4.split('_')[0],weight4_index,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_ = readconjoined1_ugrizJHK(weightin4.split('_')[0],weight4_index,constr_weight4,increment4,med_weight4,E_w4_inf,E_w4_sup)
            else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_ = readconjoined1galinner_ugrizJHK(weightin4.split('_')[0],weight4_index,constr_weight4,increment4,med_weight4,E_w4_inf,E_w4_sup)
        if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
            constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
            if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
            else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
            constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,med_weight4,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined3_ugriz(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_ = readconjoined3_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,constr_weight2,constr_weight3,constr_weight4,increment2,increment3,increment4,med_weight2,med_weight3,med_weight4,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
            else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_ = readconjoined3galinner_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,constr_weight2,constr_weight3,constr_weight4,increment2,increment3,increment4,med_weight2,med_weight3,med_weight4,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
            constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
            if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
            else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
            constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight3,med_weight4,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined2_ugriz(weightin3.split('_')[0],weight3_index,weight4_index,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
            if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_ = readconjoined2_ugrizJHK(weightin3.split('_')[0],weight3_index,weight4_index,constr_weight3,constr_weight4,increment3,increment4,med_weight3,med_weight4,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
            else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_ = readconjoined2galinner_ugrizJHK(weightin3.split('_')[0],weight3_index,weight4_index,constr_weight3,constr_weight4,increment3,increment4,med_weight3,med_weight4,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        if only8 == True: field = 1
        else: field = 8
        for j in range(field):
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

                diff = kappa_rad1_ij - kappa_rad2_ij
                if diff.any()!=0: print "error kappa"  # testing sanity

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
        del sort1,sort2
        del weight1_ij,weight1_
        del weight2_ij,weight2_
        del weight3_ij,weight3_
        del weight4_ij,weight4_
        del id_rad1_ij,id_rad1
        del id_rad2_ij,id_rad2
        del ind1_rad1_ij,ind1_rad1
        del ind1_rad2_ij,ind1_rad2
        del ind2_rad1_ij,ind2_rad1
        del ind2_rad2_ij,ind2_rad2
        del kappa_rad1_ij,kappa_rad1
        del kappa_rad2_ij,kappa_rad2

if conjoined == 5:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        sys.exit('For a single given radius, a maximum of 4 conjoined constraints is allowed (5 given).')
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup = readconjoined4_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup)
        if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_,weight4_ = readconjoined4_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_,weight4_ = readconjoined4galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,weight4_index,constr_weight1,constr_weight2,constr_weight3,constr_weight4,increment1,increment2,increment3,increment4,med_weight1,med_weight2,med_weight3,med_weight4,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup)
        constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup,med_weight5,E_w5_inf,E_w5_sup = readconjoined1_ugriz(weightin5.split('_')[0],weight5_index,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight5_ = readconjoined1_ugrizJHK(weightin5.split('_')[0],weight5_index,constr_weight5,increment5,med_weight5,E_w5_inf,E_w5_sup)
        else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight5_ = readconjoined1galinner_ugrizJHK(weightin5.split('_')[0],weight5_index,constr_weight5,increment5,med_weight5,E_w5_inf,E_w5_sup)
    if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup)
        if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight5,E_w5_inf,E_w5_sup)
        else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_ = readconjoined1galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight5,E_w5_inf,E_w5_sup)
        constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,med_weight2,med_weight3,med_weight4,med_weight5,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined4_ugriz(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,weight5_index,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_,weight5_ = readconjoined4_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,weight5_index,constr_weight2,constr_weight3,constr_weight4,constr_weight5,increment2,increment3,increment4,increment5,med_weight2,med_weight3,med_weight4,med_weight5,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)
        else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight2_,weight3_,weight4_,weight5_ = readconjoined4galinner_ugrizJHK(weightin2.split('_')[0],weight2_index,weight3_index,weight4_index,weight5_index,constr_weight2,constr_weight3,constr_weight4,constr_weight5,increment2,increment3,increment4,increment5,med_weight2,med_weight3,med_weight4,med_weight5,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
        if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_,weight3_ = readconjoined3galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup,med_weight4,med_weight5,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined2_ugriz(weightin4.split('_')[0],weight4_index,weight5_index,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_,weight5_ = readconjoined2_ugrizJHK(weightin4.split('_')[0],weight4_index,weight5_index,constr_weight4,constr_weight5,increment4,increment5,med_weight4,med_weight5,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)
        else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight4_,weight5_ = readconjoined2galinner_ugrizJHK(weightin4.split('_')[0],weight4_index,weight5_index,constr_weight4,constr_weight5,increment4,increment5,med_weight4,med_weight5,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup = readconjoined2_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup)
        if empty != 'empty': id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
        else: id_rad1,ind1_rad1,ind2_rad1,kappa_rad1,weight1_,weight2_ = readconjoined2galinner_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,constr_weight1,constr_weight2,increment1,increment2,med_weight1,med_weight2,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup)
        constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,med_weight3,med_weight4,med_weight5,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup = readconjoined3_ugriz(weightin3.split('_')[0],weight3_index,weight4_index,weight5_index,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup,constr_weight4,constrwidth_weight4_inf,constrwidth_weight4_sup,constr_weight5,constrwidth_weight5_inf,constrwidth_weight5_sup)
        if empty != 'empty': id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_,weight5_ = readconjoined3_ugrizJHK(weightin3.split('_')[0],weight3_index,weight4_index,weight5_index,constr_weight3,constr_weight4,constr_weight5,increment3,increment4,increment5,med_weight3,med_weight4,med_weight5,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)
        else: id_rad2,ind1_rad2,ind2_rad2,kappa_rad2,weight3_,weight4_,weight5_ = readconjoined3galinner_ugrizJHK(weightin3.split('_')[0],weight3_index,weight4_index,weight5_index,constr_weight3,constr_weight4,constr_weight5,increment3,increment4,increment5,med_weight3,med_weight4,med_weight5,E_w3_inf,E_w3_sup,E_w4_inf,E_w4_sup,E_w5_inf,E_w5_sup)

    if only8 == True: field = 1
    else: field = 8
    for j in range(field):
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

                diff = kappa_rad1_ij - kappa_rad2_ij
                if diff.any()!=0: print "error kappa"  # testing sanity

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
    del sort1,sort2
    del weight1_ij,weight1_
    del weight2_ij,weight2_
    del weight3_ij,weight3_
    del weight4_ij,weight4_
    del weight5_ij,weight5_
    del id_rad1_ij,id_rad1
    del id_rad2_ij,id_rad2
    del ind1_rad1_ij,ind1_rad1
    del ind1_rad2_ij,ind1_rad2
    del ind2_rad1_ij,ind2_rad1
    del ind2_rad2_ij,ind2_rad2
    del kappa_rad1_ij,kappa_rad1
    del kappa_rad2_ij,kappa_rad2

print(" Read in %s seconds" % (time.time() - start_time))

gauss = sp.stats.norm(0, 1)
start1 = time.time()
LOS = 0
print np.shape(kappa)

if conjoined == 5:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1): # use as specific value
        for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2):
            for E3 in np.arange(-limsigma * E_w3_inf, limsigma * E_w3_sup + 1, increment3):
                for E4 in np.arange(-limsigma * E_w4_inf, limsigma * E_w4_sup + 1, increment4):
                    for E5 in np.arange(-limsigma * E_w5_inf, limsigma * E_w5_sup + 1, increment5):
                        print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") ", "E3 = ", E3, "in (", -limsigma * E_w3_inf, ",", limsigma * E_w3_sup, ") ", "E4 = ", E4, "in (", -limsigma * E_w4_inf, ",", limsigma * E_w4_sup, "E5 = ", E5, "in (", -limsigma * E_w5_inf, ",", limsigma * E_w5_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] != weightin5.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight1 >= round(constr_weight4 * med_weight1) + E4 - increment4/2.0) & (weight4 * med_weight1 < round(constr_weight4 * med_weight1) + E4 + increment4/2.0) & (weight5 * med_weight5 >= round(constr_weight5 * med_weight5) + E5 - increment5/2.0) & (weight5 * med_weight5 < round(constr_weight5 * med_weight5) + E5 + increment5/2.0)] # this is equation 3 in Greene et al.
                        if (weightin1.split('_')[0] != weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight2 >= round(constr_weight2 * med_weight2) + E2 - increment2/2.0) & (weight2 * med_weight2 < round(constr_weight2 * med_weight2) + E2 + increment2/2.0) & (weight3 * med_weight2 >= round(constr_weight3 * med_weight2) + E3 - increment3/2.0) & (weight3 * med_weight2 < round(constr_weight3 * med_weight2) + E3 + increment3/2.0) & (weight4 * med_weight2 >= round(constr_weight4 * med_weight2) + E4 - increment4/2.0) & (weight4 * med_weight2 < round(constr_weight4 * med_weight2) + E4 + increment4/2.0) & (weight5 * med_weight2 >= round(constr_weight5 * med_weight2) + E5 - increment5/2.0) & (weight5 * med_weight2 < round(constr_weight5 * med_weight2) + E5 + increment5/2.0)]
                        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] != weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight4 >= round(constr_weight4 * med_weight4) + E4 - increment4/2.0) & (weight4 * med_weight4 < round(constr_weight4 * med_weight4) + E4 + increment4/2.0) & (weight5 * med_weight4 >= round(constr_weight5 * med_weight4) + E5 - increment5/2.0) & (weight5 * med_weight4 < round(constr_weight5 * med_weight4) + E5 + increment5/2.0)]
                        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] != weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]) and (weightin4.split('_')[0] == weightin5.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight3 >= round(constr_weight3 * med_weight3) + E3 - increment3/2.0) & (weight3 * med_weight3 < round(constr_weight3 * med_weight3) + E3 + increment3/2.0) & (weight4 * med_weight3 >= round(constr_weight4 * med_weight3) + E4 - increment4/2.0) & (weight4 * med_weight3 < round(constr_weight4 * med_weight3) + E4 + increment4/2.0) & (weight5 * med_weight3 >= round(constr_weight5 * med_weight3) + E5 - increment5/2.0) & (weight5 * med_weight3 < round(constr_weight5 * med_weight3) + E5 + increment5/2.0)]
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
                    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight1 >= round(constr_weight4 * med_weight1) + E4 - increment4/2.0) & (weight4 * med_weight1 < round(constr_weight4 * med_weight1) + E4 + increment4/2.0)] # this is equation 3 in Greene et al.
                    else:
                        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0) & (weight4 * med_weight4 >= round(constr_weight4 * med_weight4) + E4 - increment4/2.0) & (weight4 * med_weight4 < round(constr_weight4 * med_weight4) + E4 + increment4/2.0)]
                        if (weightin2.split('_')[0] == weightin3.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight2 >= round(constr_weight2 * med_weight2) + E2 - increment2/2.0) & (weight2 * med_weight2 < round(constr_weight2 * med_weight2) + E2 + increment2/2.0) & (weight3 * med_weight2 >= round(constr_weight3 * med_weight2) + E3 - increment3/2.0) & (weight3 * med_weight2 < round(constr_weight3 * med_weight2) + E3 + increment3/2.0) & (weight4 * med_weight2 >= round(constr_weight4 * med_weight2) + E4 - increment4/2.0) & (weight4 * med_weight2 < round(constr_weight4 * med_weight2) + E4 + increment4/2.0)]
                        if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin3.split('_')[0] == weightin4.split('_')[0]):
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight3 >= round(constr_weight3 * med_weight3) + E3 - increment3/2.0) & (weight3 * med_weight3 < round(constr_weight3 * med_weight3) + E3 + increment3/2.0) & (weight4 * med_weight3 >= round(constr_weight4 * med_weight3) + E4 - increment4/2.0) & (weight4 * med_weight3 < round(constr_weight4 * med_weight3) + E4 + increment4/2.0)]
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
                    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
                    else:
                        if weightin1.split('_')[0] == weightin2.split('_')[0]:
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight3 >= round(constr_weight3 * med_weight3) + E3 - increment3/2.0) & (weight3 * med_weight3 < round(constr_weight3 * med_weight3) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
                        if weightin2.split('_')[0] == weightin3.split('_')[0]:
                            data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight2 >= round(constr_weight2 * med_weight2) + E2 - increment2/2.0) & (weight2 * med_weight2 < round(constr_weight2 * med_weight2) + E2 + increment2/2.0) & (weight3 * med_weight2 >= round(constr_weight3 * med_weight2) + E3 - increment3/2.0) & (weight3 * med_weight2 < round(constr_weight3 * med_weight2) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
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
                    if weightin1.split('_')[0] == weightin2.split('_')[0]:
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0)] # this is equation 3 in Greene et al.
                    else:
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight2 >= round(constr_weight2 * med_weight2) + E2 - increment2/2.0) & (weight2 * med_weight2 < round(constr_weight2 * med_weight2) + E2 + increment2/2.0)]
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

#head = 'LOS: %d' % np.array([LOS])
head = 'LOS: %d' % np.array([len(kappa)])
np.savetxt(output,unbiased_kappa_constrained,header=head,fmt='%s',delimiter='\t',newline='\n')
print(" time for computing kappa %s seconds" % (time.time() - start1))

if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4) | (conjoined == 5):
    print "increment1 = ", increment1
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4) | (conjoined == 5):
    print "increment2 = ", increment2
if (conjoined == 3) | (conjoined == 4) | (conjoined == 5):
    print "increment3 = ", increment3
if (conjoined == 4) | (conjoined == 5):
    print "increment4 = ", increment4
if conjoined == 5:
    print "increment5 = ", increment5

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
