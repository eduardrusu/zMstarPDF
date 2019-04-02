# CE Rusu Mar 9 2019
# only use for HE0435 and for a single aperture radius
# uses the old MS files from 2016, converted to FITS
# Run as python /lfs08/rusucs/code/inferkappa_unbiasedwithshear45and120FITSioHE0435_only3conjoinedorsingleshear.py nohandpicked powerlawG1 5 24 45_gal 45_gamma 45_oneoverr
# shear can be used as a single constraint

import sys
import os
from os import system
import scipy as sp
from scipy.stats import norm
import numpy as np
import time
import fitsio # https://github.com/esheldon/fitsio

start_time=time.time()
handpickedstr = str(sys.argv[1])
other = str(sys.argv[2]) # refers to an optional suffix for the shear constraint
innermask = str(sys.argv[3])
mag = str(sys.argv[4])
conjoined = len(sys.argv) - 5 # total number of arguments including code name, minus the number of ones that are not weights

if innermask == "5": only8 = False
if innermask == "12": only8 = True # in this case run only 8/64 MS fields
shearwithoutprior = True # if True, do not divide by N_LOS on the shear constraint

if conjoined == 1:
    weightin1 = str(sys.argv[5])
if conjoined == 3:
    weightin1 = str(sys.argv[5])
    weightin2 = str(sys.argv[6])
    weightin3 = str(sys.argv[7])

print "conjoined:", conjoined
lens = "HE0435"
root = "/lfs08/rusucs/HE0435/MSwghtratios/"
rootcode = "/lfs08/rusucs/code/"
rootout = "/lfs08/rusucs/HE0435/MSkapparesults/"
#rootout = "/Volumes/LaCieSubaru/kapparesults/"
#rootout = "/mnt/scratch/rusucs/%s/MSkapparesults/" % lens
if innermask == '5': weightsfile = np.loadtxt(rootcode+'weightedcounts_%s_meds_%s_%sinner_%s.cat' %(lens,mag,innermask,handpickedstr),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
if innermask == '12': weightsfile = np.loadtxt(rootcode+'weightedcounts_%s_meds_%s_%sinner_%s.cat' %(lens,mag,innermask,handpickedstr),usecols=[1,2,3],unpack=True)
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

mediangal = 35 # I need this to compute the E interval for shear

increment1 = 2 # refers to the E interval from Greene et al. 2014
increment2 = 4
increment3 = 2

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
        constr_gamma = 0.041
        constrwidth_gamma_inf = 0.023
        constrwidth_gamma_sup = 0.059
    if other == 'powerlaw5pert' and innermask == '12':
        constr_gamma = 0.056
        constrwidth_gamma_inf = 0.030
        constrwidth_gamma_sup = 0.082
    if other == 'compositeG1' and innermask == '5':
        constr_gamma = 0.026
        constrwidth_gamma_inf = 0.000
        constrwidth_gamma_sup = 0.052
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

if innermask == "5":
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

if conjoined == 3: constr_weight3, constrwidth_weight3_inf, constrwidth_weight3_sup = declareweight(weightin3)
if (conjoined == 2) | (conjoined == 3): constr_weight2, constrwidth_weight2_inf, constrwidth_weight2_sup = declareweight(weightin2)
if (conjoined == 1) | (conjoined == 2) | (conjoined == 3): constr_weight1, constrwidth_weight1_inf, constrwidth_weight1_sup = declareweight(weightin1)

print "Reading..."

if conjoined == 3:
    if ((type(weight1_index) != int) | (type(weight2_index) != int) | (type(weight3_index) != int)) & (shearwithoutprior == True): stringshearprior = '_shearwithoutprior'
    else: stringshearprior = ''
if conjoined == 1:
    if (type(weight1_index) != int) & (shearwithoutprior == True): stringshearprior = '_shearwithoutprior'
    else: stringshearprior = ''
if conjoined == 3:
    output = '%skappahist_HE0435_%sinnermask_nobeta%s_%s_%s_%s_%s_%s_increments%s%s%s%s.cat' % (rootout,innermask,handpickedstr,other,weightin1,weightin2,weightin3,mag,increment1,increment2,increment3,stringshearprior)
if conjoined == 1:
    output = '%skappahist_HE0435_%sinnermask_nobeta%s_%s_%s_%s_increments%s_%s_gamma%s.cat' % (rootout,innermask,handpickedstr,other,weightin1,mag,increment1,constr_gamma,stringshearprior)

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
    E_w1_inf = np.max([1, round(mediangal * (constr_weight1 - constrwidth_weight1_inf))]) # absolute number, e.g. of galaxies within the lower width; med_weight1=mediangal for N_gal
    E_w1_sup = np.max([1, round(mediangal * (-constr_weight1 + constrwidth_weight1_sup))])
    return constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup

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
        id_ = id_[(weight * mediangal >= round(constr_weight1 * 35) - limsigma * E_w1_inf - increment1/2.0) & (weight * mediangal < round(constr_weight1 * 35) + limsigma * E_w1_sup + increment1/2.0) ]
        kappa_ = kappa_[(weight * mediangal >= round(constr_weight1 * 35) - limsigma * E_w1_inf - increment1/2.0) & (weight * mediangal < round(constr_weight1 * 35) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
        weight1_ = weight1_[(weight * mediangal >= round(constr_weight1 * 35) - limsigma * E_w1_inf - increment1/2.0) & (weight * mediangal < round(constr_weight1 * 35) + limsigma * E_w1_sup + increment1/2.0) ]
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
    return id,ind1,ind2,kappa,weight1

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

if conjoined == 1:
    constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,med_weight1,E_w1_inf,E_w1_sup = readconjoined1_ugriz(weightin1.split('_')[0],weight1_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup) # the constr argument representing 'gamma' (if 'gamma' is used) needs to be used here, so that it can be updated by the code. In addition, all constr arguments must be used because, as I am using multiple radii, global arguments such as constr_weight1 for one radius will correspond to a different global argument for another radius
    id,ind1,ind2,kappa,weight1 = readconjoined1_ugrizJHK(weightin1.split('_')[0],weight1_index,constr_weight1,increment1,med_weight1,E_w1_inf,E_w1_sup)
    del id,ind1,ind2

if conjoined == 3:
    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
        constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup = readconjoined3_ugriz(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constrwidth_weight1_inf,constrwidth_weight1_sup,constr_weight2,constrwidth_weight2_inf,constrwidth_weight2_sup,constr_weight3,constrwidth_weight3_inf,constrwidth_weight3_sup)
        id,ind1,ind2,kappa,weight1,weight2,weight3 = readconjoined3_ugrizJHK(weightin1.split('_')[0],weight1_index,weight2_index,weight3_index,constr_weight1,constr_weight2,constr_weight3,increment1,increment2,increment3,med_weight1,med_weight2,med_weight3,E_w1_inf,E_w1_sup,E_w2_inf,E_w2_sup,E_w3_inf,E_w3_sup)
        del id,ind1,ind2
    else:
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

print(" Read in %s seconds" % (time.time() - start_time))

gauss = sp.stats.norm(0, 1)
start1 = time.time()
LOS = 0
print np.shape(kappa)

if conjoined == 3:
    for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1):
        for E3 in np.arange(-limsigma * E_w3_inf, limsigma * E_w3_sup + 1, increment3):
            if ((type(weight1_index) != int) | (type(weight2_index) != int) | (type(weight3_index) != int)) & shearwithoutprior == True:
                datamarginalizedshear = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
            for E2 in np.arange(-limsigma * E_w2_inf, limsigma * E_w2_sup + 1, increment2): # I am integrating over the dimension with the shear first, so it must be nested deepest
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") ", "E2 = ", E2, "in (", -limsigma * E_w2_inf, ",", limsigma * E_w2_sup, ") ", "E3 = ", E3, "in (", -limsigma * E_w3_inf, ",", limsigma * E_w3_sup, ") "#, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    if (weightin1.split('_')[0] == weightin2.split('_')[0]) and (weightin2.split('_')[0] == weightin3.split('_')[0]):
                        data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0) & (weight2 * med_weight1 >= round(constr_weight2 * med_weight1) + E2 - increment2/2.0) & (weight2 * med_weight1 < round(constr_weight2 * med_weight1) + E2 + increment2/2.0) & (weight3 * med_weight1 >= round(constr_weight3 * med_weight1) + E3 - increment3/2.0) & (weight3 * med_weight1 < round(constr_weight3 * med_weight1) + E3 + increment3/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf)
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        if E2 < 0: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_inf)
                        else: gauss_factorE2 = gauss.pdf(float(E2)/E_w2_sup)
                        if E3 < 0: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_inf)
                        else: gauss_factorE3 = gauss.pdf(float(E3)/E_w3_sup)
                        if ((type(weight1_index) != int) | (type(weight2_index) != int) | (type(weight3_index) != int)) & shearwithoutprior == True:
                            kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 * gauss_factorE3 / datamarginalizedshear.shape[0]
                        else: kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 * gauss_factorE2 * gauss_factorE3 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained
                    LOS = LOS + data.size

if conjoined == 1:
        med_weight1 = mediangal # determined from N_gal
        for E1 in np.arange(-limsigma * E_w1_inf, limsigma * E_w1_sup + 1, increment1):
                    print "E1 = ", E1, "in (", -limsigma * E_w1_inf, ",", limsigma * E_w1_sup, ") " #, "gauss_weight4 = ", gauss.pdf(float(E4)/E_w4)
                    data = kappa[(weight1 * med_weight1 >= round(constr_weight1 * med_weight1) + E1 - increment1/2.0) & (weight1 * med_weight1 < round(constr_weight1 * med_weight1) + E1 + increment1/2.0)] # this is equation 3 in Greene et al.
                    if data.size > 0:
                        if E1 < 0: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_inf) # for asymmetric limits, implement a gaussian on each side
                        else: gauss_factorE1 = gauss.pdf(float(E1)/E_w1_sup)
                        if (type(weight1_index) != int) & (shearwithoutprior == True): kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1
                        else: kappa_constrained = np.histogram(data, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float) * gauss_factorE1 / data.shape[0]
                        if LOS == 0:
                            unbiased_kappa_constrained = kappa_constrained
                        else:
                            unbiased_kappa_constrained = unbiased_kappa_constrained + kappa_constrained # I tested that this addition works correctly
                    LOS = LOS + data.size

#head = 'LOS: %d' % np.array([LOS])
head = 'LOS: %d' % np.array([len(kappa)])
np.savetxt(output,unbiased_kappa_constrained,header=head,fmt='%s',delimiter='\t',newline='\n')
print(" time for computing kappa %s seconds" % (time.time() - start1))

if (conjoined == 1) | (conjoined == 2) | (conjoined == 3):
    print "increment1 = ", increment1
if (conjoined == 2) | (conjoined == 3):
    print "increment2 = ", increment2
if (conjoined == 3):
    print "increment3 = ", increment3

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
