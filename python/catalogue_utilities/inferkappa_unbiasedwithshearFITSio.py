# CE Rusu July 21 2018
# NEED MAKE CHANGES WHEN RUNNING ALL BUT J1206 BECAUSE I WILL HAVE DIFFERENT INPUT FILES AND COLUMNS FOR 23 and 24
# Run as python /lfs08/rusucs/code/inferkappa_unbiasedwithshearFITSio.py WFI2033 -1.0 -1.0 nohandpicked fiducial 5 45 23 measured med gal gamma oneoverr mass
# It does not accept mixed radii or selecting empty inner radii. When a single radius is used (not mixing different radii constraints) this code is faster than inferkappa_unbiasedwithshear45and120_23or24_allowsemptymskandJ1206 because it doesn't read the id column
# The input weight files have to be FITS files. In case of multiple data extensions, data is combined from every extension
# J1206 is considered separately because their the input weight files do not include columns which use mass
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
import fitsio # https://github.com/esheldon/fitsio

start_time=time.time()

lens = str(sys.argv[1])
zinf = str(sys.argv[2])
zsup = str(sys.argv[3])
handpicked = str(sys.argv[4])
other = str(sys.argv[5]) # refers to an optional suffix for the shear constraint
innermask = str(sys.argv[6])
radius = str(sys.argv[7])
mag = str(sys.argv[8])
compmeas = str(sys.argv[9])
mode = str(sys.argv[10])
conjoined = len(sys.argv) - 11 # total number of arguments including code name, minus the number of ones that are not weights

if handpicked == 'nohandpicked': handpickedstr = ''
else: handpickedstr = '_'+str(sys.argv[4])

if conjoined == 1:
    weightin1 = str(sys.argv[11])
if conjoined == 2:
    weightin1 = str(sys.argv[11])
    weightin2 = str(sys.argv[12])
if conjoined == 3:
    weightin1 = str(sys.argv[11])
    weightin2 = str(sys.argv[12])
    weightin3 = str(sys.argv[13])
if conjoined == 4:
    weightin1 = str(sys.argv[11])
    weightin2 = str(sys.argv[12])
    weightin3 = str(sys.argv[13])
    weightin4 = str(sys.argv[14])

print "conjoined:", conjoined
#root = "/lfs08/rusucs/%s/MSwghtratios/" % lens
root = "/mnt/scratch/rusucs/%s/MSwghtratios/" % lens
#root = "/Volumes/LaCieSubaru/MSweights/"
rootcode = "/mnt/scratch/rusucs/code/"
#rootcode = "/lfs08/rusucs/code/"
#rootcode = "/Users/cerusu/Dropbox/Davis_work/code/J1206/"
#rootout = "/lfs08/rusucs/%s/MSkapparesults/" % lens
#rootout = "/Volumes/LaCieSubaru/kapparesults/"
rootout = "/mnt/scratch/rusucs/%s/kapparesults/" % lens
#weightsfile = np.loadtxt(root+'weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(lens,mode,innermask,handpickedstr,zinf,zsup),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
weightsfile = np.loadtxt(rootcode+'weightedcounts_%s_%ss_%s_%sinner%s_zgap%s_%s.cat' %(lens,mode,mag,innermask,handpickedstr,zinf,zsup),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
#weightsfile = np.loadtxt('/Users/cerusu/Dropbox/Davis_work/code/%s/weightedcounts_%s_%s_%s_%sinner%s_zgap%s_%s.cat' %(lens,lens,mode,mag,innermask,handpickedstr,zinf,zsup),usecols=[1,2,3,4,5,6],unpack=True) # the file where I recorded the overdensities which I measured for the real lens
limsigma = 2 # sigma limits on either side of the assumed gaussians
bin_stat = 2000
min_kappa = -0.10
max_kappa = 1

increment1 = 4 # refers to the E interval from Greene et al. 2014
increment2 = 10
increment3 = 4
increment4 = 2

# define the shear constraints
if lens == "WFI2033":
    if other == 'fiducial' and handpicked == 'nohandpicked' and float(zsup) < 0 and innermask == '5':
        constr_gamma = 0.154
        constrwidth_gamma_inf = 0.139
        constrwidth_gamma_sup = 0.169
    if other == 'chameleon' and handpicked == 'nohandpicked' and float(zsup) < 0 and innermask == '5':
        constr_gamma = 0.128
        constrwidth_gamma_inf = 0.143
        constrwidth_gamma_sup = 0.113
    if other == 'fiducial' and (handpicked == 'handpicked' or innermask == '15' or float(zsup) > 0):
        constr_gamma = 0.10
        constrwidth_gamma_inf = 0.085
        constrwidth_gamma_sup = 0.115
    if other == 'fiducial' and (handpicked == 'removegrouphandpicked' or float(zsup) > 0):
        constr_gamma = 0.095
        constrwidth_gamma_inf = 0.070
        constrwidth_gamma_sup = 0.110
    filters = "ugrizJHK"
    print 'shear: ',constr_gamma
if lens == "J1206":
    filters = "griK"
    plane = 34
    constr_gamma = 0.04
    constrwidth_gamma_inf = 0.03
    constrwidth_gamma_sup = 0.05

# declare which weights to read
if radius == "45":
    measured_index = 0 # specifies the column index in weightsfile
    measured_index_inf = 1
    measured_index_sup = 2
if radius == "120":
    measured_index = 3
    measured_index_inf = 4
    measured_index_sup = 5

if lens != "J1206":
    def declareweight(weightin):
        if weightin == "gal": weight_index = 4
        if weightin == "z": weight_index = 6
        if weightin == "mass": weight_index = 8
        if weightin == "mass2": weight_index = 10
        if weightin == "mass3": weight_index = 12
        if weightin == "oneoverr": weight_index = 14
        if weightin == "zoverr": weight_index = 16
        if weightin == "massoverr": weight_index = 18
        if weightin == "mass2overr": weight_index = 20
        if weightin == "mass3overr": weight_index = 22
        if weightin == "mass2rms": weight_index = 24
        if weightin == "mass3rms": weight_index = 26
        if weightin == "mass2overrrms": weight_index = 28
        if weightin == "mass3overrrms": weight_index = 30
        if weightin == "flexion": weight_index = 32
        if weightin == "tidal": weight_index = 34
        if weightin == "SIS": weight_index = 36
        if weightin == "SIShalo": weight_index = 38
        if weightin == "gamma": weight_index = None
        return weight_index
if lens == "J1206":
    def declareweight(weightin):
        if weightin == "gal": weight_index = 4
        if weightin == "z": weight_index = 5
        if weightin == "oneoverr": weight_index = 6
        if weightin == "zoverr": weight_index = 7
        if weightin == "gamma": weight_index = None
        return weight_index

weight1_index = declareweight(weightin1)
if conjoined >= 2:
    weight2_index = declareweight(weightin2)
    if conjoined >= 3:
        weight3_index = declareweight(weightin3)
        if conjoined == 4:
            weight4_index = declareweight(weightin4)

# read weight constraints
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

def declareweight(weightin):
    if weightin == "gal": constr_weight = constr_gal_meds; constrwidth_weight_inf = constrwidth_gal_meds_inf; constrwidth_weight_sup = constrwidth_gal_meds_sup
    if weightin == "z": constr_weight = constr_z_meds; constrwidth_weight_inf = constrwidth_z_meds_inf; constrwidth_weight_sup = constrwidth_z_meds_sup
    if weightin == "mass": constr_weight = constr_mass_meds; constrwidth_weight_inf = constrwidth_mass_meds_inf; constrwidth_weight_sup = constrwidth_mass_meds_sup
    if weightin == "mass2": constr_weight = constr_mass2_meds; constrwidth_weight_inf = constrwidth_mass2_meds_inf; constrwidth_weight_sup = constrwidth_mass2_meds_sup
    if weightin == "mass3": constr_weight = constr_mass3_meds; constrwidth_weight_inf = constrwidth_mass3_meds_inf; constrwidth_weight_sup = constrwidth_mass3_meds_sup
    if weightin == "oneoverr": constr_weight = constr_oneoverr_meds; constrwidth_weight_inf = constrwidth_oneoverr_meds_inf; constrwidth_weight_sup = constrwidth_oneoverr_meds_sup
    if weightin == "zoverr": constr_weight = constr_zoverr_meds; constrwidth_weight_inf = constrwidth_zoverr_meds_inf; constrwidth_weight_sup = constrwidth_zoverr_meds_sup
    if weightin == "massoverr": constr_weight = constr_massoverr_meds; constrwidth_weight_inf = constrwidth_massoverr_meds_inf; constrwidth_weight_sup = constrwidth_massoverr_meds_sup
    if weightin == "mass2overr": constr_weight = constr_mass2overr_meds; constrwidth_weight_inf = constrwidth_mass2overr_meds_inf; constrwidth_weight_sup = constrwidth_mass2overr_meds_sup
    if weightin == "mass3overr": constr_weight = constr_mass3overr_meds; constrwidth_weight_inf = constrwidth_mass3overr_meds_inf; constrwidth_weight_sup = constrwidth_mass3overr_meds_sup
    if weightin == "mass2rms": constr_weight = constr_mass2rms_meds; constrwidth_weight_inf = constrwidth_mass2rms_meds_inf; constrwidth_weight_sup = constrwidth_mass2rms_meds_sup
    if weightin == "mass3rms": constr_weight = constr_mass3rms_meds; constrwidth_weight_inf = constrwidth_mass3rms_meds_inf; constrwidth_weight_sup = constrwidth_mass3rms_meds_sup
    if weightin == "mass2overrrms": constr_weight = constr_mass2overrrms_meds; constrwidth_weight_inf = constrwidth_mass2overrrms_meds_inf; constrwidth_weight_sup = constrwidth_mass2overrrms_meds_sup
    if weightin == "mass3overrrms": constr_weight = constr_mass3overrrms_meds; constrwidth_weight_inf = constrwidth_mass3overrrms_meds_inf; constrwidth_weight_sup = constrwidth_mass3overrrms_meds_sup
    if weightin == "flexion": constr_weight = constr_flexion_meds; constrwidth_weight_inf = constrwidth_flexion_meds_inf; constrwidth_weight_sup = constrwidth_flexion_meds_sup
    if weightin == "tidal": constr_weight = constr_tidal_meds; constrwidth_weight_inf = constrwidth_tidal_meds_inf; constrwidth_weight_sup = constrwidth_tidal_meds_sup
    if weightin == "SIS": constr_weight = constr_SIS_meds; constrwidth_weight_inf = constrwidth_SIS_meds_inf; constrwidth_weight_sup = constrwidth_SIS_meds_sup
    if weightin == "SIShalo": constr_weight = constr_SIShalo_meds; constrwidth_weight_inf = constrwidth_SIShalo_meds_inf; constrwidth_weight_sup = constrwidth_SIShalo_meds_sup
    if weightin == "gamma": constr_weight = constr_gamma; constrwidth_weight_inf = constrwidth_gamma_inf; constrwidth_weight_sup = constrwidth_gamma_sup
    return constr_weight, constrwidth_weight_inf, constrwidth_weight_sup

if conjoined == 4: constr_weight4, constrwidth_weight4_inf, constrwidth_weight4_sup = declareweight(weightin4)
if (conjoined == 3) | (conjoined == 4): constr_weight3, constrwidth_weight3_inf, constrwidth_weight3_sup = declareweight(weightin3)
if (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight2, constrwidth_weight2_inf, constrwidth_weight2_sup = declareweight(weightin2)
if (conjoined == 1) | (conjoined == 2) | (conjoined == 3) | (conjoined == 4): constr_weight1, constrwidth_weight1_inf, constrwidth_weight1_sup = declareweight(weightin1)

print "Reading..."

if conjoined == 4:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s_%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,weightin4,mag,radius,mode,increment1,increment2,increment3,increment4)
if conjoined == 3:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s_%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,weightin3,mag,radius,mode,increment1,increment2,increment3)
if conjoined == 2:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_%s_increments%s_%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,weightin2,mag,radius,mode,increment1,increment2)
if conjoined == 1:
    output = '%skappahist_%s_%s_%sinnermask_nobeta%s_zgap%s_%s_%s_%s_%s_%s_%s_increments%s.cat' % (rootout,lens,compmeas,innermask,handpickedstr,zinf,zsup,other,weightin1,mag,radius,mode,increment1)

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

if conjoined == 1:
    ''' Here I only read the columns of interest, without kappa, for ugriz, in order to find the medians of their values over the whole MS.'''
    med1 = np.zeros(8)
    filters1 = "ugriz"
    for j in range(8):
      for i in range(8):
        if weightin1 != "gamma":
            weight1_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index])
            if i == 0:
                weight1 = weight1_
            else:
                weight1 = np.append(weight1,weight1_)
        else:
            weight1_1_,weight1_2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[2,3])
            if i == 0:
                weight1_1 = weight1_1_
                weight1_2 = weight1_2_
            else:
                weight1_1 = np.append(weight1_1,weight1_1_)
                weight1_2 = np.append(weight1_2,weight1_2_)
        #print j,i
      if weightin1 != "gamma":
        med1[j] = np.median(weight1)
      else:
        med1[j] = np.median(np.sqrt(weight1_1**2 + weight1_2**2))
    med_weight1 = np.mean(med1) # throughout the code I use med_weight1 when computing intervals, following Green et al. For this, weight1 should always refer to simple galaxy number counts
    if weightin1 == "gamma":
        constr_weight1 = constr_weight1 / med_weight1 # for gamma, measured shear divided by the median value of shear in MS; this turns it into an overdensity, like the other weights, so that it is meaningful to multiply by med_weight1
        constrwidth_weight1_inf = constrwidth_weight1_inf / med_weight1
        constrwidth_weight1_sup = constrwidth_weight1_sup / med_weight1
        del weight1_1
        del weight1_1_
        del weight1_2
        del weight1_2_
    else:
        del weight1
        del weight1_
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))]) # absolute number, e.g. of galaxies within the lower width
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    
    ''' Here I read ugrizJHK, converting weighted counts into overdensities, and recording the kappa values only for overdensities satisfying the constraint. I consider the full range of the constraint.'''
    filters1 = filters
    for j in range(8):
      for i in range(8):
        if weightin1 != "gamma":
            kappa_, weight1_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index))
            weight1_ = weight1_ / med_weight1
        else:
            kappa_, gamma1_,gamma2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = gamma / med_weight1
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ] # convert overdensities into absolute counts
        print np.shape(kappa_)
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        if (i == 0) and (j == 0):
            kappa = kappa_
            weight1 = weight1_
        else:
            kappa = np.append(kappa,kappa_)
            weight1 = np.append(weight1,weight1_)
        #print j,i
    if weightin1 == "gamma":
        del weight1_
        del gamma1
        del gamma1_
        del gamma2
        del gamma2_
        del gamma
    else:
        del weight1_

if conjoined == 2:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    filters1 = "ugriz"
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index))
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
        else:
            weight1_,weight2_1_,weight2_2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=[weight1_index,2,3])
            if i == 0:
                weight1 = weight1_
                weight2_1 = weight2_1_
                weight2_2 = weight2_2_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2_1 = np.append(weight2_1,weight2_1_)
                weight2_2 = np.append(weight2_2,weight2_2_)
        #print j,i
      if weightin2 != "gamma":
        med1[j] = np.median(weight1)
        med2[j] = np.median(weight2)
      else:
        med1[j] = np.median(weight1)
        med2[j] = np.median(np.sqrt(weight2_1**2 + weight2_2**2))
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    if weightin2 == "gamma":
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    if weightin2 == "gamma":
        del weight1
        del weight1_
        del weight2_1
        del weight2_1_
        del weight2_2
        del weight2_2_
    else:
        del weight1
        del weight1_
        del weight2
        del weight2_

    filters1 = filters
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
        else:
            kappa_, weight1_,gamma1_,gamma2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        print np.shape(kappa_)
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
        #print j,i
    if weightin1 == "gamma":
        del weight1_
        del weight2_
        del gamma1
        del gamma1_
        del gamma2
        del gamma2_
        del gamma
    else:
        del weight1_
        del weight2_


if conjoined == 3:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    filters1 = "ugriz"
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_,weight3_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index))
            if i == 0:
                weight1 = weight1_
                weight2 = weight2_
                weight3 = weight3_
            else:
                weight1 = np.append(weight1,weight1_)
                weight2 = np.append(weight2,weight2_)
                weight3 = np.append(weight3,weight3_)
        else:
            weight1_,weight2_1_,weight2_2_,weight3_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,2,3,weight3_index))
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
      if weightin2 != "gamma":
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
    if weightin2 == "gamma":
        constr_weight2 = constr_weight2 / med_weight2
        constrwidth_weight2_inf = constrwidth_weight2_inf / med_weight2
        constrwidth_weight2_sup = constrwidth_weight2_sup / med_weight2
    E_w1_inf = np.max([1, round(med_weight1 * (constr_weight1 - constrwidth_weight1_inf))])
    E_w1_sup = np.max([1, round(med_weight1 * (-constr_weight1 + constrwidth_weight1_sup))])
    E_w2_inf = np.max([1, round(med_weight1 * (constr_weight2 - constrwidth_weight2_inf))])
    E_w2_sup = np.max([1, round(med_weight1 * (-constr_weight2 + constrwidth_weight2_sup))])
    E_w3_inf = np.max([1, round(med_weight1 * (constr_weight3 - constrwidth_weight3_inf))])
    E_w3_sup = np.max([1, round(med_weight1 * (-constr_weight3 + constrwidth_weight3_sup))])
    if weightin2 == "gamma":
        del weight1
        del weight1_
        del weight2_1
        del weight2_1_
        del weight2_2
        del weight2_2_
        del weight3
        del weight3_
    else:
        del weight1
        del weight1_
        del weight2
        del weight2_
        del weight3
        del weight3_

    filters1 = filters
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_,weight3_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index,weight3_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
        else:
            kappa_, weight1_,weight3_,gamma1_,gamma2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight3_index,2,3))
            gamma1 = gamma1_
            gamma2 = gamma2_
            gamma = gamma1 # just so that the array has the correct shape
            gamma = np.sqrt(gamma1**2 + gamma2**2)
            weight1_ = weight1_ / med_weight1
            weight2_ = gamma / med_weight2
            weight3_ = weight3_ / med_weight3
        weight = np.copy(weight1_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        del weight
        weight = np.copy(weight3_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        print np.shape(kappa_)
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
        #print j,i
    if weightin1 == "gamma":
        del weight1_
        del weight2_
        del weight3_
        del gamma1
        del gamma1_
        del gamma2
        del gamma2_
        del gamma
    else:
        del weight1_
        del weight2_
        del weight3_

if conjoined == 4:
    med1 = np.zeros(8)
    med2 = np.zeros(8)
    med3 = np.zeros(8)
    med4 = np.zeros(8)
    filters1 = "ugriz"
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            weight1_,weight2_,weight3_,weight4_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,weight2_index,weight3_index,weight4_index))
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
            weight1_,weight2_1_,weight2_2_,weight3_,weight4_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(weight1_index,2,3,weight3_index,weight4_index))
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
    med_weight1 = np.mean(med1)
    med_weight2 = np.mean(med2)
    med_weight3 = np.mean(med3)
    med_weight4 = np.mean(med4)
    if weightin2 == "gamma":
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
    if weightin2 == "gamma":
        del weight1
        del weight1_
        del weight2_1
        del weight2_1_
        del weight2_2
        del weight2_2_
        del weight3
        del weight3_
        del weight4
        del weight4_
    else:
        del weight1
        del weight1_
        del weight2
        del weight2_
        del weight3
        del weight3_
        del weight4
        del weight4_

    filters1 = filters
    for j in range(8):
      for i in range(8):
        if weightin2 != "gamma":
            kappa_, weight1_,weight2_,weight3_,weight4_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight2_index,weight3_index,weight4_index))
            weight1_ = weight1_ / med_weight1
            weight2_ = weight2_ / med_weight2
            weight3_ = weight3_ / med_weight3
            weight4_ = weight4_ / med_weight4
        else:
            kappa_, weight1_,weight3_,weight4_,gamma1_,gamma2_ = readfile("%snobeta%s%s%sinject_%s_%s_GGL_los_8_%s_%s_%s_%s_%sarcsecinner_gap_%s_%s.fits" % (root,str(plane),compmeas,mode,filters1,lens,str(j),str(i),mag,radius,innermask,zinf,zsup), usecols=(1,weight1_index,weight3_index,weight4_index,2,3))
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
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight1 * med_weight1) - limsigma * E_w1_inf - increment1/2.0) & (weight * med_weight1 < round(constr_weight1 * med_weight1) + limsigma * E_w1_sup + increment1/2.0) ]
        del weight
        weight = np.copy(weight2_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight2 * med_weight1) - limsigma * E_w2_inf - increment2/2.0) & (weight * med_weight1 < round(constr_weight2 * med_weight1) + limsigma * E_w2_sup + increment2/2.0) ]
        del weight
        weight = np.copy(weight3_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight3 * med_weight1) - limsigma * E_w3_inf - increment3/2.0) & (weight * med_weight1 < round(constr_weight3 * med_weight1) + limsigma * E_w3_sup + increment3/2.0) ]
        del weight
        weight = np.copy(weight4_)
        print np.shape(kappa_)
        kappa_ = kappa_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        print np.shape(kappa_)
        weight1_ = weight1_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight2_ = weight2_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight3_ = weight3_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
        weight4_ = weight4_[(weight * med_weight1 >= round(constr_weight4 * med_weight1) - limsigma * E_w4_inf - increment4/2.0) & (weight * med_weight1 < round(constr_weight4 * med_weight1) + limsigma * E_w4_sup + increment4/2.0) ]
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
        #print j,i
    if weightin1 == "gamma":
        del weight1_
        del weight2_
        del weight3_
        del weight4_
        del gamma1
        del gamma1_
        del gamma2
        del gamma2_
        del gamma
    else:
        del weight1_
        del weight2_
        del weight3_
        del weight4_

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

head = 'LOS: %d' % np.array([LOS])
np.savetxt(output,unbiased_kappa_constrained,header=head,fmt='%s',delimiter='\t',newline='\n')
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
