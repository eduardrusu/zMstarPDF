# Run this code in order to sample from P(z) for the lens catalogue, with or without spectroscopic information. Requires that a necessary .probs file exists in the BPZ folfer, containing P(z) computed with BPZ on a grid z=arange(0.0000,3.5100,0.0100). Also requires that the .pz files exist for each object with z=arange(0.0100,4.0000,0.0100) from a previous EAzY run on a grid z=arange(0.0100,4.00,0.0100). The code then uses a modified version of converttolephare_WFI2033.py and converttolephare_noobs_WFI2033.py in order to produce the input expected by Lephare.

import numpy as np
import scipy
from scipy import stats
import sys
import os
from os import system
import time

useeazy = 1
samples = 20
#file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
#filebpz = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.probs"
#if useeazy == 1: fileeazy = "/Users/cerusu/GITHUB/eazy-photoz/inputs/OUTPUT/sample_ir/"
file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
filebpz = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.probs"
if useeazy == 1: fileeazy = "/Users/cerusu/GITHUB/eazy-photoz/inputs/OUTPUT/sample_i/"

# zeropoint corrections suggested by lephare:
u_corr = +1.40
g_corr = -0.03
r_corr = -0.00
i_corr =  0.00
z_corr = -0.01
Y_corr = -0.07
J_corr = -0.10
H_corr = -0.01
Ks_corr = +0.06

# conversion from Vega to AB; assuming that the input mags are in Vega
J_corr += 0.94 # computed by LePhare
H_corr += 1.35
Ks_corr += 1.83

# use the modified version of converttolephare_WFI2033.py:
def lephare_noobs(data,fileout,irac):

    # a small number of objects in BPZ have good mags, but error on the mag 99; those objects should be good, and Ideally I would fix the errors one by one through closer examination. Here I just replace their errors with 1 mag
    data[u_err][np.abs(data[u_err]) == 99.00] = 1.00
    data[g_err][np.abs(data[g_err]) == 99.00] = 1.00
    data[r_err][np.abs(data[r_err]) == 99.00] = 1.00
    data[i_err][np.abs(data[i_err]) == 99.00] = 1.00
    data[z_err][np.abs(data[z_err]) == 99.00] = 1.00
    data[Y_err][np.abs(data[Y_err]) == 99.00] = 1.00
    data[J_err][np.abs(data[J_err]) == 99.00] = 1.00
    data[H_err][np.abs(data[H_err]) == 99.00] = 1.00
    data[Ks_err][np.abs(data[Ks_err]) == 99.00] = 1.00

    if irac == False:
    # If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
        data[irac1] = -99
        data[irac1_err] = -99
        data[irac2] = -99
        data[irac2_err] = -99
        data[irac3] = -99
        data[irac3_err] = -99
        data[irac4] = -99
        data[irac4_err] = -99

    data[u_err][np.abs(data[u]) == 99.0] = -99.0
    data[g_err][np.abs(data[g]) == 99.0] = -99.0
    data[r_err][np.abs(data[r]) == 99.0] = -99.0
    data[i_err][np.abs(data[i]) == 99.0] = -99.0
    data[z_err][np.abs(data[z]) == 99.0] = -99.0
    data[Y_err][np.abs(data[Y]) == 99.0] = -99.0
    data[J_err][np.abs(data[J]) == 99.0] = -99.0
    data[H_err][np.abs(data[H]) == 99.0] = -99.0
    data[Ks_err][np.abs(data[Ks]) == 99.0] = -99.0
    data[u][np.abs(data[u]) == 99.0] = -99.0
    data[g][np.abs(data[g]) == 99.0] = -99.0
    data[r][np.abs(data[r]) == 99.0] = -99.0
    data[i][np.abs(data[i]) == 99.0] = -99.0
    data[z][np.abs(data[z]) == 99.0] = -99.0
    data[Y][np.abs(data[Y]) == 99.0] = -99.0
    data[J][np.abs(data[J]) == 99.0] = -99.0
    data[H][np.abs(data[H]) == 99.0] = -99.0
    data[Ks][np.abs(data[Ks]) == 99.0] = -99.0

    # apply the corrections
    data[u][np.abs(data[u]) != 99.0] += u_corr
    data[g][np.abs(data[g]) != 99.0] += g_corr
    data[r][np.abs(data[r]) != 99.0] += r_corr
    data[i][np.abs(data[i]) != 99.0] += i_corr
    data[z][np.abs(data[z]) != 99.0] += z_corr
    data[Y][np.abs(data[Y]) != 99.0] += Y_corr
    data[J][np.abs(data[J]) != 99.0] += J_corr
    data[H][np.abs(data[H]) != 99.0] += H_corr
    data[Ks][np.abs(data[Ks]) != 99.0] += Ks_corr

    # LePhare thinks error bars = 0 means non-detection, so fix this
    data[u_err][data[u_err] == 0.00] = 0.01
    data[g_err][data[g_err] == 0.00] = 0.01
    data[r_err][data[r_err] == 0.00] = 0.01
    data[i_err][data[i_err] == 0.00] = 0.01
    data[z_err][data[z_err] == 0.00] = 0.01
    data[Y_err][data[Y_err] == 0.00] = 0.01
    data[J_err][data[J_err] == 0.00] = 0.01
    data[H_err][data[H_err] == 0.00] = 0.01
    data[Ks_err][data[Ks_err] == 0.00] = 0.01

    str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t Ks_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
    dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[Ks],data[Ks_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[redshift]]
    np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')

# use the modified version of converttolephare_noobs_WFI2033.py:
def lephare(data,fileout,irac):
    if irac == False:
    # If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
        data[irac1] = -99
        data[irac1_err] = -99
        data[irac2] = -99
        data[irac2_err] = -99
        data[irac3] = -99
        data[irac3_err] = -99
        data[irac4] = -99
        data[irac4_err] = -99

    # a small number of objects in BPZ have good mags, but error on the mag 99; those objects should be good, and Ideally I would fix the errors one by one through closer examination. Here I just replace their errors with 1 mag
    data[u_err][np.abs(data[u_err]) == 99.00] = 1.00
    data[g_err][np.abs(data[g_err]) == 99.00] = 1.00
    data[r_err][np.abs(data[r_err]) == 99.00] = 1.00
    data[i_err][np.abs(data[i_err]) == 99.00] = 1.00
    data[z_err][np.abs(data[z_err]) == 99.00] = 1.00
    data[Y_err][np.abs(data[Y_err]) == 99.00] = 1.00
    data[J_err][np.abs(data[J_err]) == 99.00] = 1.00
    data[H_err][np.abs(data[H_err]) == 99.00] = 1.00
    data[Ks_err][np.abs(data[Ks_err]) == 99.00] = 1.00

    # apply the corrections
    data[u][np.abs(data[u]) != 99.0] += u_corr
    data[g][np.abs(data[g]) != 99.0] += g_corr
    data[r][np.abs(data[r]) != 99.0] += r_corr
    data[i][np.abs(data[i]) != 99.0] += i_corr
    data[z][np.abs(data[z]) != 99.0] += z_corr
    data[Y][np.abs(data[Y]) != 99.0] += Y_corr
    data[J][np.abs(data[J]) != 99.0] += J_corr
    data[H][np.abs(data[H]) != 99.0] += H_corr
    data[Ks][np.abs(data[Ks]) != 99.0] += Ks_corr

    # correct the limiting mags
    data[u_err][np.abs(data[u]) == 99.0] += u_corr
    data[g_err][np.abs(data[g]) == 99.0] += g_corr
    data[r_err][np.abs(data[r]) == 99.0] += r_corr
    data[i_err][np.abs(data[i]) == 99.0] += i_corr
    data[z_err][np.abs(data[z]) == 99.0] += z_corr
    data[Y_err][np.abs(data[Y]) == 99.0] += Y_corr
    data[J_err][np.abs(data[J]) == 99.0] += J_corr
    data[H_err][np.abs(data[H]) == 99.0] += H_corr
    data[Ks_err][np.abs(data[Ks]) == 99.0] += Ks_corr

    # the format for nondetections is error=-1.0 and magnitude at 1-sigma
    data[u][np.abs(data[u]) == 99.0] = data[u_err][np.abs(data[u]) == 99.0]
    data[g][np.abs(data[g]) == 99.0] = data[g_err][np.abs(data[g]) == 99.0]
    data[r][np.abs(data[r]) == 99.0] = data[r_err][np.abs(data[r]) == 99.0]
    data[i][np.abs(data[i]) == 99.0] = data[i_err][np.abs(data[i]) == 99.0]
    data[z][np.abs(data[z]) == 99.0] = data[z_err][np.abs(data[z]) == 99.0]
    data[Y][np.abs(data[Y]) == 99.0] = data[Y_err][np.abs(data[Y]) == 99.0]
    data[J][np.abs(data[J]) == 99.0] = data[J_err][np.abs(data[J]) == 99.0]
    data[H][np.abs(data[H]) == 99.0] = data[H_err][np.abs(data[H]) == 99.0]
    data[Ks][np.abs(data[Ks]) == 99.0] = data[Ks_err][np.abs(data[Ks]) == 99.0]
    data[u_err][np.abs(data[u_err]) > 20] = -1.0
    data[g_err][np.abs(data[g_err]) > 20] = -1.0
    data[r_err][np.abs(data[r_err]) > 20] = -1.0
    data[i_err][np.abs(data[i_err]) > 20] = -1.0
    data[z_err][np.abs(data[z_err]) > 20] = -1.0
    data[Y_err][np.abs(data[Y_err]) > 20] = -1.0
    data[J_err][np.abs(data[J_err]) > 20] = -1.0
    data[H_err][np.abs(data[H_err]) > 20] = -1.0
    data[Ks_err][np.abs(data[Ks_err]) > 20] = -1.0

    # LePhare thinks error bars = 0 means non-detection, so fix this
    data[u_err][data[u_err] == 0.00] = 0.01
    data[g_err][data[g_err] == 0.00] = 0.01
    data[r_err][data[r_err] == 0.00] = 0.01
    data[i_err][data[i_err] == 0.00] = 0.01
    data[z_err][data[z_err] == 0.00] = 0.01
    data[Y_err][data[Y_err] == 0.00] = 0.01
    data[J_err][data[J_err] == 0.00] = 0.01
    data[H_err][data[H_err] == 0.00] = 0.01
    data[Ks_err][data[Ks_err] == 0.00] = 0.01

    str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t Ks_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
    dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[Ks],data[Ks_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[redshift]]
    np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')


###################
# sampling from BPZ:

pdz = np.loadtxt(filebpz, unpack=False)
samplez = np.zeros((pdz.shape[0],samples)) # sample 9 times, the 10th will be the original best-fit
zgridint = np.arange(350) # integers corresponding to z=arange(0.0000,3.5100,0.0100)

id = 8
spec = 40 # -1: no spec data; -2: spec star; >0: available spec
u = 9
u_err = 10
g = 11
g_err = 12
r = 13
r_err = 14
i = 15
i_err = 16
z = 17
z_err = 18
Y = 19
Y_err = 20
J = 21
J_err = 22
H = 23
H_err = 24
Ks = 25
Ks_err = 26
irac_1 = 74
irac_1_err = 75
irac_2 = 76
irac_2_err = 77
irac_3 = 78
irac_3_err = 79
irac_4 = 80
irac_4_err = 81

phot = np.loadtxt(file, usecols = [id,spec,u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,Ks,Ks_err,irac_1,irac_1_err,irac_2,irac_2_err,irac_3,irac_3_err,irac_4,irac_4_err], unpack = False)

# now relabel the columns:
id = 0
spec = 1
u = 2
u_err = 3
g = 4
g_err = 5
r = 6
r_err = 7
i = 8
i_err = 9
z = 10
z_err = 11
Y = 12
Y_err = 13
J = 14
J_err = 15
H = 16
H_err = 17
Ks = 18
Ks_err = 19
irac1 = 20
irac1_err = 21
irac2 = 22
irac2_err = 23
irac3 = 24
irac3_err = 25
irac4 = 26
irac4_err = 27
redshift = 28 # the redshift needed by lephare

for k in range(pdz.shape[0]): # for each galaxy
    if np.sum(pdz[k][1:]) != 1: # not all BPZ probabilities are perfectly normalized
        l = 1
        while pdz[k][l] != 0: l += 1
        pdz[k][l] = 1 - np.sum(pdz[k][1:]) # the first instance where the probability is zero, add the require offset to have a perfect normalization
    custm = stats.rv_discrete(name='custm', values=(zgridint, pdz[k][1:])) # ignore the first column, which is the ID
    sample = custm.rvs(size = samples)
    if phot[k][spec] < 0: # if no spectrum is available or if the object is a spectroscopic star
        samplez[k] = np.array([np.max([0.01 * sample[0],0.01]),np.max([0.01 * sample[1],0.01]),np.max([0.01 * sample[2],0.01]),np.max([0.01 * sample[3],0.01]),np.max([0.01 * sample[4],0.01]),np.max([0.01 * sample[5],0.01]),np.max([0.01 * sample[6],0.01]),np.max([0.01 * sample[7],0.01]),np.max([0.01 * sample[8],0.01]),np.max([0.01 * sample[9],0.01]),np.max([0.01 * sample[10],0.01]),np.max([0.01 * sample[11],0.01]),np.max([0.01 * sample[12],0.01]),np.max([0.01 * sample[13],0.01]),np.max([0.01 * sample[14],0.01]),np.max([0.01 * sample[15],0.01]),np.max([0.01 * sample[16],0.01]),np.max([0.01 * sample[17],0.01]),np.max([0.01 * sample[18],0.01]),np.max([0.01 * sample[19],0.01])])
        # because 0.01 is the redshift step; the minimum redshift accepted by LePhare is 0.01
    else:
        samplez[k] = np.array([phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec]])
    for j in range(phot.shape[0]): # match the two files by ID
        if pdz[k][0] == phot[j][0]:
            if k != 0:
                lephdata = np.c_[lephdata,phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j]]
            else:
                lephdata = np.c_[phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j],phot[j]]

samplez = samplez.reshape(samplez.size)
#print np.shape(lephdata)
lephdata = np.r_[lephdata,samplez.reshape(1,samplez.size)]
#print np.shape(lephdata)

data = np.copy(lephdata)
lephare(data,file[:-4] + "_bpzsample.cat",True)
data = np.copy(lephdata)
lephare(data,file[:-4] + "_noIRACbpzsample.cat",False)
data = np.copy(lephdata)
lephare_noobs(data,file[:-4] + "_bpzsample_noobs.cat",True)
data = np.copy(lephdata)
lephare_noobs(data,file[:-4] + "_noIRACbpzsample_noobs.cat",False)


###################
# sampling from EAzY:
if useeazy == 1:
    del lephdata
    lst = [x for x in os.listdir('%s' %fileeazy) if ('.pz' in x)]
    samplez = np.zeros((len(lst),samples)) # sample 9 times

    gal = np.zeros(len(lst)) # array of all the ID of galaxies
    for k in range(len(lst)):
        gal[k] = int(lst[k].strip('.pz'))
    gal = gal.astype(int)

    zgridint = np.arange(400)
    nr = 0
    for k in range(phot.shape[0]):
        if int(phot[k][0]) in gal:
            pdz = np.loadtxt('%s%s.pz' % (fileeazy,int(phot[k][0])), usecols = [3], unpack=True)
            pdz = pdz/np.sum(pdz) # it needs to be normalized
            custm = stats.rv_discrete(name='custm', values=(zgridint, pdz)) # ignore the first column, which is the ID
            sample = custm.rvs(size = samples)
            if phot[k][spec] < 0: # if no spectrum is available or if the object is a spectroscopic star
                samplez[nr] = np.array([0.01 + 0.01 * sample[0],0.01 + 0.01 * sample[1],0.01 + 0.01 * sample[2],0.01 + 0.01 * sample[3],0.01 + 0.01 * sample[4],0.01 + 0.01 * sample[5],0.01 + 0.01 * sample[6],0.01 + 0.01 * sample[7],0.01 + 0.01 * sample[8],0.01 + 0.01 * sample[9],0.01 + 0.01 * sample[10],0.01 + 0.01 * sample[11],0.01 + 0.01 * sample[12],0.01 + 0.01 * sample[13],0.01 + 0.01 * sample[14],0.01 + 0.01 * sample[15],0.01 + 0.01 * sample[16],0.01 + 0.01 * sample[17],0.01 + 0.01 * sample[18],0.01 + 0.01 * sample[19]]) # because 0.01 is the redshift step
            else:
                samplez[nr] = np.array([phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec],phot[k][spec]])
            if nr != 0:
                lephdata = np.c_[lephdata,phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k]]
            else:
                lephdata = np.c_[phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k],phot[k]]
            nr = nr + 1

    samplez = samplez.reshape(samplez.size)
    print np.shape(lephdata)
    lephdata = np.r_[lephdata,samplez.reshape(1,samplez.size)]

    data = np.copy(lephdata)
    lephare(data,file[:-4] + "_eazysample.cat",True)
    data = np.copy(lephdata)
    lephare(data,file[:-4] + "_noIRACeazysample.cat",False)
    data = np.copy(lephdata)
    lephare_noobs(data,file[:-4] + "_eazysample_noobs.cat",True)
    data = np.copy(lephdata)
    lephare_noobs(data,file[:-4] + "_noIRACeazysample_noobs.cat",False)
