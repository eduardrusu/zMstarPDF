# run this in order to sample from P(Mstar,z) for the lens catalogue. It requires that the output files produced by photozsampling.py first be used as input for LePhare, and then combined with combinelephare_withredshift.py

import numpy as np
import scipy
from scipy import stats
import sys
import os
from os import system

# Behroozi et al 2010 parameters for z < 1:
M10_  = np.array([12.35,-0.16,+0.07])
M1a_  = np.array([0.28,-0.97,+0.19])
Ms00_ = np.array([10.72,-0.29,+0.22])
Ms0a_ = np.array([0.55,-0.79,+0.18])
b0_   = np.array([0.44,-0.06,+0.04])
ba_   = np.array([0.18,-0.34,+0.08])
d0_   = np.array([0.57,-0.06,+0.15])
da_   = np.array([0.17,-0.41,+0.42])
g0_   = np.array([1.56,-0.38,+0.12])
ga_   = np.array([2.51,-1.83,+0.15])
# z >= 1:
M10   = np.array([12.27,-0.27,+0.59])
M1a   = np.array([-0.84,-0.58,+0.87])
Ms00  = np.array([11.09,-0.31,+0.54])
Ms0a  = np.array([0.56,-0.44,+0.89])
b0    = np.array([0.65,-0.20,+0.26])
ba    = np.array([0.31,-0.47,+0.38])
d0    = np.array([0.56,-0.29,+1.33])
da    = np.array([-0.12,-0.50,+0.76])
g0    = np.array([1.12,-0.36,+7.47])
ga    = np.array([-0.53,-2.50,+7.87])

def sample(median,stdlow,stdhigh): # samples from different standard deviation gaussians on each side of the mean
    rand = np.random.uniform(0,1,1)[0]
    if rand <= 0.5: # since MED divides the distribution in two equal parts
        return median - np.abs(np.random.normal(0, np.abs(stdlow), 1)[0])
    else:
        return median + np.abs(np.random.normal(0, stdhigh, 1)[0])

samples = 20
masterfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_bpzsample_combined.cat.MAG_BC03_I09.lephareout"
massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_eazysample_combined.cat.MAG_BC03_I09.lephareout"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_noIRACbpzsample_combined.cat.MAG_BC03_I09.lephareout"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_noIRACeazysample_combined.cat.MAG_BC03_I09.lephareout"
#masterfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_bpzsample_combined.cat.MAG_BC03_I09.lephareout"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_eazysample_combined.cat.MAG_BC03_I09.lephareout"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_noIRACbpzsample_combined.cat.MAG_BC03_I09.lephareout"
#massfile = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_noIRACeazysample_combined.cat.MAG_BC03_I09.lephareout"


x = 0
y = 1
ra = 2
dec = 3
i_auto = 4
i_auto_err = 5
id = 8
z_b = 28
z_m2 = 48
spec = 40
#mBEST_ind_bpz = 60 # no IRAC
#mINF_ind_bpz = 61 # no IRAC
#mMED_ind_bpz = 62 # no IRAC
#mSUP_ind_bpz = 63 # no IRAC
#mBEST_ind_eazy = 67 # no IRAC
#mINF_ind_eazy = 68 # no IRAC
#mMED_ind_eazy = 69 # no IRAC
#mSUP_ind_eazy = 70 # no IRAC
#class_bpz = 72 # no IRAC
#class_eazy = 73 # no IRAC
mBEST_ind_bpz = 85
mINF_ind_bpz = 86
mMED_ind_bpz = 87
mSUP_ind_bpz = 88
mBEST_ind_eazy = 92
mINF_ind_eazy = 93
mMED_ind_eazy = 94
mSUP_ind_eazy = 95
class_bpz = 97
class_eazy = 98
if "bpzsample" in massfile: classify = class_bpz
else: classify = class_eazy

master = np.loadtxt(masterfile,usecols=[x,y,ra,dec,i_auto,i_auto_err,id,classify,spec,z_b,z_m2,mBEST_ind_bpz,mINF_ind_bpz,mMED_ind_bpz,mSUP_ind_bpz,mBEST_ind_eazy,mINF_ind_eazy,mMED_ind_eazy,mSUP_ind_eazy],unpack=False)
mass = np.loadtxt(massfile,usecols=[0,1,4,5,6,7]) # index, redshift and masses

# Re-label columns
x = 0
y = 1
ra = 2
dec = 3
i_auto = 4
i_auto_err = 5
id = 6
classify = 7
spec = 8
z_b = 9
z_m2 = 10
mBEST_ind_bpz = 11
mINF_ind_bpz = 12
mMED_ind_bpz = 13
mSUP_ind_bpz = 14
mBEST_ind_eazy = 15
mINF_ind_eazy = 16
mMED_ind_eazy = 17
mSUP_ind_eazy = 18

if "bpzsample" in massfile:
    z_ind = z_b
    mBEST_ind = mBEST_ind_bpz
    mINF_ind = mINF_ind_bpz
    mMED_ind = mMED_ind_bpz
    mSUP_ind = mSUP_ind_bpz
else:
    z_ind = z_m2
    mBEST_ind = mBEST_ind_eazy
    mINF_ind = mINF_ind_eazy
    mMED_ind = mMED_ind_eazy
    mSUP_ind = mSUP_ind_eazy

# for missing mMED take its error bars to be typical for that mag range
errbar18 = np.median(master[:,mMED_ind][(master[:,i_auto]>15) & (master[:,i_auto]<18) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>15) & (master[:,i_auto]<18) & (master[:,mMED_ind]!=-99)])
errbar19 = np.median(master[:,mMED_ind][(master[:,i_auto]>18) & (master[:,i_auto]<19) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>18) & (master[:,i_auto]<19) & (master[:,mMED_ind]!=-99)])
errbar20 = np.median(master[:,mMED_ind][(master[:,i_auto]>19) & (master[:,i_auto]<20) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>19) & (master[:,i_auto]<20) & (master[:,mMED_ind]!=-99)])
errbar21 = np.median(master[:,mMED_ind][(master[:,i_auto]>20) & (master[:,i_auto]<21) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>20) & (master[:,i_auto]<21) & (master[:,mMED_ind]!=-99)])
errbar22 = np.median(master[:,mMED_ind][(master[:,i_auto]>21) & (master[:,i_auto]<22) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>21) & (master[:,i_auto]<22) & (master[:,mMED_ind]!=-99)])
errbar23 = np.median(master[:,mMED_ind][(master[:,i_auto]>21) & (master[:,mMED_ind]!=-99)] - master[:,mINF_ind][(master[:,i_auto]>21) & (master[:,mMED_ind]!=-99)])

for i in range(master.shape[0]):
        line = np.zeros(8 + (samples + 1) * 3) # x,y,ra,dec,i_auto,i_auto_err,id,classify + 10 * (z,Mstar,Mhalo)
        line[0:8] = master[i][0:8] # x,y,ra,dec,i_auto,i_auto_err,id,classify are the first 8 columns from master
        z = np.zeros(samples + 1)
        if master[i][spec] > 0: z[0] = master[i][spec] # if there is a galaxy spectrum
        else: z[0] = master[i][z_ind]

        mBEST = master[i][mBEST_ind]
        mINF = master[i][mINF_ind]
        mMED = master[i][mMED_ind]
        mSUP = master[i][mSUP_ind]
        if mBEST < 0:
            mBEST = 9
        if mMED < 0:
            mMED = mBEST
            if (master[i][i_auto] > 15) & (master[i][i_auto] <= 18): errbar = errbar18
            if (master[i][i_auto] > 18) & (master[i][i_auto] <= 19): errbar = errbar19
            if (master[i][i_auto] > 19) & (master[i][i_auto] <= 20): errbar = errbar20
            if (master[i][i_auto] > 20) & (master[i][i_auto] <= 21): errbar = errbar21
            if (master[i][i_auto] > 21) & (master[i][i_auto] <= 22): errbar = errbar22
            if master[i][i_auto] > 22: errbar = errbar23
            mINF = mBEST - errbar
            mSUP = mBEST + errbar
        m = np.zeros(samples + 1)
        mhalo = np.zeros(samples + 1)
        m[0] = mMED
        a = 1 / (1 + z[0])
        if z[0] <= 1:
            logM1a = M10_[0] + M1a_[0] * (a - 1)
            logMs0a = Ms00_[0] + Ms0a_[0] * (a - 1)
            notlogMs0a = 10 ** logMs0a
            b = b0_[0] + ba_[0] * (a-1)
            d = d0_[0] + da_[0] * (a-1)
            g = g0_[0] + ga_[0] * (a-1)
            mhalo[0] = logM1a + b * (m[0] - logMs0a) + ((10 ** m[0]/notlogMs0a)**d)/(1+(10 ** m[0]/notlogMs0a)**(-g)) - 1/2
        else:
            logM1a = M10[0] + M1a[0] * (a - 1)
            logMs0a = Ms00[0] + Ms0a[0] * (a - 1)
            notlogMs0a = 10 ** logMs0a
            b = b0[0] + ba[0] * (a-1)
            d = d0[0] + da[0] * (a-1)
            g = g0[0] + ga[0] * (a-1)
            mhalo[0] = logM1a + b * (m[0] - logMs0a) + ((10 ** m[0]/notlogMs0a)**d)/(1+(10 ** m[0]/notlogMs0a)**(-g)) - 1/2
        line[8:11] = np.array([z[0],m[0],mhalo[0]])

        ind = 0
        for j in range(mass.shape[0]):
                if mass[j][0] == master[i][id]:
                    ind = ind + 1 # this goes from 0 to 8
                    z[ind] = mass[j][1]

                    rand = np.random.uniform(0,1,1)[0]
                    mBEST = mass[j][2]
                    mINF = mass[j][3]
                    mMED = mass[j][4]
                    mSUP = mass[j][5]
                    if mBEST < 0:
                        mBEST = 9
                    if mMED < 0:
                        mMED = mBEST
                        mINF = mBEST - errbar
                        mSUP = mBEST + errbar
                    m[ind] = sample(mMED,mINF-mMED,mSUP - mMED)
                    if z[ind] <= 1:
                        logM1a = sample(M10_[0],M10_[1],M10_[2]) + sample(M1a_[0],M1a_[1],M1a_[2]) * (a - 1)
                        logMs0a = sample(Ms00_[0],Ms00_[1],Ms00_[2]) + sample(Ms0a_[0],Ms0a_[1],Ms0a_[2]) * (a - 1)
                        notlogMs0a = 10 ** logMs0a
                        b = sample(b0_[0],b0_[1],b0_[2]) + sample(ba_[0],ba_[1],ba_[2]) * (a - 1)
                        d = sample(d0_[0],d0_[1],d0_[2]) + sample(da_[0],da_[1],da_[2]) * (a - 1)
                        g = sample(g0_[0],g0_[1],g0_[2]) + sample(ga_[0],ga_[1],ga_[2]) * (a - 1)
                    else:
                        logM1a = sample(M10[0],M10[1],M10[2]) + sample(M1a[0],M1a[1],M1a[2]) * (a - 1)
                        logMs0a = sample(Ms00[0],Ms00[1],Ms00[2]) + sample(Ms0a[0],Ms0a[1],Ms0a[2]) * (a - 1)
                        notlogMs0a = 10 ** logMs0a
                        b = sample(b0[0],b0[1],b0_[2]) + sample(ba[0],ba[1],ba[2]) * (a - 1)
                        d = sample(d0[0],d0[1],d0_[2]) + sample(da[0],da[1],da[2]) * (a - 1)
                        g = sample(g0[0],g0[1],g0_[2]) + sample(ga[0],ga[1],ga[2]) * (a - 1)
                    mhalo[ind] = logM1a + b * (m[ind] - logMs0a) + ((10 ** m[ind]/notlogMs0a)**d)/(1+(10 ** m[ind]/notlogMs0a)**(-g)) - 1/2

        for j in range(samples):
            line[8+(j+1)*3:8+(j+2)*3] = np.array([z[j+1],m[j+1],mhalo[j+1]])
        if i == 0:
            data = line
            data[0:8] = master[0][0:8]
        else:
            data = np.c_[data,line]

#np.savetxt(masterfile[:-4] + "_WFI2033noIRACbpz_nobeta.cat",data.T,fmt='%s %s %s %s %.2f %.2f %d %d %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f')
#np.savetxt(masterfile[:-4] + "_WFI2033noIRACeazy_nobeta.cat",data.T,fmt='%s %s %s %s %.2f %.2f %d %d %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f')
#np.savetxt(masterfile[:-4] + "_WFI2033IRACbpz_nobeta.cat",data.T,fmt='%s %s %s %s %.2f %.2f %d %d %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f')
np.savetxt(masterfile[:-4] + "_WFI2033IRACeazy_nobeta.cat",data.T,fmt='%s %s %s %s %.2f %.2f %d %d %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f')
