# The code calculates realistic BPZ photoz for galaxies in the Millenium Simulation. It then creates the necessary files to run LePhare on SLAC. The code expects input files created by extractMillenium.py
# run as: python photozMillenium.py GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_J1206.images.txt

import numpy as np
import sys
import os
from os import system
import time

start_timefield = time.time()

# grep -v 99.00 GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_J1206.images.txt > GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_J1206noupperlimits.images.txt # this is for the case when I am calibrating the zpt, so I should not use upper limits
root_bpz = "/Users/cerusu/bpz-1.99.3/test/"
root_original = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/0408_SA_gal_sampledphot/"
file = str(sys.argv[1])
file1 = root_bpz+file[:-4]+'.cat'
os.system("cp %s %s" % (root_original+file,file1))
if "griK" in file: os.system("cp %smillennium_griK.columns %s" % (root_bpz,file1[:-4]+'.columns'))
if "griz" in file: os.system("cp %smillennium_griz.columns %s" % (root_bpz,file1[:-4]+'.columns'))

os.system("python $BPZPATH/bpz.py %s -INTERP 2" % file1)
#os.system("python $BPZPATH/bpz.py %s -INTERP 2 -ONLY_TYPE yes" % file1)
os.system("python $BPZPATH/bpzfinalize.py %s" % file1[:-4])

id = 0
#u = 1
#u_err = 2
g = 1
g_err = 2
r = 3
r_err = 4
i = 5
i_err = 6
z = 7
z_err = 8
#J = 11
#J_err = 12
#H = 13
#H_err = 14
#K = 7
#K_err = 8

if "griK" in file:
    data = np.loadtxt(root_original+file,usecols=[0,7,8,9,10,11,12,13,14],unpack=True) # ID + mags
    data_bpz = np.loadtxt(file1[:-4]+"_bpz.cat",usecols=[9,1],unpack=True)
    ''' Write the mags into LePhare-expected format, and assuming no observations for mags below detection threshold '''
    #u_corr =  0.0 # for griK these are the offsets suggested by running BPZ on one simulation field; scaled so that I apply the same i-band offset as for ugrizJHK
    g_corr =  0.010
    r_corr =  0.010
    i_corr =  0.016
    #z_corr =  0.0
    #J_corr =  0.0
    #H_corr =  0.0
    K_corr = -0.012
    # If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
    #data[u_err] = -99.0
    data[g_err][np.abs(data[g]) == 99.0] = -99.0
    data[r_err][np.abs(data[r]) == 99.0] = -99.0
    data[i_err][np.abs(data[i]) == 99.0] = -99.0
    #data[z_err] = -99.0
    #data[J_err] = -99.0
    #data[H_err] = -99.0
    data[K_err][np.abs(data[K]) == 99.0] = -99.0
    #data[u] = -99.0
    data[g][np.abs(data[g]) == 99.0] = -99.0
    data[r][np.abs(data[r]) == 99.0] = -99.0
    data[i][np.abs(data[i]) == 99.0] = -99.0
    #data[z] = -99.0
    #data[J] = -99.0
    #data[H] = -99.0
    data[K][np.abs(data[K]) == 99.0] = -99.0
    # apply the corrections suggested by BPZ
    #data[u][np.abs(data[u]) != 99.0] += u_corr
    data[g][np.abs(data[g]) != 99.0] += g_corr
    data[r][np.abs(data[r]) != 99.0] += r_corr
    data[i][np.abs(data[i]) != 99.0] += i_corr
    #data[z][np.abs(data[z]) != 99.0] += z_corr
    #data[J][np.abs(data[J]) != 99.0] += J_corr
    #data[H][np.abs(data[H]) != 99.0] += H_corr
    data[K][np.abs(data[K]) != 99.0] += K_corr
    # LePhare thinks error bars = 0 means non-detection, so fix this
    #data[u_err][data[u_err] == 0.00] = 0.01
    data[g_err][data[g_err] == 0.00] = 0.01
    data[r_err][data[r_err] == 0.00] = 0.01
    data[i_err][data[i_err] == 0.00] = 0.01
    #data[z_err][data[z_err] == 0.00] = 0.01
    #data[J_err][data[J_err] == 0.00] = 0.01
    #data[H_err][data[H_err] == 0.00] = 0.01
    data[K_err][data[K_err] == 0.00] = 0.01
    fileout = root_original + file[:-4] + "_forlephare.txt"
    #str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t K \t K_err \t specz \t photoz"
    str = "ID \t g \t g_err \t r \t r_err \t i \t i_err \t K \t K_err \t specz \t photoz"
    dataout = np.c_[data[id],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[K],data[K_err],data_bpz[0],data_bpz[1]]
    #np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.3f \t %.2f')

# If I don't require to use Lephare
#data = np.loadtxt(root_original+file,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],unpack=True) # ID + mags
data = np.loadtxt(root_original+file,usecols=[0,1,2,3,8],unpack=True) # ID + mags
data_bpz = np.loadtxt(file1[:-4]+"_bpz.cat",usecols=[1],unpack=True)
fileout = root_original + file[:-4] + "_forNAOJ.txt"
#str = "GalID \t z_spec \t pos0 \t pos_1 \t M_Halo \t M_Stellar \t mag_SDSS_iorig \t mag_SDSS_i \t photoz"
str = "GalID \t z_spec \t pos0 \t pos_1 \t mag_SDSS_i \t photoz"
#dataout = np.c_[data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[11],data_bpz]
#np.savetxt(fileout,dataout,header=str,fmt='%d \t %.3f \t %.7f \t %.7f \t %.3e \t %.3e \t %.2f \t %.2f \t %.2f')
dataout = np.c_[data[0],data[1],data[2],data[3],data[4],data_bpz]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.3f \t %.7f \t %.7f \t %.2f %.2f')

os.system("rm %s" % (file1[:-4]+".bpz"))
os.system("rm %s" % (file1[:-4]+".bpz.bak"))
os.system("rm %s" % (file1[:-4]+"_bpz.cat"))
os.system("rm %s" % (file1[:-4]+".flux_comparison"))
os.system("rm %s" % (file1[:-4]+".probs"))
os.system("rm %s" % file1)
if "gri" in file: os.system("rm %s" % (file1[:-4]+'.columns'))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'
