# The code calculates realistic BPZ photoz for galaxies in the Millenium Simulation, based either on griK fluxes. It then creates the necessary files to run LePhare on SLAC. The code expects input files created by extractMillenium.py
# run as: python photozMillenium_WFI2033_griK.py /Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images.txt

import numpy as np
import sys
import os
from os import system
import time

start_timefield = time.time()

root_bpz = "/Users/perseus/bpz-1.99.3/test/"
root_original = "/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/txt/"
file = str(sys.argv[1])
file1 = file.replace(root_original, root_bpz)
file2 = file1.replace("ugrizJHK", "griK")
os.system("cp %s %s" % (file,file2))
if "ugrizJHK" in file:
    os.system("cp %smillennium_griK.columns %s" % (root_bpz,file2[:-4]+'.columns'))

os.system("python $BPZPATH/bpz.py %s -INTERP 2" % file2)
os.system("python $BPZPATH/bpzfinalize.py %s" % file2[:-4])

id = 0
u = 1
u_err = 2
g = 3
g_err = 4
r = 5
r_err = 6
i = 7
i_err = 8
z = 9
z_err = 10
J = 11
J_err = 12
H = 13
H_err = 14
K = 15
K_err = 16

if "ugrizJHK" in file:
    data = np.loadtxt(file,usecols=[0,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],unpack=True)
    data_bpz = np.loadtxt(file2[:-4]+"_bpz.cat",usecols=[9,1],unpack=True)
    ''' Write the mags into LePhare-expected format, and assuming no observations for mags below detection threshold '''
    u_corr =  0.0 # for griK these are the offsets suggested by running BPZ on one simulation field; scaled so that I apply the same i-band offset as for ugrizJHK
    g_corr =  0.010
    r_corr =  0.010
    i_corr =  0.016
    z_corr =  0.0
    J_corr =  0.0
    H_corr =  0.0
    K_corr = -0.012
    # If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
    data[u_err] = -99.0
    data[g_err][np.abs(data[g]) == 99.0] = -99.0
    data[r_err][np.abs(data[r]) == 99.0] = -99.0
    data[i_err][np.abs(data[i]) == 99.0] = -99.0
    data[z_err] = -99.0
    data[J_err] = -99.0
    data[H_err] = -99.0
    data[K_err][np.abs(data[K]) == 99.0] = -99.0
    data[u] = -99.0
    data[g][np.abs(data[g]) == 99.0] = -99.0
    data[r][np.abs(data[r]) == 99.0] = -99.0
    data[i][np.abs(data[i]) == 99.0] = -99.0
    data[z] = -99.0
    data[J] = -99.0
    data[H] = -99.0
    data[K][np.abs(data[K]) == 99.0] = -99.0
    # apply the corrections suggested by BPZ
    data[u][np.abs(data[u]) != 99.0] += u_corr
    data[g][np.abs(data[g]) != 99.0] += g_corr
    data[r][np.abs(data[r]) != 99.0] += r_corr
    data[i][np.abs(data[i]) != 99.0] += i_corr
    data[z][np.abs(data[z]) != 99.0] += z_corr
    data[J][np.abs(data[J]) != 99.0] += J_corr
    data[H][np.abs(data[H]) != 99.0] += H_corr
    data[K][np.abs(data[K]) != 99.0] += K_corr
    # LePhare thinks error bars = 0 means non-detection, so fix this
    data[u_err][data[u_err] == 0.00] = 0.01
    data[g_err][data[g_err] == 0.00] = 0.01
    data[r_err][data[r_err] == 0.00] = 0.01
    data[i_err][data[i_err] == 0.00] = 0.01
    data[z_err][data[z_err] == 0.00] = 0.01
    data[J_err][data[J_err] == 0.00] = 0.01
    data[H_err][data[H_err] == 0.00] = 0.01
    data[K_err][data[K_err] == 0.00] = 0.01
    fileout = file[:-4] + "_forlephare.txt"
    str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t K \t K_err \t specz \t photoz"
    dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data_bpz[0],data_bpz[1]]
    np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.3f \t %.2f')

os.system("rm %s" % (file2[:-4]+".bpz"))
os.system("rm %s" % (file2[:-4]+".bpz.bak"))
os.system("rm %s" % (file2[:-4]+"_bpz.cat"))
os.system("rm %s" % (file2[:-4]+".flux_comparison"))
os.system("rm %s" % (file2[:-4]+".probs"))
os.system("rm %s" % file2)
if "ugrizJHK" in file:
    os.system("rm %s" % (file2[:-4]+'.columns'))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

