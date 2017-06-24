# The code computes stellar masses with LePhare for the galaxies in the Millennium Simulation, using a set of griK filters. It does this for the photometric redshift estimated in a previous run with photozMillenium_WFI2033.py, and for the "true" catalogue redshift as well.
# Run on SLAC batch as python /scratch/cerusu/MstarMilleniumSLAC_WFI2033.py /scratch/cerusu/GGL_los_8_7_7_3_3_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_WFI2033.images_forlephare.txt
# Needed to remove all headers because SLAC is running an old version of numpy

import numpy as np
import scipy
import sys
import os
from os import system
import time

start_timefield = time.time()

file = str(sys.argv[1])
if "griK" in file:
    fileinspecz = file[:-4] + "_mstarspecz.txt"
    fileinphotoz = file[:-4] + "_mstarphotoz.txt"
fileout = file[:-4] + "_mstar.txt"

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

if "griK" in file:
    specz = 17
    photoz = 18

print file
data = np.loadtxt(file,unpack=True)
if "griK" in file:
    #str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t Ks \t Ks_err \t context \t z-spec \t string"
    dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[specz]]
    #np.savetxt(fileinspecz,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 255 \t %.3f')
    np.savetxt(fileinspecz,dataout,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 255 \t %.3f')
    os.system("/scratch/cerusu/runlephare_phys_para_millenium.sh %s" % fileinspecz) # run LePhare
    os.system("rm %s" % fileinspecz[16:]) # to remove /scratch/cerusu/
    mstar_specz = np.loadtxt(fileinspecz + ".MAG_BC03_I09.lephareout",usecols=[33,35],unpack=True)
    mstar_specz[0][mstar_specz[1] > 0] = mstar_specz[1][mstar_specz[1] > 0] # replace good MASS_BEST with good MASS_MED
    mstar_specz[0][mstar_specz[0] < 0] = 9 # replace bad MASS_BEST with 9
    os.system("rm %s" % fileinspecz[16:] + ".MAG_BC03_I09.lephareout")
    
    #str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t Ks \t Ks_err \t context \t z-spec \t string"
    dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[photoz]]
    #np.savetxt(fileinphotoz,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 255 \t %.2f')
    np.savetxt(fileinphotoz,dataout,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 255 \t %.2f')
    os.system("/scratch/cerusu/runlephare_phys_para_millenium.sh %s" % fileinphotoz)
    os.system("rm %s" % fileinphotoz[16:])
    mstar_photoz = np.loadtxt(fileinphotoz + ".MAG_BC03_I09.lephareout",usecols=[33,35],unpack=True)
    mstar_photoz[0][mstar_photoz[1] > 0] = mstar_photoz[1][mstar_photoz[1] > 0]
    mstar_photoz[0][mstar_photoz[0] < 0] = 9
    #str = "ID \t mstar_specz \t mstar_photoz"
    #np.savetxt(fileout,np.c_[data[id],mstar_specz[0],mstar_photoz[0]],header=str,fmt='%d \t %.3f \t %.3f')
    np.savetxt(fileout,np.c_[data[id],mstar_specz[0],mstar_photoz[0]],fmt='%d \t %.3f \t %.3f')
    os.system("rm %s" % fileinphotoz[16:] + ".MAG_BC03_I09.lephareout")

os.system("rm %s" % file[16:])

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))
                               
print 'Done!'
