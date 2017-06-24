# This code uses the output from MstarMilleniumSLAC_WFI2033.py, photozMillenium_WFI2033.py and extractMillenium.py in order to prepare a final file from the Millennium Simulation to send to the NAOJ server. It also computes halo masses following Behroozi et al 2010
# Run as python finalprepMilleniumNAOJ_griK.py GGL_los_8_4_6_2_1_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images.txt

import numpy as np
import sys

fileextract = str(sys.argv[1])

#fileextract = "GGL_los_8_4_6_2_1_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images.txt"
#fileforlephare = "GGL_los_8_1_2_3_2_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_WFI2033.images_forlephare.txt"
#filemstar = "GGL_los_8_4_5_3_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_WFI2033.images_forlephare_mstar.txt"
fileforlephare = fileextract[:-27] + "griK_WFI2033.images_forlephare.txt"
filemstar = fileextract[:-27] + "griK_WFI2033.images_forlephare_mstar.txt"
fileout = fileextract[:-27] + "griK_WFI2033.images_forNAOJ.txt"

id = 0
zspec = 1
posx = 2
posy = 3
mhalo = 4
mstar = 5
iorig = 6
i_o = 13

photoz = 18
mstar_specz = 1
mstar_photoz = 2
dataextract = np.loadtxt(fileextract,usecols=[id,zspec,posx,posy,mhalo,mstar,iorig,i_o],unpack=True)
dataforlephare = np.loadtxt(fileforlephare,usecols=[photoz],unpack=True)
datamstar = np.loadtxt(filemstar,usecols=[mstar_specz,mstar_photoz],unpack=True)

id = 0
zspec = 1
posx = 2
posy = 3
mhalo = 4
mstar = 5
iorig = 6
i_o = 7
mstar_specz = 0
mstar_photoz = 1

# Behroozi et al 2010 parameters for z < 1:
M10_ = 12.35
M1a_ = 0.28
Ms00_ = 10.72
Ms0a_ = 0.55
b0_ = 0.44
ba_ = 0.18
d0_ = 0.57
da_ = 0.17
g0_ = 1.56
ga_ = 2.51
# z >= 1:
M10 = 12.27
M1a = -0.84
Ms00 = 11.09
Ms0a = 0.56
b0 = 0.65
ba = 0.31
d0 = 0.56
da = -0.12
g0 = 1.12
ga = -0.53

datahalo = np.zeros([len(dataextract[id]),2])

if "JHK" in fileextract:
    a = 1 / (1 + dataextract[zspec][dataextract[zspec] <= 1])
    logM1a = M10_ + M1a_ * (a - 1)
    logMs0a = Ms00_ + Ms0a_ * (a - 1)
    notlogMs0a = 10 ** logMs0a
    b = b0_ + ba_ * (a - 1)
    d = d0_ + da_ * (a - 1)
    g = g0_ + ga_ * (a - 1)
    datahalo[:,0][dataextract[zspec] <= 1] = logM1a + b * (datamstar[mstar_specz][dataextract[zspec] <= 1] - logMs0a) + ((10 ** datamstar[mstar_specz][dataextract[zspec] <= 1]/notlogMs0a)**d)/(1+(10 ** datamstar[mstar_specz][dataextract[zspec] <= 1]/notlogMs0a)**(-g)) - 1/2
    del logM1a
    del logMs0a
    del notlogMs0a
    del b
    del d
    del g

    a = 1 / (1 + dataextract[zspec][dataextract[zspec] > 1])
    logM1a = M10 + M1a * (a-1)
    logMs0a = Ms00 + Ms0a * (a-1)
    notlogMs0a = 10 ** logMs0a
    b = b0 + ba * (a-1)
    d = d0 + da * (a-1)
    g = g0 + ga * (a-1)
    datahalo[:,0][dataextract[zspec] > 1] = logM1a + b * (datamstar[mstar_specz][dataextract[zspec] > 1] - logMs0a) + ((10 ** datamstar[mstar_specz][dataextract[zspec] > 1]/notlogMs0a)**d)/(1+(10 ** datamstar[mstar_specz][dataextract[zspec] > 1]/notlogMs0a)**(-g)) - 1/2
    del logM1a
    del logMs0a
    del notlogMs0a
    del b
    del d
    del g

    a = 1 / (1 + dataforlephare[dataforlephare <= 1])
    logM1a = M10_ + M1a_ * (a - 1)
    logMs0a = Ms00_ + Ms0a_ * (a - 1)
    notlogMs0a = 10 ** logMs0a
    b = b0_ + ba_ * (a - 1)
    d = d0_ + da_ * (a - 1)
    g = g0_ + ga_ * (a - 1)
    datahalo[:,1][dataforlephare <= 1] = logM1a + b * (datamstar[mstar_photoz][dataforlephare <= 1] - logMs0a) + ((10 ** datamstar[mstar_photoz][dataforlephare <= 1]/notlogMs0a)**d)/(1+(10 ** datamstar[mstar_photoz][dataforlephare <= 1]/notlogMs0a)**(-g)) - 1/2
    del logM1a
    del logMs0a
    del notlogMs0a
    del b
    del d
    del g
    
    a = 1 / (1 + dataforlephare[dataforlephare > 1])
    logM1a = M10 + M1a * (a-1)
    logMs0a = Ms00 + Ms0a * (a-1)
    notlogMs0a = 10 ** logMs0a
    b = b0 + ba * (a-1)
    d = d0 + da * (a-1)
    g = g0 + ga * (a-1)
    datahalo[:,1][dataforlephare > 1] = logM1a + b * (datamstar[mstar_photoz][dataforlephare > 1] - logMs0a) + ((10 ** datamstar[mstar_photoz][dataforlephare > 1]/notlogMs0a)**d)/(1+(10 ** datamstar[mstar_photoz][dataforlephare > 1]/notlogMs0a)**(-g)) - 1/2
    del logM1a
    del logMs0a
    del notlogMs0a
    del b
    del d
    del g

dataout = np.c_[dataextract[id],dataextract[zspec],dataextract[posx],dataextract[posy],dataextract[mhalo],dataextract[mstar],dataextract[iorig],dataextract[i_o],dataforlephare,datamstar[mstar_specz],datamstar[mstar_photoz],datahalo[:,0],datahalo[:,1]]
head = "GalID \t z_spec \t pos_0 \t pos_1 \t M_Halo \t M_Stellar \t mag_SDSS_iorig \t mag_SDSS_i \t photoz \t mstar_specz \t mstar_photoz \t mhalo_specz \t mhalo_photoz"
np.savetxt(fileout,dataout,header=head,fmt='%d \t %.3f \t %.7f \t %.7f \t %.3e \t %.3e \t %.2f \t %.2f \t %.2f \t %.3f \t %.3f \t %.3f \t %.3f')

print str(sys.argv[1]) + ' Done!'
