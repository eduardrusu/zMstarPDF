# This code samples from the Millenium Simulation (MS) photometry, assuming observed CFHTLenS-like or WFI2033-like uncertainties.The files it creates will be used by photozMillenium_WFI2033.py
# run as: python extractMillenium.py /Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt

import numpy as np
import sys

filein = str(sys.argv[1])
fileout_ugriz = str(sys.argv[1])[:-11] + "_ugriz.images.txt"
fileout_ugrizJHK = str(sys.argv[1])[:-11] + "_ugrizJHK_WFI2033.images.txt"

id = 0
zspec = 5
posx = 6
posy = 7
mhalo = 9
mstar = 11
u = 12
g = 13
r = 14
i = 15
z = 16
J = 17
H = 18
K = 19

data = np.loadtxt(str(sys.argv[1]),usecols=[id,zspec,posx,posy,mhalo,mstar,u,g,r,i,z,J,H,K],comments='GalID',unpack=True)
#data = np.loadtxt("GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt",usecols=[id,zspec,posx,posy,mhalo,mstar,u,g,r,i,z,J,H,K],comments='GalID',unpack=True)

# reposition labels
id = 0
zspec = 1
posx = 2
posy = 3
mhalo = 4
mstar = 5
u = 6
g = 7
r = 8
i = 9
z = 10
J = 11
H = 12
K = 13

def positive(x):
    x[x<0.01] = 0.01
    return x

dataout = np.zeros([np.shape(data)[0]+9,np.shape(data)[1]])

# name dataout labels
id_o = 0
zspec_o = 1
posx_o = 2
posy_o = 3
mhalo_o = 4
mstar_o = 5
iorig_o = 6
u_o = 7
uerr_o = 8
g_o = 9
gerr_o = 10
r_o = 11
rerr_o = 12
i_o = 13
ierr_o = 14
z_o = 15
zerr_o = 16
J_o = 17
Jerr_o = 18
H_o = 19
Herr_o = 20
K_o = 21
Kerr_o = 22

# o stands for out
dataout[id_o] = np.copy(data[id])
dataout[zspec_o] = np.copy(data[zspec])
dataout[posx_o] = np.copy(data[posx])
dataout[posy_o] = np.copy(data[posy])
dataout[mhalo_o] = np.copy(data[mhalo])
dataout[mstar_o] = np.copy(data[mstar])
dataout[iorig_o] = np.copy(data[i])
dataout[u_o] = np.copy(data[u])
dataout[g_o] = np.copy(data[g])
dataout[r_o] = np.copy(data[r])
dataout[i_o] = np.copy(data[i])
dataout[z_o] = np.copy(data[z])
dataout[J_o] = np.copy(data[J])
dataout[H_o] = np.copy(data[H])
dataout[K_o] = np.copy(data[K])

''' Sampling the photometry by assuming error bars from the observations. First assign the non-fixed error bar, then sample a new photometry point'''
# ugriz (CFHTLenS-like)
dataout[uerr_o][data[u]<23]                        = positive(np.random.normal(0.0076, 0.0021, len(data[u][data[u]<23])))
dataout[uerr_o][(data[u]>=23) & (data[u]<24)]      = positive(np.random.normal(0.0136, 0.0046, len(data[u][(data[u]>=23) & (data[u]<24)])))
dataout[uerr_o][(data[u]>=24) & (data[u]<24.5)]    = positive(np.random.normal(0.0212, 0.0072, len(data[u][(data[u]>=24) & (data[u]<24.5)])))
dataout[uerr_o][(data[u]>=24.5) & (data[u]<25)]    = positive(np.random.normal(0.0308, 0.0113, len(data[u][(data[u]>=24.5) & (data[u]<25)])))
dataout[uerr_o][(data[u]>=25) & (data[u]<25.5)]    = positive(np.random.normal(0.0451, 0.0166, len(data[u][(data[u]>=25) & (data[u]<25.5)])))
dataout[uerr_o][(data[u]>=25.5) & (data[u]<26)]    = positive(np.random.normal(0.0689, 0.0230, len(data[u][(data[u]>=25.5) & (data[u]<26)])))
dataout[uerr_o][(data[u]>=26) & (data[u]<26.5)]    = positive(np.random.normal(0.1054, 0.0308, len(data[u][(data[u]>=26) & (data[u]<26.5)])))
dataout[uerr_o][(data[u]>=26.5) & (data[u]<26.75)] = positive(np.random.normal(0.1430, 0.0357, len(data[u][(data[u]>=26.5) & (data[u]<26.75)])))
dataout[uerr_o][(data[u]>=26.75) & (data[u]<27)]   = positive(np.random.normal(0.1735, 0.0413, len(data[u][(data[u]>=26.75) & (data[u]<27)])))
dataout[uerr_o][(data[u]>=27) & (data[u]<27.25)]   = positive(np.random.normal(0.2100, 0.0498, len(data[u][(data[u]>=27) & (data[u]<27.25)])))
dataout[uerr_o][(data[u]>=27.25) & (data[u]<27.5)] = positive(np.random.normal(0.2576, 0.0584, len(data[u][(data[u]>=27.25) & (data[u]<27.5)])))
dataout[uerr_o][(data[u]>=27.5) & (data[u]<27.75)] = positive(np.random.normal(0.3070, 0.0661, len(data[u][(data[u]>=27.5) & (data[u]<27.75)])))
dataout[uerr_o][data[u]>=27.75]                    = positive(np.random.normal(0.4078, 0.0959, len(data[u][data[u]>=27.75])))

dataout[gerr_o][data[g]<24]                        = positive(np.random.normal(0.0098, 0.0038, len(data[g][data[g]<24])))
dataout[gerr_o][(data[g]>=24) & (data[g]<24.5)]    = positive(np.random.normal(0.0156, 0.0062, len(data[g][(data[g]>=24) & (data[g]<24.5)])))
dataout[gerr_o][(data[g]>=24.5) & (data[g]<25)]    = positive(np.random.normal(0.0227, 0.0100, len(data[g][(data[g]>=24.5) & (data[g]<25)])))
dataout[gerr_o][(data[g]>=25) & (data[g]<25.5)]    = positive(np.random.normal(0.0339, 0.0141, len(data[g][(data[g]>=25) & (data[g]<25.5)])))
dataout[gerr_o][(data[g]>=25.5) & (data[g]<26)]    = positive(np.random.normal(0.0514, 0.0213, len(data[g][(data[g]>=25.5) & (data[g]<26)])))
dataout[gerr_o][(data[g]>=26) & (data[g]<26.5)]    = positive(np.random.normal(0.0778, 0.0237, len(data[g][(data[g]>=26) & (data[g]<26.5)])))
dataout[gerr_o][(data[g]>=26.5) & (data[g]<26.75)] = positive(np.random.normal(0.1080, 0.0300, len(data[g][(data[g]>=26.5) & (data[g]<26.75)])))
dataout[gerr_o][(data[g]>=26.75) & (data[g]<27)]   = positive(np.random.normal(0.1329, 0.0308, len(data[g][(data[g]>=26.75) & (data[g]<27)])))
dataout[gerr_o][(data[g]>=27) & (data[g]<27.25)]   = positive(np.random.normal(0.1548, 0.0410, len(data[g][(data[g]>=27) & (data[g]<27.25)])))
dataout[gerr_o][(data[g]>=27.25) & (data[g]<27.5)] = positive(np.random.normal(0.1881, 0.0327, len(data[g][(data[g]>=27.25) & (data[g]<27.5)])))
dataout[gerr_o][(data[g]>=27.5) & (data[g]<27.75)] = positive(np.random.normal(0.2331, 0.0229, len(data[g][(data[g]>=27.5) & (data[g]<27.75)])))
dataout[gerr_o][data[g]>=27.75]                    = positive(np.random.normal(0.2817, 0.1344, len(data[g][data[g]>=27.75])))

dataout[rerr_o][data[r]<23]                        = positive(np.random.normal(0.0064, 0.0022, len(data[r][data[r]<23])))
dataout[rerr_o][(data[r]>=23) & (data[r]<24)]      = positive(np.random.normal(0.0126, 0.0055, len(data[r][(data[r]>=23) & (data[r]<24)])))
dataout[rerr_o][(data[r]>=24) & (data[r]<24.5)]    = positive(np.random.normal(0.0203, 0.0084, len(data[r][(data[r]>=24) & (data[r]<24.5)])))
dataout[rerr_o][(data[r]>=24.5) & (data[r]<25)]    = positive(np.random.normal(0.0293, 0.0147, len(data[r][(data[r]>=24.5) & (data[r]<25)])))
dataout[rerr_o][(data[r]>=25) & (data[r]<25.5)]    = positive(np.random.normal(0.0419, 0.0147, len(data[r][(data[r]>=25) & (data[r]<25.5)])))
dataout[rerr_o][(data[r]>=25.5) & (data[r]<26)]    = positive(np.random.normal(0.0647, 0.0255, len(data[r][(data[r]>=25.5) & (data[r]<26)])))
dataout[rerr_o][(data[r]>=26) & (data[r]<26.5)]    = positive(np.random.normal(0.0834, 0.0196, len(data[r][(data[r]>=26) & (data[r]<26.5)])))
dataout[rerr_o][(data[r]>=26.5) & (data[r]<26.75)] = positive(np.random.normal(0.1261, 0.0277, len(data[r][(data[r]>=26.5) & (data[r]<26.75)])))
dataout[rerr_o][(data[r]>=26.75) & (data[r]<27)]   = positive(np.random.normal(0.1403, 0.0134, len(data[r][(data[r]>=26.75) & (data[r]<27)])))
dataout[rerr_o][data[r]>=27]                       = positive(np.random.normal(0.1943, 0.0152, len(data[r][data[r]>=27])))

dataout[ierr_o][data[i]<23]                        = positive(np.random.normal(0.0049, 0.0030, len(data[i][data[i]<23])))
dataout[ierr_o][(data[i]>=23) & (data[i]<23.5)]    = positive(np.random.normal(0.0115, 0.0051, len(data[i][(data[i]>=23) & (data[i]<23.5)])))
dataout[ierr_o][(data[i]>=23.5) & (data[i]<23.75)] = positive(np.random.normal(0.0155, 0.0066, len(data[i][(data[i]>=23.5) & (data[i]<23.75)])))
dataout[ierr_o][data[i]>=23.75]                    = positive(np.random.normal(0.0193, 0.0078, len(data[i][data[i]>=23.75])))

dataout[zerr_o][data[z]<21]                        = positive(np.random.normal(0.0056, 0.0026, len(data[z][data[z]<21])))
dataout[zerr_o][(data[z]>=21) & (data[z]<21.5)]    = positive(np.random.normal(0.0102, 0.0032, len(data[z][(data[z]>=21) & (data[z]<21.5)])))
dataout[zerr_o][(data[z]>=21.5) & (data[z]<22)]    = positive(np.random.normal(0.0146, 0.0043, len(data[z][(data[z]>=21.5) & (data[z]<22)])))
dataout[zerr_o][(data[z]>=22) & (data[z]<22.5)]    = positive(np.random.normal(0.0212, 0.0084, len(data[z][(data[z]>=22) & (data[z]<22.5)])))
dataout[zerr_o][(data[z]>=22.5) & (data[z]<23)]    = positive(np.random.normal(0.0315, 0.0126, len(data[z][(data[z]>=22.5) & (data[z]<23)])))
dataout[zerr_o][(data[z]>=23) & (data[z]<23.5)]    = positive(np.random.normal(0.0480, 0.0182, len(data[z][(data[z]>=23) & (data[z]<23.5)])))
dataout[zerr_o][(data[z]>=23.5) & (data[z]<24)]    = positive(np.random.normal(0.0721, 0.0244, len(data[z][(data[z]>=23.5) & (data[z]<24)])))
dataout[zerr_o][(data[z]>=24) & (data[z]<24.25)]   = positive(np.random.normal(0.1099, 0.0304, len(data[z][(data[z]>=24) & (data[z]<24.25)])))
dataout[zerr_o][(data[z]>=24.25) & (data[z]<24.5)] = positive(np.random.normal(0.1412, 0.0359, len(data[z][(data[z]>=24.25) & (data[z]<24.5)])))
dataout[zerr_o][(data[z]>=24.5) & (data[z]<24.75)] = positive(np.random.normal(0.1714, 0.0367, len(data[z][(data[z]>=24.5) & (data[z]<24.75)])))
dataout[zerr_o][(data[z]>=24.75) & (data[z]<25)]   = positive(np.random.normal(0.2162, 0.0354, len(data[z][(data[z]>=24.75) & (data[z]<25)])))
dataout[zerr_o][data[z]>=25]                       = positive(np.random.normal(0.2525, 0.0346, len(data[z][data[z]>=25])))

dataout[u_o] = np.random.normal(dataout[u_o], dataout[uerr_o])
dataout[g_o] = np.random.normal(dataout[g_o], dataout[gerr_o])
dataout[r_o] = np.random.normal(dataout[r_o], dataout[rerr_o])
dataout[i_o] = np.random.normal(dataout[i_o], dataout[ierr_o])
dataout[z_o] = np.random.normal(dataout[z_o], dataout[zerr_o])

# eliminate objects that are i> 24 in both the original iorig and the randomized i
dataout = np.delete(dataout,np.where((dataout[i_o] > 24) & (dataout[iorig_o] > 24)),axis=1)

# applying 1-sigma detection limits from Erben et al. 2013, including its uncertainty
dataout[u_o][dataout[u_o]>=np.random.normal(25.24+1.75, 0.17, len(dataout[u_o]))] = 99.0
dataout[g_o][dataout[g_o]>=np.random.normal(25.58+1.75, 0.15, len(dataout[g_o]))] = 99.0
dataout[r_o][dataout[r_o]>=np.random.normal(24.88+1.75, 0.16, len(dataout[r_o]))] = 99.0
dataout[z_o][dataout[z_o]>=np.random.normal(23.46+1.75, 0.20, len(dataout[z_o]))] = 99.0
dataout[uerr_o][dataout[u_o]==99.0] = 25.24+1.75
dataout[gerr_o][dataout[g_o]==99.0] = 25.58+1.75
dataout[rerr_o][dataout[r_o]==99.0] = 24.88+1.75
dataout[zerr_o][dataout[z_o]==99.0] = 23.46+1.75

#head = "GalID \t z_spec \t pos_0[rad] \t pos_1[rad] \t M_Halo[M_sol/h] \t M_Stellar[M_sol/h] \t mag_SDSS_iorig \t mag_SDSS_u \t mag_SDSS_uerr \t mag_SDSS_g \t mag_SDSS_gerr \t mag_SDSS_r \t mag_SDSS_rerr \t mag_SDSS_i \t mag_SDSS_ierr \t mag_SDSS_z \t mag_SDSS_zerr" NOT USING THIS ONE BECAUSE THE SPECIAL CHARACTERS ARE NOT RECOGNIZED PROPERLY BY BPZ
head = "GalID \t z_spec \t pos0 \t pos1 \t M_Halo \t M_Stellar \t mag_SDSS_iorig \t mag_SDSS_u \t mag_SDSS_uerr \t mag_SDSS_g \t mag_SDSS_gerr \t mag_SDSS_r \t mag_SDSS_rerr \t mag_SDSS_i \t mag_SDSS_ierr \t mag_SDSS_z \t mag_SDSS_zerr"
np.savetxt(fileout_ugriz,np.c_[dataout[id_o],dataout[zspec_o],dataout[posx_o],dataout[posy_o],dataout[mhalo_o],dataout[mstar_o],dataout[iorig_o],dataout[u_o],dataout[uerr_o],dataout[g_o],dataout[gerr_o],dataout[r_o],dataout[rerr_o],dataout[i_o],dataout[ierr_o],dataout[z_o],dataout[zerr_o]],header=head,fmt='%d \t %.3f \t %.7f \t %.7f \t %.3e \t %.3e \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f')

dataout = np.zeros([np.shape(data)[0]+9,np.shape(data)[1]])

dataout[id_o] = np.copy(data[id])
dataout[zspec_o] = np.copy(data[zspec])
dataout[posx_o] = np.copy(data[posx])
dataout[posy_o] = np.copy(data[posy])
dataout[mhalo_o] = np.copy(data[mhalo])
dataout[mstar_o] = np.copy(data[mstar])
dataout[iorig_o] = np.copy(data[i])
dataout[u_o] = np.copy(data[u])
dataout[g_o] = np.copy(data[g])
dataout[r_o] = np.copy(data[r])
dataout[i_o] = np.copy(data[i])
dataout[z_o] = np.copy(data[z])
dataout[J_o] = np.copy(data[J])
dataout[H_o] = np.copy(data[H])
dataout[K_o] = np.copy(data[K])

''' Sampling the photometry by assuming error bars from the observations. '''
# ugriJHK (WFI2033-like)
dataout[uerr_o][data[u]<20]                        = positive(np.random.normal(0.01, 0.006, len(data[u][data[u]<20])))
dataout[uerr_o][(data[u]>=20) & (data[u]<20.5)]    = positive(np.random.normal(0.02, 0.009, len(data[u][(data[u]>=20) & (data[u]<20.5)])))
dataout[uerr_o][(data[u]>=20.5) & (data[u]<21)]    = positive(np.random.normal(0.03, 0.008, len(data[u][(data[u]>=20.5) & (data[u]<21)])))
dataout[uerr_o][(data[u]>=21) & (data[u]<21.5)]    = positive(np.random.normal(0.04, 0.013, len(data[u][(data[u]>=21) & (data[u]<21.5)])))
dataout[uerr_o][(data[u]>=21.5) & (data[u]<22)]    = positive(np.random.normal(0.06, 0.015, len(data[u][(data[u]>=21.5) & (data[u]<22)])))
dataout[uerr_o][(data[u]>=22) & (data[u]<22.5)]    = positive(np.random.normal(0.07, 0.016, len(data[u][(data[u]>=22) & (data[u]<22.5)])))
dataout[uerr_o][(data[u]>=22.5) & (data[u]<22.75)] = positive(np.random.normal(0.09, 0.014, len(data[u][(data[u]>=22.5) & (data[u]<22.75)])))
dataout[uerr_o][(data[u]>=22.75) & (data[u]<23)]   = positive(np.random.normal(0.11, 0.022, len(data[u][(data[u]>=22.75) & (data[u]<23)])))
dataout[uerr_o][(data[u]>=23) & (data[u]<23.25)]   = positive(np.random.normal(0.12, 0.013, len(data[u][(data[u]>=23) & (data[u]<23.25)])))
dataout[uerr_o][(data[u]>=23.25) & (data[u]<23.5)] = positive(np.random.normal(0.14, 0.018, len(data[u][(data[u]>=23.25) & (data[u]<23.5)])))
dataout[uerr_o][(data[u]>=23.5) & (data[u]<23.75)] = positive(np.random.normal(0.16, 0.034, len(data[u][(data[u]>=23.5) & (data[u]<23.75)])))
dataout[uerr_o][(data[u]>=23.75) & (data[u]<24)]   = positive(np.random.normal(0.17, 0.029, len(data[u][(data[u]>=23.75) & (data[u]<24)])))
dataout[uerr_o][(data[u]>=24) & (data[u]<24.25)]   = positive(np.random.normal(0.17, 0.017, len(data[u][(data[u]>=24) & (data[u]<24.25)])))
dataout[uerr_o][(data[u]>=24.25) & (data[u]<24.5)] = positive(np.random.normal(0.19, 0.026, len(data[u][(data[u]>=24.25) & (data[u]<24.5)])))
dataout[uerr_o][(data[u]>=24.5) & (data[u]<24.75)] = positive(np.random.normal(0.20, 0.046, len(data[u][(data[u]>=24.5) & (data[u]<24.75)])))
dataout[uerr_o][(data[u]>=24.75) & (data[u]<25)]   = positive(np.random.normal(0.26, 0.038, len(data[u][(data[u]>=24.75) & (data[u]<25)])))
dataout[uerr_o][(data[u]>=25) & (data[u]<25.25)]   = positive(np.random.normal(0.32, 0.032, len(data[u][(data[u]>=25) & (data[u]<25.25)])))
dataout[uerr_o][data[u]>=25.25]                    = positive(np.random.normal(0.36, 0.039, len(data[u][data[u]>=25.25])))

dataout[gerr_o][data[g]<20]                        = positive(np.random.normal(0.005,0.004, len(data[g][data[g]<20])))
dataout[gerr_o][(data[g]>=20) & (data[g]<21)]      = positive(np.random.normal(0.01, 0.002, len(data[g][(data[g]>=20) & (data[g]<21)])))
dataout[gerr_o][(data[g]>=21) & (data[g]<21.5)]    = positive(np.random.normal(0.02, 0.011, len(data[g][(data[g]>=21) & (data[g]<21.5)])))
dataout[gerr_o][(data[g]>=21.5) & (data[g]<22)]    = positive(np.random.normal(0.03, 0.007, len(data[g][(data[g]>=21.5) & (data[g]<22)])))
dataout[gerr_o][(data[g]>=22) & (data[g]<22.5)]    = positive(np.random.normal(0.04, 0.010, len(data[g][(data[g]>=22) & (data[g]<22.5)])))
dataout[gerr_o][(data[g]>=22.5) & (data[g]<22.75)] = positive(np.random.normal(0.05, 0.007, len(data[g][(data[g]>=22.5) & (data[g]<22.75)])))
dataout[gerr_o][(data[g]>=22.75) & (data[g]<23)]   = positive(np.random.normal(0.07, 0.013, len(data[g][(data[g]>=22.75) & (data[g]<23)])))
dataout[gerr_o][(data[g]>=23) & (data[g]<23.25)]   = positive(np.random.normal(0.08, 0.019, len(data[g][(data[g]>=23) & (data[g]<23.25)])))
dataout[gerr_o][(data[g]>=23.25) & (data[g]<23.5)] = positive(np.random.normal(0.09, 0.024, len(data[g][(data[g]>=23.25) & (data[g]<23.5)])))
dataout[gerr_o][(data[g]>=23.5) & (data[g]<23.75)] = positive(np.random.normal(0.12, 0.032, len(data[g][(data[g]>=23.5) & (data[g]<23.75)])))
dataout[gerr_o][(data[g]>=23.75) & (data[g]<24)]   = positive(np.random.normal(0.14, 0.026, len(data[g][(data[g]>=23.75) & (data[g]<24)])))
dataout[gerr_o][(data[g]>=24) & (data[g]<24.25)]   = positive(np.random.normal(0.18, 0.040, len(data[g][(data[g]>=24) & (data[g]<24.25)])))
dataout[gerr_o][(data[g]>=24.25) & (data[g]<24.5)] = positive(np.random.normal(0.22, 0.056, len(data[g][(data[g]>=24.25) & (data[g]<24.5)])))
dataout[gerr_o][(data[g]>=24.5) & (data[g]<24.75)] = positive(np.random.normal(0.25, 0.053, len(data[g][(data[g]>=24.5) & (data[g]<24.75)])))
dataout[gerr_o][data[g]>=24.75]                    = positive(np.random.normal(0.29, 0.062, len(data[g][data[g]>=24.75])))

dataout[rerr_o][data[r]<20]                        = positive(np.random.normal(0.005,0.007, len(data[r][data[r]<20])))
dataout[rerr_o][(data[r]>=20) & (data[r]<20.5)]    = positive(np.random.normal(0.01, 0.008, len(data[r][(data[r]>=20) & (data[r]<20.5)])))
dataout[rerr_o][(data[r]>=20.5) & (data[r]<21)]    = positive(np.random.normal(0.02, 0.008, len(data[r][(data[r]>=20.5) & (data[r]<21)])))
dataout[rerr_o][(data[r]>=21) & (data[r]<21.5)]    = positive(np.random.normal(0.03, 0.007, len(data[r][(data[r]>=21) & (data[r]<21.5)])))
dataout[rerr_o][(data[r]>=21.5) & (data[r]<22)]    = positive(np.random.normal(0.03, 0.009, len(data[r][(data[r]>=21.5) & (data[r]<22)])))
dataout[rerr_o][(data[r]>=22) & (data[r]<22.5)]    = positive(np.random.normal(0.05, 0.016, len(data[r][(data[r]>=22) & (data[r]<22.5)])))
dataout[rerr_o][(data[r]>=22.5) & (data[r]<22.75)] = positive(np.random.normal(0.06, 0.009, len(data[r][(data[r]>=22.5) & (data[r]<22.75)])))
dataout[rerr_o][(data[r]>=22.75) & (data[r]<23)]   = positive(np.random.normal(0.08, 0.016, len(data[r][(data[r]>=22.75) & (data[r]<23)])))
dataout[rerr_o][(data[r]>=23) & (data[r]<23.25)]   = positive(np.random.normal(0.08, 0.019, len(data[r][(data[r]>=23) & (data[r]<23.25)])))
dataout[rerr_o][(data[r]>=23.25) & (data[r]<23.5)] = positive(np.random.normal(0.10, 0.018, len(data[r][(data[r]>=23.25) & (data[r]<23.5)])))
dataout[rerr_o][(data[r]>=23.5) & (data[r]<23.75)] = positive(np.random.normal(0.11, 0.018, len(data[r][(data[r]>=23.5) & (data[r]<23.75)])))
dataout[rerr_o][(data[r]>=23.75) & (data[r]<24)]   = positive(np.random.normal(0.14, 0.021, len(data[r][(data[r]>=23.75) & (data[r]<24)])))
dataout[rerr_o][(data[r]>=24) & (data[r]<24.25)]   = positive(np.random.normal(0.15, 0.020, len(data[r][(data[r]>=24) & (data[r]<24.25)])))
dataout[rerr_o][data[r]>=24.25]                    = positive(np.random.normal(0.16, 0.026, len(data[r][data[r]>=24.25])))

dataout[ierr_o][data[i]<19]                        = positive(np.random.normal(0.005,0.002, len(data[i][data[i]<19])))
dataout[ierr_o][(data[i]>=19) & (data[i]<20)]      = positive(np.random.normal(0.01, 0.005, len(data[i][(data[i]>=19) & (data[i]<20)])))
dataout[ierr_o][(data[i]>=20) & (data[i]<20.5)]    = positive(np.random.normal(0.02, 0.005, len(data[i][(data[i]>=20) & (data[i]<20.5)])))
dataout[ierr_o][(data[i]>=20.5) & (data[i]<21)]    = positive(np.random.normal(0.03, 0.011, len(data[i][(data[i]>=20.5) & (data[i]<21)])))
dataout[ierr_o][(data[i]>=21) & (data[i]<21.5)]    = positive(np.random.normal(0.04, 0.015, len(data[i][(data[i]>=21) & (data[i]<21.5)])))
dataout[ierr_o][(data[i]>=21.5) & (data[i]<22)]    = positive(np.random.normal(0.07, 0.019, len(data[i][(data[i]>=21.5) & (data[i]<22)])))
dataout[ierr_o][(data[i]>=22) & (data[i]<22.5)]    = positive(np.random.normal(0.10, 0.023, len(data[i][(data[i]>=22) & (data[i]<22.5)])))
dataout[ierr_o][(data[i]>=22.5) & (data[i]<22.75)] = positive(np.random.normal(0.13, 0.034, len(data[i][(data[i]>=22.5) & (data[i]<22.75)])))
dataout[ierr_o][(data[i]>=22.75) & (data[i]<23)]   = positive(np.random.normal(0.15, 0.038, len(data[i][(data[i]>=22.75) & (data[i]<23)])))
dataout[ierr_o][(data[i]>=23) & (data[i]<23.25)]   = positive(np.random.normal(0.19, 0.039, len(data[i][(data[i]>=23) & (data[i]<23.25)])))
dataout[ierr_o][data[i]>=23.25]                    = positive(np.random.normal(0.23, 0.031, len(data[i][data[i]>=23.25])))

dataout[zerr_o][data[z]<18]                        = positive(np.random.normal(0.005,0.005, len(data[z][data[z]<18])))
dataout[zerr_o][(data[z]>=18) & (data[z]<19)]      = positive(np.random.normal(0.01, 0.009, len(data[z][(data[z]>=18) & (data[z]<19)])))
dataout[zerr_o][(data[z]>=19) & (data[z]<20)]      = positive(np.random.normal(0.02, 0.008, len(data[z][(data[z]>=19) & (data[z]<20)])))
dataout[zerr_o][(data[z]>=20) & (data[z]<20.5)]    = positive(np.random.normal(0.04, 0.015, len(data[z][(data[z]>=20) & (data[z]<20.5)])))
dataout[zerr_o][(data[z]>=20.5) & (data[z]<21)]    = positive(np.random.normal(0.06, 0.014, len(data[z][(data[z]>=20.5) & (data[z]<21)])))
dataout[zerr_o][(data[z]>=21) & (data[z]<21.5)]    = positive(np.random.normal(0.07, 0.017, len(data[z][(data[z]>=21) & (data[z]<21.5)])))
dataout[zerr_o][(data[z]>=21.5) & (data[z]<22)]    = positive(np.random.normal(0.10, 0.027, len(data[z][(data[z]>=21.5) & (data[z]<22)])))
dataout[zerr_o][(data[z]>=22) & (data[z]<22.5)]    = positive(np.random.normal(0.14, 0.027, len(data[z][(data[z]>=22) & (data[z]<22.5)])))
dataout[zerr_o][(data[z]>=22.5) & (data[z]<22.75)] = positive(np.random.normal(0.17, 0.035, len(data[z][(data[z]>=22.5) & (data[z]<22.75)])))
dataout[zerr_o][(data[z]>=22.75) & (data[z]<23)]   = positive(np.random.normal(0.18, 0.034, len(data[z][(data[z]>=22.75) & (data[z]<23)])))
dataout[zerr_o][(data[z]>=23) & (data[z]<23.25)]   = positive(np.random.normal(0.21, 0.037, len(data[z][(data[z]>=23) & (data[z]<23.25)])))
dataout[zerr_o][data[z]>=23.25]                    = positive(np.random.normal(0.30, 0.082, len(data[z][data[z]>=23.25])))

dataout[Jerr_o][data[J]<16]                        = positive(np.random.normal(np.sqrt(0.005**2+0.04**2),0.004, len(data[J][data[J]<16])))
dataout[Jerr_o][(data[J]>=16) & (data[J]<17)]      = positive(np.random.normal(np.sqrt(0.01**2+0.04**2), 0.006, len(data[J][(data[J]>=16) & (data[J]<17)])))
dataout[Jerr_o][(data[J]>=17) & (data[J]<18)]      = positive(np.random.normal(np.sqrt(0.02**2+0.04**2), 0.014, len(data[J][(data[J]>=17) & (data[J]<18)])))
dataout[Jerr_o][(data[J]>=18) & (data[J]<18.5)]    = positive(np.random.normal(np.sqrt(0.03**2+0.04**2), 0.009, len(data[J][(data[J]>=18) & (data[J]<18.5)])))
dataout[Jerr_o][(data[J]>=18.5) & (data[J]<19)]    = positive(np.random.normal(np.sqrt(0.04**2+0.04**2), 0.010, len(data[J][(data[J]>=18.5) & (data[J]<19)])))
dataout[Jerr_o][(data[J]>=19) & (data[J]<19.5)]    = positive(np.random.normal(np.sqrt(0.06**2+0.04**2), 0.020, len(data[J][(data[J]>=19) & (data[J]<19.5)])))
dataout[Jerr_o][(data[J]>=19.5) & (data[J]<20)]    = positive(np.random.normal(np.sqrt(0.07**2+0.04**2), 0.022, len(data[J][(data[J]>=19.5) & (data[J]<20)])))
dataout[Jerr_o][(data[J]>=20) & (data[J]<20.5)]    = positive(np.random.normal(np.sqrt(0.10**2+0.04**2), 0.031, len(data[J][(data[J]>=20) & (data[J]<20.5)])))
dataout[Jerr_o][(data[J]>=20.5) & (data[J]<20.75)] = positive(np.random.normal(np.sqrt(0.11**2+0.04**2), 0.037, len(data[J][(data[J]>=20.5) & (data[J]<20.75)])))
dataout[Jerr_o][(data[J]>=20.75) & (data[J]<21)]   = positive(np.random.normal(np.sqrt(0.15**2+0.04**2), 0.028, len(data[J][(data[J]>=20.75) & (data[J]<21)])))
dataout[Jerr_o][(data[J]>=21) & (data[J]<21.25)]   = positive(np.random.normal(np.sqrt(0.16**2+0.04**2), 0.035, len(data[J][(data[J]>=21) & (data[J]<21.25)])))
dataout[Jerr_o][(data[J]>=21.25) & (data[J]<21.5)] = positive(np.random.normal(np.sqrt(0.18**2+0.04**2), 0.028, len(data[J][(data[J]>=21.25) & (data[J]<21.5)])))
dataout[Jerr_o][data[J]>=21.5]                     = positive(np.random.normal(np.sqrt(0.19**2+0.04**2), 0.029, len(data[J][data[J]>=21.5])))

dataout[Herr_o][data[H]<15]                        = positive(np.random.normal(np.sqrt(0.005**2+0.04**2),0.003, len(data[H][data[H]<15])))
dataout[Herr_o][(data[H]>=15) & (data[H]<16)]      = positive(np.random.normal(np.sqrt(0.01**2+0.04**2), 0.003, len(data[H][(data[H]>=15) & (data[H]<16)])))
dataout[Herr_o][(data[H]>=16) & (data[H]<17)]      = positive(np.random.normal(np.sqrt(0.02**2+0.04**2), 0.013, len(data[H][(data[H]>=16) & (data[H]<17)])))
dataout[Herr_o][(data[H]>=17) & (data[H]<17.5)]    = positive(np.random.normal(np.sqrt(0.03**2+0.04**2), 0.014, len(data[H][(data[H]>=17) & (data[H]<17.5)])))
dataout[Herr_o][(data[H]>=17.5) & (data[H]<18)]    = positive(np.random.normal(np.sqrt(0.04**2+0.04**2), 0.008, len(data[H][(data[H]>=17.5) & (data[H]<18)])))
dataout[Herr_o][(data[H]>=18) & (data[H]<18.5)]    = positive(np.random.normal(np.sqrt(0.06**2+0.04**2), 0.015, len(data[H][(data[H]>=18) & (data[H]<18.5)])))
dataout[Herr_o][(data[H]>=18.5) & (data[H]<19)]    = positive(np.random.normal(np.sqrt(0.07**2+0.04**2), 0.025, len(data[H][(data[H]>=18.5) & (data[H]<19)])))
dataout[Herr_o][(data[H]>=19) & (data[H]<19.5)]    = positive(np.random.normal(np.sqrt(0.09**2+0.04**2), 0.021, len(data[H][(data[H]>=19) & (data[H]<19.5)])))
dataout[Herr_o][(data[H]>=19.5) & (data[H]<19.75)] = positive(np.random.normal(np.sqrt(0.13**2+0.04**2), 0.048, len(data[H][(data[H]>=19.5) & (data[H]<19.75)])))
dataout[Herr_o][(data[H]>=19.75) & (data[H]<20)]   = positive(np.random.normal(np.sqrt(0.12**2+0.04**2), 0.034, len(data[H][(data[H]>=19.75) & (data[H]<20)])))
dataout[Herr_o][(data[H]>=20) & (data[H]<20.25)]   = positive(np.random.normal(np.sqrt(0.17**2+0.04**2), 0.038, len(data[H][(data[H]>=20) & (data[H]<20.25)])))
dataout[Herr_o][(data[H]>=20.25) & (data[H]<20.5)] = positive(np.random.normal(np.sqrt(0.17**2+0.04**2), 0.025, len(data[H][(data[H]>=20.25) & (data[H]<20.5)])))
dataout[Herr_o][data[H]>=20.5]                     = positive(np.random.normal(np.sqrt(0.19**2+0.04**2), 0.029, len(data[H][data[H]>=20.5])))

dataout[Kerr_o][data[K]<15]                        = positive(np.random.normal(np.sqrt(0.005**2+0.08**2),0.004, len(data[K][data[K]<15])))
dataout[Kerr_o][(data[K]>=15) & (data[K]<16)]      = positive(np.random.normal(np.sqrt(0.02**2+0.08**2), 0.013, len(data[K][(data[K]>=15) & (data[K]<16)])))
dataout[Kerr_o][(data[K]>=16) & (data[K]<17)]      = positive(np.random.normal(np.sqrt(0.03**2+0.08**2), 0.007, len(data[K][(data[K]>=16) & (data[K]<17)])))
dataout[Kerr_o][(data[K]>=17) & (data[K]<17.5)]    = positive(np.random.normal(np.sqrt(0.04**2+0.08**2), 0.016, len(data[K][(data[K]>=17) & (data[K]<17.5)])))
dataout[Kerr_o][(data[K]>=17.5) & (data[K]<18)]    = positive(np.random.normal(np.sqrt(0.06**2+0.08**2), 0.018, len(data[K][(data[K]>=17.5) & (data[K]<18)])))
dataout[Kerr_o][(data[K]>=18) & (data[K]<18.5)]    = positive(np.random.normal(np.sqrt(0.07**2+0.08**2), 0.014, len(data[K][(data[K]>=18) & (data[K]<18.5)])))
dataout[Kerr_o][(data[K]>=18.5) & (data[K]<19)]    = positive(np.random.normal(np.sqrt(0.10**2+0.08**2), 0.026, len(data[K][(data[K]>=18.5) & (data[K]<19)])))
dataout[Kerr_o][(data[K]>=19) & (data[K]<19.5)]    = positive(np.random.normal(np.sqrt(0.14**2+0.08**2), 0.032, len(data[K][(data[K]>=19) & (data[K]<19.5)])))
dataout[Kerr_o][(data[K]>=19.5) & (data[K]<19.75)] = positive(np.random.normal(np.sqrt(0.17**2+0.08**2), 0.035, len(data[K][(data[K]>=19.5) & (data[K]<19.75)])))
dataout[Kerr_o][(data[K]>=19.75) & (data[K]<20)]   = positive(np.random.normal(np.sqrt(0.19**2+0.08**2), 0.070, len(data[K][(data[K]>=19.75) & (data[K]<20)])))
dataout[Kerr_o][data[K]>=20]                       = positive(np.random.normal(np.sqrt(0.19**2+0.08**2), 0.037, len(data[K][data[K]>=20])))

dataout = np.delete(dataout,np.where(dataout[i_o] > 23),axis=1)

dataout[u_o] = np.random.normal(dataout[u_o], dataout[uerr_o])
dataout[g_o] = np.random.normal(dataout[g_o], dataout[gerr_o])
dataout[r_o] = np.random.normal(dataout[r_o], dataout[rerr_o])
dataout[i_o] = np.random.normal(dataout[i_o], dataout[ierr_o])
dataout[z_o] = np.random.normal(dataout[z_o], dataout[zerr_o])
dataout[J_o] = np.random.normal(dataout[J_o], dataout[Jerr_o])
dataout[H_o] = np.random.normal(dataout[H_o], dataout[Herr_o])
dataout[K_o] = np.random.normal(dataout[K_o], dataout[Kerr_o])

# eliminate objects that are i> 24 in both the original iorig and the randomized i
dataout = np.delete(dataout,np.where((dataout[i_o] > 23) & (dataout[iorig_o] > 23)),axis=1)

# applying 1-sigma detection limits including their uncertainty; these are final values I used as input for lephare for the real lens field; all values are in AB
dataout[u_o][dataout[u_o]>=np.random.normal(27.39, 0.05, len(dataout[u_o]))] = 99.0
dataout[g_o][dataout[g_o]>=np.random.normal(25.48, 0.15, len(dataout[g_o]))] = 99.0
dataout[r_o][dataout[r_o]>=np.random.normal(25.55, 0.12, len(dataout[r_o]))] = 99.0
dataout[z_o][dataout[z_o]>=np.random.normal(24.64, 0.45, len(dataout[z_o]))] = 99.0
dataout[J_o][dataout[J_o]>=np.random.normal(23.21, 0.12, len(dataout[J_o]))] = 99.0
dataout[H_o][dataout[H_o]>=np.random.normal(22.60, 0.08, len(dataout[H_o]))] = 99.0
dataout[K_o][dataout[K_o]>=np.random.normal(22.50, 0.04, len(dataout[K_o]))] = 99.0
dataout[uerr_o][dataout[u_o]==99.0] = 27.39
dataout[gerr_o][dataout[g_o]==99.0] = 25.48
dataout[rerr_o][dataout[r_o]==99.0] = 25.55
dataout[zerr_o][dataout[z_o]==99.0] = 24.64
dataout[Jerr_o][dataout[J_o]==99.0] = 23.21
dataout[Herr_o][dataout[H_o]==99.0] = 22.60
dataout[Kerr_o][dataout[K_o]==99.0] = 22.50

#head = "GalID \t z_spec \t pos_0[rad] \t pos_1[rad] \t M_Halo[M_sol/h] \t M_Stellar[M_sol/h] \t mag_SDSS_iorig \t mag_SDSS_u \t mag_SDSS_uerr \t mag_SDSS_g \t mag_SDSS_gerr \t mag_SDSS_r \t mag_SDSS_rerr \t mag_SDSS_i \t mag_SDSS_ierr \t mag_SDSS_z \t mag_SDSS_zerr \t mag_J \t mag_Jerr \t mag_H \t mag_Herr \t mag_K \t mag_Kerr" NOT USING THIS ONE BECAUSE THE SPECIAL CHARACTERS ARE NOT RECOGNIZED PROPERLY BY BPZ
head = "GalID \t z_spec \t pos0 \t pos_1 \t M_Halo \t M_Stellar \t mag_SDSS_iorig \t mag_SDSS_u \t mag_SDSS_uerr \t mag_SDSS_g \t mag_SDSS_gerr \t mag_SDSS_r \t mag_SDSS_rerr \t mag_SDSS_i \t mag_SDSS_ierr \t mag_SDSS_z \t mag_SDSS_zerr \t mag_J \t mag_Jerr \t mag_H \t mag_Herr \t mag_K \t mag_Kerr"
np.savetxt(fileout_ugrizJHK,dataout.T,header=head,fmt='%d \t %.3f \t %.7f \t %.7f \t %.3e \t %.3e \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t')

print str(sys.argv[1]) + ' Done!'

