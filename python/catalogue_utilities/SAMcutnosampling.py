# Selects desired fields from the SAM galaxy catalogues and implements cuts

import numpy as np
dirin = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/original/"
dirout = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/SA_gal_i225_redshift2375/"

z_s = 2.375 # DES0408
lim_i = 22.5

for i in range(8):
    for j in range(8):
        for k in range(4):
            for l in range(4):
                print i, j, k, l
                input = dirin + "GGL_los_8_%s_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt" %(i,j,k,l)
                output = dirout + "SA_gal_i225_redshift2375_8_%s_%s_%s_%s_N_4096_ang_4.images.txt" %(i,j,k,l)
                zspec = 5
                pos0 = 6
                pos1 = 7
                mag_g = 13
                mag_r = 14
                mag_i = 15
                mag_z = 16
                data = np.loadtxt(input,comments="GalID",usecols=[zspec,pos0,pos1,mag_g,mag_r,mag_i,mag_z],unpack=True)
                zspec = 0
                pos0 = 1
                pos1 = 2
                mag_g = 3
                mag_r = 4
                mag_i = 5
                mag_z = 6
                dataout1 = data[zspec][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout2 = data[pos0][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout3 = data[pos1][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout4 = data[mag_g][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout5 = data[mag_r][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout6 = data[mag_i][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                dataout7 = data[mag_z][(data[zspec]<=z_s) & (data[mag_i]<=lim_i)]
                np.savetxt(output,np.c_[dataout1,dataout2,dataout3,dataout4,dataout5,dataout6,dataout7],fmt='%1.2f %1.7f %1.7f %1.2f %1.2f %1.2f %1.2f', header = 'z_spec pos_0[rad] pos_1[rad] mag_SDSS_g mag_SDSS_r mag_SDSS_i mag_SDSS_z')
