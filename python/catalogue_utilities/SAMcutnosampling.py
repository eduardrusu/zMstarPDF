# Selects desired fields from the SAM galaxy catalogues and implements cuts

import numpy as np

#sam = "SA"
sam = "Henriques"
z_s = 2.375 # DES0408
lim_i = 22.5

if sam == "SA":
  dirin = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/original/"
  dirout = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/SA_gal_i225_redshift2375/"
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

if sam == "Henriques":
    import class_Henriques2014 # this is class_Henriques2014.py
    dirin = "/lfs08/rusucs/0408/completegalcat/"
    dirout = "/lfs08/rusucs/0408/completegalcat/Henriques_gal_i225_redshift2375/"
    for i in range(8):
        for j in range(8):
            print i, j
            pl = np.linspace(27,63,63 - 27 + 1) # plane 27 redshift 3.06
            out = np.empty(5)
            for k in range(len(pl)):
                #print int(pl[k])
                with open("%sGGL_los_8_%d_%d_N_4096_ang_4_Henriques2014_galaxies_on_plane_%d_f.images" % (dirin,i,j,pl[k]), mode = 'rb') as file:
                    lower_bound = np.fromfile(file, 'f8', 2)
                    upper_bound = np.fromfile(file, 'f8', 2)
                    plane_angle, = np.fromfile(file, 'f8', 1)
                    redshift, = np.fromfile(file, 'f8', 1)
                    n_galaxies, = np.fromfile(file, 'i8', 1)
                    n_cells = np.fromfile(file, 'i4', 2)
                    gal_struct = class_Henriques2014.Henriques2014()
                    galaxy = np.fromfile(file, gal_struct.galaxy_struct, n_galaxies)

                    id = galaxy['galaxy_id']
                    z = galaxy['redshift']
                    pos0 = galaxy['position'][:,0]
                    pos1 = galaxy['position'][:,1]
                    imag = galaxy['mag'][:,gal_struct.filter_number_for_i_band_trans]
                    gal = np.c_[id,z,pos0,pos1,imag].T
                    ind_z = 1
                    ind_imag = 4
                    gal = np.delete(gal,np.where(gal[ind_z] > z_s),axis=1)
                    gal = np.delete(gal,np.where(gal[ind_imag] > lim_i),axis=1)
                    gal = np.delete(gal,np.where(gal[ind_imag] < 0),axis=1) # because there are some -inf in 8_7_7 and 8_3_7
                out = np.c_[out,gal]
            np.savetxt("%sGGL_los_8_%d_%d_N_4096_ang_4_Henriques2014_galaxies.txt" % (dirout,i,j),out.T,fmt='%d %1.7f %1.7f %1.7f %1.2f', header = 'galaxy_id redshift	pos_0	pos_1	i')
