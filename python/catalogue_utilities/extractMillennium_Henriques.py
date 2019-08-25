# This code samples from the Millenium Simulation (MS) photometry, assuming observed CFHTLenS-like or observed lens-like uncertainties. The files it creates will be used by photozMillenium.py
# run as: python extractMillennium_Henriques.py GGL_los_8_0_0_N_4096_ang_4_Henriques2014_galaxies_on_plane

import numpy as np
import sys
import class_Henriques2014 # this is class_Henriques2014.py

ilim = 22.5
filters = "griz"
rootin = "/lfs08/rusucs/0408/completegalcat/"
rootout = "/lfs08/rusucs/0408/completegalcat/Henriques_gal_i225_sampled/"
filein = str(sys.argv[1])
fileout = rootout + filein + "_%s.images.txt" % filters

pl = np.linspace(27,63,63 - 27 + 1) # plane 27 redshift 3.06
data = np.empty(8)
for i in range(len(pl)):
    with open("%s%s_%d_f.images" % (rootin,filein,pl[i]), mode = 'rb') as file:
        lower_bound = np.fromfile(file, 'f8', 2)
        upper_bound = np.fromfile(file, 'f8', 2)
        plane_angle, = np.fromfile(file, 'f8', 1)
        redshift, = np.fromfile(file, 'f8', 1)
        n_galaxies, = np.fromfile(file, 'i8', 1)
        n_cells = np.fromfile(file, 'i4', 2)
        gal_struct = class_Henriques2014.Henriques2014()
        galaxy = np.fromfile(file, gal_struct.galaxy_struct, n_galaxies)

        id = galaxy['galaxy_id']
        zspec = galaxy['redshift']
        pos0 = galaxy['position'][:,0]
        pos1 = galaxy['position'][:,1]
        gmag = galaxy['mag'][:,gal_struct.filter_number_for_g_band_trans]
        rmag = galaxy['mag'][:,gal_struct.filter_number_for_r_band_trans]
        imag = galaxy['mag'][:,gal_struct.filter_number_for_i_band_trans]
        zmag = galaxy['mag'][:,gal_struct.filter_number_for_z_band_trans]
        gal = np.c_[id,zspec,pos0,pos1,gmag,rmag,imag,zmag].T
    data = np.c_[data,gal]

id = 0
zspec = 1
posx = 2
posy = 3
g = 4
r = 5
i = 6
z = 7
dataout = np.zeros([np.shape(data)[0]+4,np.shape(data)[1]]) # use this instead of dataout since I am modifying the content but I will still need the original content when working with multiple sets of filters
id_o = 0
zspec_o = 1
posx_o = 2
posy_o = 3
g_o = 4
gerr_o = 5
r_o = 6
rerr_o = 7
i_o = 8
ierr_o = 9
z_o = 10
zerr_o = 11

dataout[id_o] = np.copy(data[id])
dataout[zspec_o] = np.copy(data[zspec])
dataout[posx_o] = np.copy(data[posx])
dataout[posy_o] = np.copy(data[posy])
dataout[g_o] = np.copy(data[g])
dataout[r_o] = np.copy(data[r])
dataout[i_o] = np.copy(data[i])
dataout[z_o] = np.copy(data[z])

''' Sampling the photometry by assuming error bars from the observations. First assign the non-fixed error bar, then sample a new photometry point'''
if filters == "griz":
# griz (DES-like)
    file_errorbars = "/lfs08/rusucs/code/median_errors_hlin_12Aug2019edited.txt"
    err = np.loadtxt(file_errorbars,unpack=True)
    err_maginf = 0
    err_magsup = 1
    err_g_lnmed = 2
    err_g_lnsig = 3
    err_r_lnmed = 5
    err_r_lnsig = 6
    err_i_lnmed = 8
    err_i_lnsig = 9
    err_z_lnmed = 11
    err_z_lnsig = 12

    dataout[gerr_o][dataout[g_o] < err[err_maginf][0]]                                            = 2.718 ** np.random.normal(err[err_g_lnmed][0], err[err_g_lnsig][0], len(dataout[gerr_o][dataout[g_o] < err[err_maginf][0]]))
    for j in range(len(err[0])):
        dataout[gerr_o][(dataout[g_o] >= err[err_maginf][j]) & (dataout[g_o] < err[err_magsup][j])] = 2.718 ** np.random.normal(err[err_g_lnmed][j], err[err_g_lnsig][j], len(dataout[gerr_o][(dataout[g_o] >= err[err_maginf][j]) & (dataout[g_o] < err[err_magsup][j])]))
    dataout[gerr_o][dataout[g_o] >= err[err_magsup][-1]]                                          = 2.718 ** np.random.normal(err[err_g_lnmed][-1], err[err_g_lnsig][-1], len(dataout[gerr_o][dataout[g_o] >= err[err_magsup][-1]]))

    dataout[rerr_o][dataout[r_o] < err[err_maginf][0]]                                            = 2.718 ** np.random.normal(err[err_r_lnmed][0], err[err_r_lnsig][0], len(dataout[rerr_o][dataout[r_o] < err[err_maginf][0]]))
    for j in range(len(err[0])):
        dataout[rerr_o][(dataout[r_o] >= err[err_maginf][j]) & (dataout[r_o] < err[err_magsup][j])] = 2.718 ** np.random.normal(err[err_r_lnmed][j], err[err_r_lnsig][j], len(dataout[rerr_o][(dataout[r_o] >= err[err_maginf][j]) & (dataout[r_o] < err[err_magsup][j])]))
    dataout[rerr_o][dataout[r_o] >= err[err_magsup][-1]]                                          = 2.718 ** np.random.normal(err[err_r_lnmed][-1], err[err_r_lnsig][-1], len(dataout[rerr_o][dataout[r_o] >= err[err_magsup][-1]]))

    dataout[ierr_o][dataout[i_o] < err[err_maginf][0]]                                            = 2.718 ** np.random.normal(err[err_i_lnmed][0], err[err_i_lnsig][0], len(dataout[ierr_o][dataout[i_o] < err[err_maginf][0]]))
    for j in range(len(err[0])):
        dataout[ierr_o][(dataout[i_o] >= err[err_maginf][j]) & (dataout[i_o] < err[err_magsup][j])] = 2.718 ** np.random.normal(err[err_i_lnmed][j], err[err_i_lnsig][j], len(dataout[ierr_o][(dataout[i_o] >= err[err_maginf][j]) & (dataout[i_o] < err[err_magsup][j])]))
    dataout[ierr_o][dataout[i_o] >= err[err_magsup][-1]]                                          = 2.718 ** np.random.normal(err[err_i_lnmed][-1], err[err_i_lnsig][-1], len(dataout[ierr_o][dataout[i_o] >= err[err_magsup][-1]]))

    dataout[zerr_o][dataout[z_o] < err[err_maginf][0]]                                            = 2.718 ** np.random.normal(err[err_z_lnmed][0], err[err_z_lnsig][0], len(dataout[zerr_o][dataout[z_o] < err[err_maginf][0]]))
    for j in range(len(err[0])):
        dataout[zerr_o][(dataout[z_o] >= err[err_maginf][j]) & (dataout[z_o] < err[err_magsup][j])] = 2.718 ** np.random.normal(err[err_z_lnmed][j], err[err_z_lnsig][j], len(dataout[zerr_o][(dataout[z_o] >= err[err_maginf][j]) & (dataout[z_o] < err[err_magsup][j])]))
    dataout[zerr_o][dataout[z_o] >= err[err_magsup][-1]]                                          = 2.718 ** np.random.normal(err[err_z_lnmed][-1], err[err_z_lnsig][-1], len(dataout[zerr_o][dataout[z_o] >= err[err_magsup][-1]]))

    dataout[g_o] = np.random.normal(dataout[g_o], dataout[gerr_o])
    dataout[r_o] = np.random.normal(dataout[r_o], dataout[rerr_o])
    dataout[i_o] = np.random.normal(dataout[i_o], dataout[ierr_o])
    dataout[z_o] = np.random.normal(dataout[z_o], dataout[zerr_o])

    dataout[gerr_o][dataout[gerr_o] <= 0.01] = 0.01
    dataout[rerr_o][dataout[rerr_o] <= 0.01] = 0.01
    dataout[ierr_o][dataout[ierr_o] <= 0.01] = 0.01
    dataout[zerr_o][dataout[zerr_o] <= 0.01] = 0.01

    dataout = np.delete(dataout,np.where(dataout[i_o] > ilim),axis=1)
    dataout = np.delete(dataout,np.where(dataout[i_o] < 0),axis=1) # because there are some -inf in 8_7_7 and 8_3_7

    head = "GalID \t z_spec \t pos0 \t pos1 \t mag_SDSS_g \t mag_SDSS_gerr \t mag_SDSS_r \t mag_SDSS_rerr \t mag_SDSS_i \t mag_SDSS_ierr \t mag_SDSS_z \t mag_SDSS_zerr"
    np.savetxt(fileout,np.c_[dataout[id_o],dataout[zspec_o],dataout[posx_o],dataout[posy_o],dataout[g_o],dataout[gerr_o],dataout[r_o],dataout[rerr_o],dataout[i_o],dataout[ierr_o],dataout[z_o],dataout[zerr_o]],header=head,fmt='%d \t %.7f \t %.7f \t %.7f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f')

print filein + ' Done!'
