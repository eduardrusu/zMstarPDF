# This code samples from a single-band the Millenium Simulation (MS) photometry

import numpy as np

rootin = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/original/"
rootout = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/PG1115/"

fileout_ugriz = rootout+str(sys.argv[1])[:-11] + "_ugriz.images.txt"
fileout_ugrizJHK = rootout+str(sys.argv[1])[:-11] + "_%s_%s.images.txt" % (filters,lens)

id = 0
zspec = 5
posx = 6
posy = 7
r = 14

rRstd = 0.07
head = "GalID \t z_spec \t pos0 \t pos_1 mag_SDSS_r"
for i in range(1):
    for j in range(1):
        for k in range(1):
            for l in range(1):
                filein = rootin+'GGL_los_8_%s_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images' %(i,j,k,l)
                fileout = rootout+'GGL_los_8_%s_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_originalr.images' %(i,j,k,l)
                data = np.loadtxt(filein,usecols=[id,zspec,posx,posy,r],comments='GalID')
                data = data[data[:,4]<=24]
                np.savetxt(fileout,data.T,header=head,fmt='%d \t %.3f \t %.7f \t %.7f \t %.2f')
                fileout = rootout+'GGL_los_8_%s_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_sampledr.images' %(i,j,k,l)
                data = np.loadtxt(filein,usecols=[id,zspec,posx,posy,r],comments='GalID')
                data[:,4] = np.random.normal(dataout[:,4], 0.07)
                data = data[data[:,4]<=24]
                np.savetxt(fileout,data.T,header=head,fmt='%d \t %.3f \t %.7f \t %.7f \t %.2f')
