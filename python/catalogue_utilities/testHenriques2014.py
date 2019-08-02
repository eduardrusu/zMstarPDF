# extract Henriques 2014 galaxies from *.images.ugriz.txt
import numpy as np
pl = np.linspace(30,62,62-30 + 1)
out = np.empty(6)
for i in range(len(pl)):
    print int(pl[i])
    x = np.loadtxt("GGL_los_8_0_0_N_4096_ang_4_Henriques2014_galaxies_on_plane_%s_f.images.ugriz.txt" % int(pl[i]),usecols=[4,5,6,7,8,9])
    j = 0
    ind = 0
    while (j < 1000) and (ind < len(x)):
        if x[ind][4]<23:
            out = np.c_[out,x[ind]]
            j += 1
        ind += 1
nr = np.linspace(1,len(out[0]),len(out[0]))
np.savetxt("Henriques2014selectforBPZ.cat",np.c_[nr,out[0],out[1],out[2],out[3],out[4],out[5]],fmt='%d %1.4f %1.2f 0.00 %1.2f 0.00 %1.2f 0.00 %1.2f 0.00 %1.2f 0.00', header = 'id z_s u u_err g g_err r r_err i i_err z z_err')

# now also extract from SA
import numpy as np
out = np.loadtxt("/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/original/GGL_los_8_0_0_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt",usecols=[5,12,13,14,15,16], comments="GalID",unpack=True)
nr = np.linspace(1,len(out[0]),len(out[0]))
np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/0408/SAselectforBPZ.cat",np.c_[nr,out[0],out[1],out[2],out[3],out[4],out[5]],fmt='%d %1.4f %1.2f 0.00 %1.2f 0.00 %1.2f 0.00 %1.2f 0.00 %1.2f 0.00', header = 'id z_s u u_err g g_err r r_err i i_err z z_err')

# now extract from .images
