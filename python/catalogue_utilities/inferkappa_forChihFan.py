# CE Rusu July 9 2018
# kappagammaforChihFan.py should be run first to get the files kappagamma_values_GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f.dat etc.
# I am not correcting for bias by dividing by the number of lines of sight inside subintervals because the shears from Ken are very narrow

import numpy as np
import time

start_time=time.time()

#models = ['noshear','powNFWG1G2','powNFWG1','powNFW','compNFW','powSISG1']
#shear = [0,0.047,0.038,0.041,0.050,0.011]
#shearerr = [10,0.0062,0.0095,0.0099,0.0075,0.0070]
models = ['powSISG1']
shear = [0.011]
shearerr = [0.0070]

root = "/lfs08/rusucs/kappaplanes/"

print "Reading..."

for j in range(8):
    for i in range(8):
        print j,i
        kappa_,gamma_ = np.loadtxt("%skappagamma_values_GGL_los_8_%s_%s_N_4096_ang_4_rays_to_plane_34_f.dat" % (root,str(j),str(i)), unpack=True)
        if i == 0 and j == 0:
            kappa = kappa_
            gamma = gamma_
        else:
            kappa = np.append(kappa,kappa_)
            gamma = np.append(gamma,gamma_)

bin_stat = 3000
min_kappa = -0.50
max_kappa = 1

for i in range(len(models)):
    output = '%skappagamma_values_GGL_los_8_N_4096_ang_4_rays_to_plane_34_f_%s.dat' % (root,models[i])
    kappafinal = kappa[(gamma > shear[i] - shearerr[i]) & (gamma < shear[i] + shearerr[i])]
    kappahist = np.histogram(kappafinal, bins = bin_stat, range=(min_kappa,max_kappa))[0].astype(float)
    kappaindex = np.histogram(kappafinal, bins = bin_stat, range=(min_kappa,max_kappa))[1].astype(float)
    np.savetxt(output,np.c_[kappaindex[:-1],kappahist],fmt='%s',delimiter='\t',newline='\n',header="%s" % len(kappafinal))

print(" Total time --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
