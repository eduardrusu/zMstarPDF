# Given an image containing flux from a perfect model, as well as a value for the original (unsubtracted) sky, create an image containing pure Poisson noise, with null sky level
# Careful, the code needs to be updates for each particular input file

import sys
import os
import numpy as np
from astropy.io import fits
import corner

exptime = 1440 # exposure time
sim = 100 # number of simulations
model = "ilens_out_model.fits"
sky = 40000 # in case I compute this from the std of sky pixels, I need to do it on the subtracted best-fit image
file = "ilens_out_file_simulated.input"

for i in range(sim):# sim
    
    # simulate noisy image
    image = fits.open(model)
    data = image[0].data
    data = data + sky
    noise = np.random.poisson(data)
    noise = noise - sky
    image[0].data = noise.astype("float32") # Hostlens expects this
    image.writeto(model[:-5]+"_noise.fits",clobber=True)
    
    # redirect the output
    with open(file, 'r') as f:
        hostlens = f.readlines()
    hostlens[5 - 1] = "prefix        %s" % file[:-6] + str(i+1) + "\n"
    if i == 0:
        hostlens[24 - 1] = "dump_subtract 1 \n"
    else:
        hostlens[24 - 1] = "dump_subtract 0 \n"
    if i == sim - 1:
        pixscale = float(hostlens[6 - 1].split()[1])
        # update these:
        init1 = float(hostlens[43 - 1].split()[1]) * pixscale
        init2 = float(hostlens[44 - 1].split()[1]) * pixscale
        init3 = 25 - 2.5*np.log10(float(hostlens[50 - 1].split()[1])/exptime)
        init4 = float(hostlens[54 - 1].split()[1]) * pixscale
        init5 = float(hostlens[55 - 1].split()[1]) * pixscale
        init6 = 25 - 2.5*np.log10(float(hostlens[61 - 1].split()[1])/exptime)
        init7 = float(hostlens[65 - 1].split()[1]) * pixscale
        init8 = float(hostlens[66 - 1].split()[1]) * pixscale
        init9 = 25 - 2.5*np.log10(float(hostlens[72 - 1].split()[1])/exptime)
        init10 = float(hostlens[76 - 1].split()[1]) * pixscale
        init11 = float(hostlens[77 - 1].split()[1]) * pixscale
        init12 = 25 - 2.5*np.log10(float(hostlens[83 - 1].split()[1])/exptime)
        init13 = float(hostlens[87 - 1].split()[1]) * pixscale
        init14 = float(hostlens[88 - 1].split()[1]) * pixscale
        init15 = 25 - 2.5*np.log10(float(hostlens[89 - 1].split()[1])/exptime)
    with open(file, 'w') as f:
        f.writelines(hostlens)
        f.close()

    os.system("hostlens %s" % file)

# update number of parameters according to input
par1 = np.zeros(sim)
par2 = np.zeros(sim)
par3 = np.zeros(sim)
par4 = np.zeros(sim)
par5 = np.zeros(sim)
par6 = np.zeros(sim)
par7 = np.zeros(sim)
par8 = np.zeros(sim)
par9 = np.zeros(sim)
par10 = np.zeros(sim)
par11 = np.zeros(sim)
par12 = np.zeros(sim)
par13 = np.zeros(sim)
par14 = np.zeros(sim)
par15 = np.zeros(sim)

for i in range(sim):
    with open(file[:-6] + str(i+1) + "_optresult.dat", 'r') as f:
        hostlens = f.readlines()
    # updte these:
    par1[i] = float(hostlens[0 - 7].split()[1]) * pixscale # -7 if there are 5 objects of interest
    par2[i] = float(hostlens[0 - 7].split()[2]) * pixscale
    par3[i] = 25 - 2.5*np.log10(float(hostlens[0 - 7].split()[8])/exptime)
    par4[i] = float(hostlens[0 - 6].split()[1]) * pixscale
    par5[i] = float(hostlens[0 - 6].split()[2]) * pixscale
    par6[i] = 25 - 2.5*np.log10(float(hostlens[0 - 6].split()[8])/exptime)
    par7[i] = float(hostlens[0 - 5].split()[1]) * pixscale
    par8[i] = float(hostlens[0 - 5].split()[2]) * pixscale
    par9[i] = 25 - 2.5*np.log10(float(hostlens[0 - 5].split()[8])/exptime)
    par10[i] = float(hostlens[0 - 4].split()[1]) * pixscale
    par11[i] = float(hostlens[0 - 4].split()[2]) * pixscale
    par12[i] = 25 - 2.5*np.log10(float(hostlens[0 - 4].split()[8])/exptime)
    par13[i] = float(hostlens[0 - 3].split()[1]) * pixscale
    par14[i] = float(hostlens[0 - 3].split()[2]) * pixscale
    par15[i] = 25 - 2.5*np.log10(float(hostlens[0 - 3].split()[3])/exptime)
    #os.system("rm %s" % (file[:-6] + str(i+1) + "_optresult.dat"))

# updte these:
std1 = np.sqrt(np.sum((par1 - init1)**2)/(sim - 1))
std2 = np.sqrt(np.sum((par2 - init2)**2)/(sim - 1))
std3 = np.sqrt(np.sum((par3 - init3)**2)/(sim - 1))
std4 = np.sqrt(np.sum((par4 - init4)**2)/(sim - 1))
std5 = np.sqrt(np.sum((par5 - init5)**2)/(sim - 1))
std6 = np.sqrt(np.sum((par6 - init6)**2)/(sim - 1))
std7 = np.sqrt(np.sum((par7 - init7)**2)/(sim - 1))
std8 = np.sqrt(np.sum((par8 - init8)**2)/(sim - 1))
std9 = np.sqrt(np.sum((par9 - init9)**2)/(sim - 1))
std10 = np.sqrt(np.sum((par10 - init10)**2)/(sim - 1))
std11 = np.sqrt(np.sum((par11 - init11)**2)/(sim - 1))
std12 = np.sqrt(np.sum((par12 - init12)**2)/(sim - 1))
std13 = np.sqrt(np.sum((par13 - init13)**2)/(sim - 1))
std14 = np.sqrt(np.sum((par14 - init14)**2)/(sim - 1))
std15 = np.sqrt(np.sum((par15 - init15)**2)/(sim - 1))

# updte these:
print "std 1: ", std1, " [arcsec]"
print "std 2: ", std2, " [arcsec]"
print "std 3: ", std3, " [dmag]"
print "std 4: ", std4, " [arcsec]"
print "std 5: ", std5, " [arcsec]"
print "std 6: ", std6, " [dmag]"
print "std 7: ", std7, " [arcsec]"
print "std 8: ", std8, " [arcsec]"
print "std 9: ", std9, " [dmag]"
print "std 10: ", std10, " [arcsec]"
print "std 11: ", std11, " [arcsec]"
print "std 12: ", std12, " [dmag]"
print "std 13: ", std13, " [arcsec]"
print "std 14: ", std14, " [arcsec]"
print "std 15: ", std15, " [dmag]"

# update the length:
#data = np.c_[par1,par2,par3,par4,par5,par6,par7,par8,par9,par10,par11,par12,par13,par14,par15]
#figure = corner.corner(data, labels=np.linspace(1,np.shape(data)[0],np.shape(data)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
#figure.savefig(file[:-6] + "_simulations.png", dpi=100)
# ignore any WARNING:root:Too few points to create valid contours


