# Given an image containing flux from a best-fit model and its generating hostlens configuration file, simulate similar images with noise, then fit them by starting from the correct input parameter values as initial guesses, and compute the scatter from these values

import numpy as np
from astropy.io import fits
import os

model = "fixsersstars_sdfbFCSA00205330cutcutfliprotate_model.fits"
sky = 2250 # value efore sky subtraction
skyrms = 25 # measured on empty sky regions of the original image with imexam

image = fits.open(model)
data = image[0].data + sky
sims = 10
out = np.zeros((sims,np.shape(data)[0],np.shape(data)[1]))

for i in range(sims):
    out[i] = np.random.poisson(abs(data)) - sky + np.random.normal(0,25,size=(np.shape(data)[0],np.shape(data)[1]))
    image[0].data = out[i].astype("float32")
    image.writeto(model[:-5]+"_noise%s.fits" % str(i),clobber=True)

for i in range(sims):
    os.system("cp %s %s" % (model[:-11] + "_file.input",model[:-11] + str(i) + "_file.input"))
    with open(model[:-11] + str(i) + "_file.input", 'r') as f:
        hostlens = f.readlines()
    hostlens[1 - 1] = "obsfits       %s" % (model[:-5]+"_noise%s.fits" %str(i)) + "\n"
    hostlens[5 - 1] = "prefix        %s" % (model[:-11] + str(i)) + "\n"
    hostlens[22 - 1] = "dump_model    0" + "\n"
    with open(model[:-11] + str(i) + "_file.input", 'w') as f:
        f.writelines(hostlens)
        f.close()
    if i == 0: # identify the variables
        var = []
        with open(model[:-11] + str(i) + "_file.input", 'r') as g:
            hostlens = g.readline()
            while hostlens:
                hostlens = g.readline()
                if (len(hostlens.split()) >= 3) and (hostlens.split()[0].isdigit() == True) and (hostlens.split()[2] == '1'):
                    var.append(float(hostlens.split()[1]))

execute = ""
for i in range(sims):
    if i < sims - 1: execute += "hostlens %s & " % (model[:-11] + str(i) + "_file.input")
    else: execute += "hostlens %s " % (model[:-11] + str(i) + "_file.input")
os.system(execute)

fitted = np.zeros((sims,len(var)))
for i in range(sims):
    fit = []
    with open(model[:-11] + str(i) + "_file.input", 'r') as g:
        hostlens = g.readline()
        while hostlens:
            hostlens = g.readline()
            if (len(hostlens.split()) >= 3) and (hostlens.split()[0].isdigit() == True) and (hostlens.split()[2] == '1'):
                fit.append(float(hostlens.split()[1]))
    fitted[i] = fit

print "Standard deviations from given values:"
print np.sqrt(np.std(fitted,axis=0,ddof=1)**2 + (var-np.mean(fitted,axis=0))**2)
