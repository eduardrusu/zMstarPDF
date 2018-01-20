# Randomly varies the input parameters for Hostlens, then delects all outputs except for the one with minimum chi^2

import sys
import os
import numpy as np

file = "ylens_out_file.input"
prefix_out = "ylens_out"
vary = 100

def string(hostlens):
    return "  " + hostlens.split()[2] + "   " + hostlens.split()[3] + " " + hostlens.split()[4] + " " + hostlens.split()[5] + "\n"
def stringunit(hostlens): # for lines that end with e.g. "[arcsec]"
    return "  " + hostlens.split()[2] + "   " + hostlens.split()[3] + " " + hostlens.split()[4] + " " + hostlens.split()[5] + " " + hostlens.split()[6] + "\n"

for i in range(vary):
    with open(file, 'r') as f:
        hostlens = f.readlines()
    hostlens[5 - 1] = "prefix        %s%s" % (prefix_out,str(i + 1)) + "\n"
# psf + sky
    line = 30; x = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(0.8,0.2))) + stringunit(hostlens[line - 1]); hostlens[line - 1] = x # FWHM1
    line = 31; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(0.1,0.05))) + string(hostlens[line - 1]) # e1
    line = 32; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (np.random.uniform(0,360)) + stringunit(hostlens[line - 1]) # PA1
    line = 33; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(7,2))) + string(hostlens[line - 1]) # beta1
    line = 34; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(1.2,0.3))) + stringunit(hostlens[line - 1]) # FWHM2
    line = 35; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(0.1,0.05))) + string(hostlens[line - 1]) # e2
    line = 36; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (np.random.uniform(0,360)) + stringunit(hostlens[line - 1]) # PA2
    line = 37; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(2,1))) + string(hostlens[line - 1]) # beta2
    line = 38; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (np.random.uniform(0.001,0.999)) + "  " + hostlens[line - 1].split()[2] + "   " + hostlens[line - 1].split()[3] + " # flux1 / (flux1 + flux 2) \n"
    line = 39; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (np.random.normal(0,10)) + stringunit(hostlens[line - 1]) # sky
# point source
    line = 43; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(78.8,1))) + stringunit(hostlens[line - 1]) # x
    line = 44; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(41.0,1))) + stringunit(hostlens[line - 1]) # y
    line = 50; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(5.0e+06,1.0e+06))) + stringunit(hostlens[line - 1]) # flux_point
    line = 54; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(83.4,1))) + stringunit(hostlens[line - 1]) # x
    line = 55; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(43.6,1))) + stringunit(hostlens[line - 1]) # y
    line = 61; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(6.0e+06,1.0e+06))) + stringunit(hostlens[line - 1]) # flux_point
    line = 65; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(83.4,1))) + stringunit(hostlens[line - 1]) # x
    line = 66; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(43.6,1))) + stringunit(hostlens[line - 1]) # y
    line = 72; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(3.0e+06,1.0e+06))) + stringunit(hostlens[line - 1]) # flux_point
    line = 76; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(85.0,1))) + stringunit(hostlens[line - 1]) # x
    line = 77; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(30.7,1))) + stringunit(hostlens[line - 1]) # y
    line = 83; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(6.0e+05,3.0e+05))) + stringunit(hostlens[line - 1]) # flux_point
# galaxy
    line = 87; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(84.5,1))) + stringunit(hostlens[line - 1]) # x
    line = 88; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(34.5,1))) + stringunit(hostlens[line - 1]) # y
    line = 89; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(1.0e+05,1.0e+05))) + stringunit(hostlens[line - 1]) # flux_galaxy
    #line = 92; hostlens[line - 1] = hostlens[line - 1].split()[0] + "   " + "%.6e" % (abs(np.random.normal(0.5,0.2))) + stringunit(hostlens[line - 1]) # re

    with open(file, 'w') as f:
        f.writelines(hostlens)
    f.close()

    os.system("hostlens %s" % file)

chi = np.zeros(vary)

for i in range(vary):
    with open(prefix_out + str(i + 1) + "_optresult.dat", 'r') as f:
        hostlens = f.readlines()
    j = np.shape(hostlens)[0]
    while "chi^2/nu" not in hostlens[j - 1]: j -= 1
    chi[i] = float(hostlens[j - 1].split()[4])

save = np.where(chi == np.min(chi))[0][0] + 1
for i in range(vary):
    if i + 1 != save:
        os.system("rm %s%s_file.input" % (prefix_out,str(i + 1)))
        os.system("rm %s%s_optresult.dat" % (prefix_out,str(i + 1)))
        os.system("rm %s%s_subtract.fits" % (prefix_out,str(i + 1)))




