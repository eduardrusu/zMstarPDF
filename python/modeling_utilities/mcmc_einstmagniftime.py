# Uses an existent glafic MCMC chain to plot corresponding distributions for the einstein radius, time delays and magnifications
# run in a fresh terminal, because I need to save the terminal output to disk

import sys
import os
import numpy as np
import corner

file = "pointSIEgamma.input"
fileout = "out_SIEgammafield_einstmagniftime_point.dat"
chains = 10
length = 10000
z_s = 2.517
img = 4

for i in range(chains):
    mcmc = np.loadtxt(file[:-6] + str(i+1) + "_mcmc.dat",unpack=True)
    mcmci = mcmc[1:,int(np.shape(mcmc)[1]/4):np.shape(mcmc)[1]] # eliminate the first column, containing chi^2, as well as the first 25% of the chains
    if i == 0: mcmcfinal = mcmci
    else: mcmcfinal = np.append(mcmcfinal,mcmci, axis = 1)
mcmclength = np.shape(mcmcfinal)[1]
interval = mcmclength / length

for i in range(length):
    if interval * i < mcmclength:
        print 'step', i+1
        os.system("cp %s %s" % (file,file[:-6] + "_einstmagniftime.input"))
        with open(file[:-6] + "_einstmagniftime.input", 'r') as f:
            glafic = f.readlines()
        glafic[10 - 1] = "prefix      %s \n" % (fileout[:-10])
        glafic[28 - 1] = "lens   sie      %s  %s  %s  %s  %s  0.000000e+00  0.000000e+00 \n" % (mcmcfinal[:,interval * i][0],mcmcfinal[:,interval * i][1],mcmcfinal[:,interval * i][2],mcmcfinal[:,interval * i][3],mcmcfinal[:,interval * i][4])
        glafic[29 - 1] = "lens   pert     %s  %s  %s  %s  %s  0.000000e+00  0.000000e+00 \n" % (z_s,mcmcfinal[:,interval * i][5],mcmcfinal[:,interval * i][6],mcmcfinal[:,interval * i][7],mcmcfinal[:,interval * i][8])
        glafic[21 - 1] = "chi2_checknimg 1 \n"
        glafic[36 - 1] = "0 0 0 0 0 0 0 \n"
        glafic[37 - 1] = "0 0 0 0 0 0 0 \n"
        glafic[51 - 1] = "\n"
        with open(file[:-6] + "_einstmagniftime.input", 'w') as f:
            f.writelines(glafic)
            f.close()
        os.system("glafic %s" % (file[:-6] + "_einstmagniftime.input"))
        x = np.loadtxt(fileout)
        if x[0][0] == img:
            with open(file[:-6] + "_einstmagniftime_out_.dat", 'a') as f:
                f.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n" % (x[1][0],x[1][1],x[1][2],x[1][3],x[2][0],x[2][1],x[2][2],x[2][3],x[3][0],x[3][1],x[3][2],x[3][3],x[4][0],x[4][1],x[4][2],x[4][3]))
                f.close()
