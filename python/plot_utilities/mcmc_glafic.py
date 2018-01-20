# Given the best-fit input file from glafic, run 10 mcmc chains, combine them (accounting for burn-in) and plot the histogram with the 16th and 84th percentiles, using the corner package. It also computes the R_hat convergence diagnostic.

import sys
import os
import numpy as np
import corner

file = str(sys.argv[1])
length = 10000000
burnin = 0.5 # fraction of the chain to remove

for i in range(10): #10
    os.system("cp %s %s" % (file,file[:-6] + str(i+1) + ".input"))
    with open(file[:-6] + str(i+1) + ".input", 'r') as f:
        glafic = f.readlines()
    glafic[10 - 1] = "prefix        %s" % file[:-6] + str(i+1) + "\n"
    glafic[20 - 1] = "ran_seed      -%s" % str(i+1) + "\n"
    glafic[55 - 1] = "mcmc       %s" % str(length) + "\n"
    with open(file[:-6] + str(i+1) + ".input", 'w') as f:
        f.writelines(glafic)
        f.close()

os.system("glafic %s & glafic %s & glafic %s & glafic %s & glafic %s & glafic %s & glafic %s & glafic %s & glafic %s & glafic %s" % (file[:-6] + "1.input",file[:-6] + "2.input",file[:-6] + "3.input",file[:-6] + "4.input",file[:-6] + "5.input",file[:-6] + "6.input",file[:-6] + "7.input",file[:-6] + "8.input",file[:-6] + "9.input",file[:-6] + "10.input"))

minlen = 0
for i in range(10):#10 # read to find the lengths, since for stacking they have to have equal lengths
    mcmc = np.loadtxt(file[:-6] + str(i+1) + "_mcmc.dat",unpack=True)
    if (minlen == 0) or (len(mcmc[0]) < minlen): minlen = len(mcmc[0])

for i in range(10):#10
    mcmc = np.loadtxt(file[:-6] + str(i+1) + "_mcmc.dat",unpack=True)
    mcmci = mcmc[1:,np.shape(mcmc)[1] - minlen + int(minlen * burnin):] # eliminate the first column, containing chi^2, as well as the burnin
    if i == 0: mcmcfinal = mcmci
    else: mcmcfinal = np.append(mcmcfinal,mcmci, axis = 1)
    
    if i == 0: mcmc0 = mcmci
    if i == 1: mcmc1 = mcmci
    if i == 2: mcmc2 = mcmci
    if i == 3: mcmc3 = mcmci
    if i == 4: mcmc4 = mcmci
    if i == 5: mcmc5 = mcmci
    if i == 6: mcmc6 = mcmci
    if i == 7: mcmc7 = mcmci
    if i == 8: mcmc8 = mcmci
    if i == 9: mcmc9 = mcmci

# Convergence diagnostics with Gelman-Rubin 1995 R_hat
def R_hat(samples): # https://groups.google.com/forum/#!topic/hddm-users/qWzCWTz-wFQ # formulas in https://pymc-devs.github.io/pymc/modelchecking.html
    m, n = np.shape(samples) # m = chains, n = samples
    # Chain variance
    chain_var = np.var(samples, axis=1, ddof=1) # degrees of freedom = n-ddof
    # Within-chain variance (mean of variances of each chain)
    W = 1./m * np.sum(chain_var)
    # Chain means
    chain_means = np.mean(samples, axis=1)
    # mean_of_means = numpy.mean(chain_means) # all chains have same length
    # Variance of chain means
    chain_means_var = np.var(chain_means, ddof=1)
    # Between-chain variance
    B = n * chain_means_var
    # Weighted average of within and between variance
    #(marginal posterior variance)
    Var_hat = (float(n-1)/n)*W + B/n
    # Potential scale reduction factor
    R_hat = np.sqrt(Var_hat / W)
    return R_hat

for i in range(np.shape(mcmc0)[0]):
    samples = np.vstack((mcmc0[i],mcmc1[i],mcmc2[i],mcmc3[i],mcmc4[i],mcmc5[i],mcmc6[i],mcmc7[i],mcmc8[i],mcmc9[i]))
    print "[%d] R_hat = " %i, R_hat(samples)

figure = corner.corner(mcmcfinal[0:np.shape(mcmcfinal)[0]].T, labels=np.linspace(1,np.shape(mcmcfinal)[0],np.shape(mcmcfinal)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig(file[:-6] + "_mcmc.png", dpi=100)

