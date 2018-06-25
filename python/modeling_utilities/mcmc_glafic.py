# Given the best-fit input file from glafic, run 10 mcmc chains, combine them (accounting for burn-in) and plot the histogram with the 16th and 84th percentiles, using the corner package. It also computes the R_hat convergence diagnostic.

import sys
import os
import numpy as np
import corner

file = "point1pertSIEgamma.input"
length = 400000

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

for i in range(10):#10
    mcmc = np.loadtxt(file[:-6] + str(i+1) + "_mcmc.dat",unpack=True)
    mcmci = mcmc[1:,int(np.shape(mcmc)[1]/4):np.shape(mcmc)[1]] # eliminate the first column, containing chi^2, as well as the first 25% of the chains
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

minsize = np.min([np.shape(mcmc0)[1],np.shape(mcmc1)[1],np.shape(mcmc2)[1],np.shape(mcmc3)[1],np.shape(mcmc4)[1],np.shape(mcmc5)[1],np.shape(mcmc6)[1],np.shape(mcmc7)[1],np.shape(mcmc8)[1],np.shape(mcmc9)[1]])
for i in range(np.shape(mcmc)[0] - 1): # added -1
    samples = np.vstack((mcmc0[i][:minsize],mcmc1[i][:minsize],mcmc2[i][:minsize],mcmc3[i][:minsize],mcmc4[i][:minsize],mcmc5[i][:minsize],mcmc6[i][:minsize],mcmc7[i][:minsize],mcmc8[i][:minsize],mcmc9[i][:minsize]))
    print "[%d] R_hat = " %(i+1), R_hat(samples) # added +1

figure = corner.corner(mcmcfinal[0:np.shape(mcmcfinal)[0]].T, labels=np.linspace(1,np.shape(mcmcfinal)[0],np.shape(mcmcfinal)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig(file[:-6] + "_mcmc.png", dpi=100)

