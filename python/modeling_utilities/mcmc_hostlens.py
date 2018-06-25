# Given the best-fit input file from hostlens, run 10 mcmc chains, combine them (accounting for burn-in) and plot the histogram with the 16th and 84th percentiles, using the corner package. It also computes the R_hat convergence diagnostic.

import sys
import os
import numpy as np
import corner

file = "ylens_out23_file.input"
length = 1000000
for i in range(0): #10
    os.system("hostlens %s" % file) # since it takes a few runs to fully converge

for i in range(10): #10
    os.system("cp %s %s" % (file,file[:-10] + str(i+1) + "_file.input"))
    with open(file[:-10] + str(i+1) + "_file.input", 'r') as f:
        hostlens = f.readlines()
    hostlens[5 - 1] = "prefix        %s" % file[:-10] + str(i+1) + "\n"
    hostlens[17 - 1] = "ran_seed      %s" % str(i) + "\n"
    hostlens[20 - 1] = "do_optimize   0 \n"
    hostlens[21 - 1] = "dump_psf      0 \n"
    hostlens[22 - 1] = "dump_model    0 \n"
    hostlens[23 - 1] = "dump_modelori 0 \n"
    hostlens[24 - 1] = "dump_subtract 0 \n"
    hostlens[25 - 1] = "dump_infile   0 \n"
    hostlens[26 - 1] = "dump_lenseq   0 \n"
    hostlens[27 - 1] = "do_mcmc       %s" % str(length) + "\n"
    with open(file[:-10] + str(i+1) + "_file.input", 'w') as f:
        f.writelines(hostlens)
        f.close()

#os.system("hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s & hostlens %s" % (file[:-10] + "1_file.input",file[:-10] + "2_file.input",file[:-10] + "3_file.input",file[:-10] + "4_file.input",file[:-10] + "5_file.input",file[:-10] + "6_file.input",file[:-10] + "7_file.input",file[:-10] + "8_file.input",file[:-10] + "9_file.input",file[:-10] + "10_file.input"))

for i in range(10):#10
    mcmc = np.loadtxt(file[:-10] + str(i+1) + "_mcmc.dat",unpack=True)
    mcmci = mcmc[1:,int(np.shape(mcmc)[1]/5):np.shape(mcmc)[1]] # eliminate the first column, containing chi^2, as well as the first 20% of the chains
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

