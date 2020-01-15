# Given a list of MCMC chains, do corner plots and assess convergence
# run as python mcmc_multipleprogression.py frac file1 [file2] [...]
# Above, "frac" is the fraction of the chain length at which burn-in ends

import sys
import numpy as np
import pylab as plt
import corner

frac = float(sys.argv[1])
chains = len(sys.argv) - 2
files = sys.argv
print "Using %s chains..." % chains

list = []
for i in range(chains):#10
    mcmc = np.loadtxt(sys.argv[2+i],unpack=True)
    mcmci = mcmc[1:,int(np.shape(mcmc)[1]*frac):np.shape(mcmc)[1]]
    # eliminate the first column, containing chi^2, and the burn-in
    if i == 0:
        mcmcfinal = mcmci
        minsize = np.shape(mcmci)[1]
    else:
        mcmcfinal = np.append(mcmcfinal,mcmci, axis = 1)
        minsize = np.min([np.shape(mcmci)[1],minsize])
    list += [mcmci]


# Convergence diagnostics with Gelman-Rubin 1995 R_hat
def R_hat(samples): # https://groups.google.com/forum/#!topic/hddm-users/qWzCWTz-wFQ
# formulae in https://pymc-devs.github.io/pymc/modelchecking.html
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

for i in range(np.shape(mcmci)[0]):
    for j in range(chains):
        if j == 0: samples = list[0][i][:minsize]
        else: samples = np.vstack((samples,list[j][i][:minsize]))
    print "[%d] R_hat = " %(i+1), R_hat(samples)

figure = corner.corner(mcmcfinal[0:np.shape(mcmcfinal)[0]].T, labels=np.linspace(1,np.shape(mcmcfinal)[0],\
np.shape(mcmcfinal)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})

# program to find the stem of given list of words function to find the stem (longest common substring) from the string array
# This code is contributed by ita_c (https://www.geeksforgeeks.org/longest-common-substring-array-strings/)
def findstem(arr):
    # Determine size of the array
    n = len(arr)
    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)
    res = ""
    for i in range( l) :
        for j in range( i + 1, l + 1) :
            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(1, n):
                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break
            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if (k + 1 == n and len(res) < len(stem)):
                res = stem
    return res

stems = findstem(files[2:])
figure.savefig("%smcmc.png" % stems, dpi=100)
