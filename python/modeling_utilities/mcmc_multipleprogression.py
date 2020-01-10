# Given a list of MCMC chains, plot them all to assess burn-in and convergence
# run as python mcmc_multipleprogression.py file1 [file2] [...]

import sys
import numpy as np
import pylab as plt

chains = len(sys.argv) - 1
print "Using %s chains..." % chains
files = sys.argv

plt.clf()
for j in range(chains):
    x=np.loadtxt(sys.argv[j + 1],unpack=True)
    ndim = np.shape(x)[0] - 1
    if j == 0: fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(5, 1*10*ndim))
    x = x[1:]
    for i in range(ndim):
        axes[i].plot(x[i],linewidth=0.1,alpha=0.8)
    if i ==0: axes[i].set_ylabel("%s" %i)
axes[i].set_xlabel("step number")
fig.tight_layout(h_pad=0.0)

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

stems = findstem(files[1:])
fig.savefig("%stimeline.png" % stems)
