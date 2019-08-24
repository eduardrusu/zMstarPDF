# Creates both time laps and corner plots for MCMC output files in the format of glafic and hostlens
# Use as python mcmcplots.py filename

import numpy as np
import sys

file = str(sys.argv[1])
filetime = file[:-4] + "_timelapse.png"
filecorner = file[:-4] + "_cornerplot.png"

burnin = 10 # in percent, to eliminate from the head
data = np.loadtxt(file,unpack=True)
data = data[0:len(data),int(burnin/100.0*len(data[0])):len(data[0])]

import corner
import pylab as plt
from matplotlib.ticker import MaxNLocator

# plot timelapse
plt.clf()
fig, axes = plt.subplots(len(data), 1, sharex=True, figsize=(20, 20))
for i in range(len(data)):
    axes[i].plot(data[i], color="k", alpha=0.4)
    axes[i].yaxis.set_major_locator(MaxNLocator(5))
    #axes[0].axhline(chimin_z, color="#888888", lw=2)
    axes[i].set_ylabel(str(i))
axes[i].set_xlabel("step number")
fig.tight_layout(h_pad=0.0)
#fig.show()
fig.savefig(filetime, dpi=150)

# corner plot
#nburn = nsteps/2 # "burn-in" to stabilize chains
#samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
#alpha_samp = sampler.flatchain.T[0]
#beta_samp = sampler.flatchain.T[1]
#sigma_samp = sampler.flatchain.T[2]
#print("Autocorrelation time:", sampler.get_autocorr_time())
#print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
labels = np.linspace(0,len(data)-1,len(data)).astype(int).astype(str)
meds = np.zeros(len(data))
for i in range(len(data)):
    meds[i] = np.median(data[i])
fig = corner.corner(data.T, labels=labels,truths=meds)
#fig.show()
fig.savefig(filecorner, dpi=150)
