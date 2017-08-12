import emcee
import corner
import numpy as np
import pylab as plt
from matplotlib.ticker import MaxNLocator

#np.random.seed(42) # for repeatability
theta_true = (25, 0.5)
xdata = 100 * np.random.random(20)
ydata = theta_true[0] + theta_true[1] * xdata
ydata = np.random.normal(ydata, 10) # add error

def log_prior(theta):
    alpha, beta, sigma = theta
    if sigma < 0:
        return -np.inf # log(0)
    #else:
        #return (-1.5 * np.log(1 + beta**2) - np.log(sigma))
    else:
        return 0.0 # flat prior
def log_like(theta, x, y):
    alpha, beta, sigma = theta
    y_model = alpha + beta * x
    return -0.5 * np.sum(np.log(2*np.pi*sigma**2) + (y-y_model)**2 / sigma**2)
def log_posterior(theta, x, y):
    return log_prior(theta) + log_like(theta,x,y)

ndim = 3 # number of parameters in the model
nwalkers = 50 # number of MCMC walkers
nburn = 1000 # "burn-in" to stabilize chains
nsteps = 10000 # number of MCMC steps to take after burn-in
starting_guesses = np.random.rand(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim,log_posterior,args=[xdata,ydata])
#sampler = emcee.EnsembleSampler(nwalkers, ndim,log_posterior,args=[xdata,ydata],threads = 4) # multiple processors
pos, prob, state = sampler.run_mcmc(starting_guesses, nburn)
# plot the time laps, only for the burn-in

plt.clf()
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(theta_true[0], color="#888888", lw=2)
axes[0].set_ylabel("$a$")
axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(theta_true[1], color="#888888", lw=2)
axes[1].set_ylabel("$b$")
axes[2].plot(np.exp(sampler.chain[:, :, 2]).T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].set_ylabel("$s$")
axes[2].set_xlabel("step number")
fig.tight_layout(h_pad=0.0)
fig.show()
#fig.savefig("line-time.png")

sampler.reset() # save in pos the position of the walkers at the end of the burn-in, and rerun from there
#sampler.run_mcmc(pos, nsteps) # instead of this, I am doing as follows so I can display the progress:
for i, result in enumerate(sampler.sample(pos,iterations=nsteps)): # fails unless I keep the keyword *iterations*
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nsteps))

alpha_samp = sampler.flatchain.T[0] # combines all walkers
beta_samp = sampler.flatchain.T[1]
sigma_samp = sampler.flatchain.T[2]
print("Autocorrelation time:", sampler.get_autocorr_time())
print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
print "median, std 1: ", np.median(alpha_samp), np.std(alpha_samp)
print "median, std 2: ", np.median(beta_samp), np.std(beta_samp)
print "median, std 3: ", np.median(sigma_samp), np.std(sigma_samp)
fig = corner.corner(sampler.flatchain, labels=["$a$", "$b$", "$s$"],truths=[25, 0.5, 10])
fig.show()

# emcee can probe multiple minima by using different temperatures; Emcee includes a PTSampler (PT = parallel tempering) that has a method PTSampler.thermodynamic_integration_log_evidence() for performing the evidence integral after the sampler is run
# saving progress
# there is an ipython parallel example in the emcee code. how does that work? also, look into the loadbalance example; there is also subprocessing.py example to run on multiple computers
