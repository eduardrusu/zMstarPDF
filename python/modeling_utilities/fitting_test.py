# Tools for simulating a Seric profile, convolving it with the PSF, optimizing it and computing parameter uncertinties with emcee

import numpy as np
import time

from astropy.io import fits
pixels = 100
# create a blanck image
img = np.zeros([pixels,pixels])
hdu = fits.PrimaryHDU(img)
hdu.writeto('new.fits',clobber=True)

# create a fits file with a Sersic
from astropy.modeling.models import Sersic2D
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
sers = Sersic2D(amplitude = 100, r_eff = 10, n=1, x_0=50, y_0=50, ellip=.5, theta=-1)
# amplitude = Central surface brightness, within r_eff
x,y = np.meshgrid(np.arange(pixels), np.arange(pixels))
img = sers(x,y)
hdu = fits.PrimaryHDU(img)
hdu.writeto('new.fits',clobber=True)

# insert noise to the model
originalsky = 20
noise = np.random.poisson(img + originalsky)
noise_forsigma = np.random.poisson(img + originalsky) # I don't want to use the created noise as a noise estimate because that would produce an artificially low chi^2
img = noise - originalsky
hdu = fits.PrimaryHDU(img)
hdu.writeto('new_noise.fits',clobber=True)

# fitting with sherpa a sersic profile convolved with the psf
from astropy.modeling.models import Gaussian2D
gauss = Gaussian2D(amplitude = 1, x_mean = 26, y_mean = 26, x_stddev=3, y_stddev=3)
x_,y_ = np.meshgrid(np.arange(pixels/2.0+1), np.arange(pixels/2.0+1))
psf = gauss(x_,y_)
hdu = fits.PrimaryHDU(psf)
hdu.writeto('new_psfsherpa.fits',clobber=True)
img = convolve(sers(x,y), gauss(x_,y_), normalize_kernel=True)
noise = np.random.poisson(img + originalsky)
img = noise - originalsky
hdu = fits.PrimaryHDU(img)
hdu.writeto('new_convnoise.fits',clobber=True)
import sherpa.astro.ui as ui
ui.load_data('new_convnoise.fits')
ui.load_psf("psf","new_psfsherpa.fits")
ui.set_psf(psf)
ui.set_source(ui.const2d.c + ui.sersic2d.g)
c.c0 = 5
g.xpos=45
g.ypos=55
g.ellip = .4
g.theta=-1.5
g.n = 1.5
g.n.max=4
g.ampl=90.
ui.thaw(g)
ui.thaw(c.c0)
ui.set_method("neldermead")
ui.set_method_opt("iquad",0)
ui.set_method_opt("finalsimplex",0)
ui.fit()
#ui.conf()
sherpamodel = ui.get_model_image().y
hdu = fits.PrimaryHDU(img - sherpamodel)
hdu.writeto('new_residualsconvsherpa.fits',clobber=True)

# optimize in scipy by minimizing the chi^2 on the grid, including the convolution
start_timeoptimize = time.time()
img = convolve(Sersic2D(amplitude = 100, r_eff = 10, n=1, x_0=50, y_0=50, ellip=.5, theta=-1)(x,y), gauss(x_,y_), normalize_kernel=True)
noise = np.random.poisson(img + originalsky)
img = noise - originalsky
hdu = fits.PrimaryHDU(noise - originalsky)
hdu.writeto('new_convnoise.fits',clobber=True)
def chi2(arg):
    amplitude,r_eff,n,x_0,y_0,ellip,theta = arg
    sersmodel = Sersic2D(amplitude, r_eff, n, x_0, y_0, ellip, theta)
    model = convolve(sersmodel(x,y), gauss(x_,y_), normalize_kernel=True)
    return np.sum((img - model)**2/noise_forsigma**2)
from scipy.optimize import minimize
print "Optimizing..."
#optimal = minimize(chi2, [90, 8, 1.5, 45, 55, .3, -1.2],method='Powell')
optimal = minimize(chi2, [90, 8, 1.5, 45, 55, .3, -1.2],method='Nelder-Mead',options={'maxiter': 2000})
hdu = fits.PrimaryHDU(img - convolve(Sersic2D(optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6])(x,y), gauss(x_,y_), normalize_kernel=True))
hdu.writeto('new_residualsconv.fits',clobber=True)
print optimal
print "Chi^2/nu",np.sum((img - convolve(Sersic2D(optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6])(x,y), gauss(x_,y_), normalize_kernel=True) + originalsky)**2/noise_forsigma**2)/(pixels**2 - 7), "nu =", pixels**2 - 7
    # !!!!!!!!!!!!! careful what I did here, adding the sky level to match the noise level
print("Time for scipy optimization: --- %s seconds ---" % (time.time() - start_timeoptimize))

# emcee

import emcee
theta_true = (100, 10, 1, 50, 50, .5, -1)
# priors I use:
# amplitude [0,1000]
# r_eff > 0; Gauss(10,10)
# n (0.5,10]; Gauss(2,1)
# x_0 Gauss(50,2)
# y_0 Gauss(50,2)
# ellip [0.2,1)
# theta Gauss(-1,0.2)
def log_prior(arg):
    amplitude, r_eff, n, x_0, y_0, ellip, theta = arg
    if amplitude < 0 or amplitude > 1000 or r_eff <= 0 or n <= 0.5 or n > 10 or ellip < 0.2 or ellip >= 1:
        return -np.inf # log(0)
    return -0.5 * ((r_eff - 10)**2)/(10**2) -0.5 * ((n - 2)**2)/(0.5**2) -0.5 * ((x_0 - 50)**2)/(2**2) -0.5 * ((y_0 - 50)**2)/(2**2) -0.5 * ((ellip + 0.5)**2)/(0.15**2) -0.5 * ((theta + 1)**2)/(0.2**2)
def log_like(arg):
    amplitude, r_eff, n, x_0, y_0, ellip, theta = arg
    return -0.5 * chi2(arg)
def log_posterior(arg):
    if not np.isfinite(log_prior(arg)):
        return -np.inf
    return log_prior(arg) + log_like(arg)
ndim = 7 # number of parameters in the model
nwalkers = 20 # number of MCMC walkers
nsteps = 2000 # number of MCMC steps to take including burn-in
# initial guesses in a small ball around the best-fit
starting_guesses = np.c_[abs(np.random.normal(optimal.x[0],10,nwalkers)),abs(np.random.normal(optimal.x[1],2,nwalkers)),abs(np.random.normal(optimal.x[2],0.1,nwalkers)),abs(np.random.normal(optimal.x[3],2,nwalkers)),abs(np.random.normal(optimal.x[4],2,nwalkers)),abs(np.random.normal(optimal.x[5],0.1,nwalkers)),np.random.normal(optimal.x[6],0.1,nwalkers)] # use this in case I want to run emcee instead of parallel tempering

print("Running MCMC...")
#start_timeemcee1 = time.time()
#sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior)
# run while showing the progress
#for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations*
    #if (i+1) % 100 == 0:
        #print("{0:5.1%}".format(float(i) / nsteps))
#print("Time for emcee with 1 thread: --- %s seconds ---" % (time.time() - start_timeemcee1))

start_timeemcee8 = time.time()
sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior,threads = 7)
# run while showing the progress
for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations*
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nsteps))
print("Time for emcee with 7 threads: --- %s seconds ---" % (time.time() - start_timeemcee8))

# plot the time laps
import corner
import pylab as plt
from matplotlib.ticker import MaxNLocator
plt.clf()
fig, axes = plt.subplots(7, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(theta_true[0], color="#888888", lw=2)
axes[0].set_ylabel("amplitude")
axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(theta_true[1], color="#888888", lw=2)
axes[1].set_ylabel("r_eff")
axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(theta_true[2], color="#888888", lw=2)
axes[2].set_ylabel("n")
axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
axes[3].yaxis.set_major_locator(MaxNLocator(5))
axes[3].axhline(theta_true[3], color="#888888", lw=2)
axes[3].set_ylabel("x_0")
axes[4].plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)
axes[4].yaxis.set_major_locator(MaxNLocator(5))
axes[4].axhline(theta_true[4], color="#888888", lw=2)
axes[4].set_ylabel("y_0")
axes[5].plot(sampler.chain[:, :, 5].T, color="k", alpha=0.4)
axes[5].yaxis.set_major_locator(MaxNLocator(5))
axes[5].axhline(theta_true[5], color="#888888", lw=2)
axes[5].set_ylabel("ellip")
axes[6].plot(sampler.chain[:, :, 6].T, color="k", alpha=0.4)
axes[6].yaxis.set_major_locator(MaxNLocator(5))
axes[6].axhline(theta_true[6], color="#888888", lw=2)
axes[6].set_ylabel("theta")
axes[2].set_xlabel("step number")
fig.tight_layout(h_pad=0.0)
fig.show()

# print diagnostics and do corner plot
nburn = 1000 # "burn-in" to stabilize chains
samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
#alpha_samp = sampler.flatchain.T[0]
#beta_samp = sampler.flatchain.T[1]
#sigma_samp = sampler.flatchain.T[2]
#print("Autocorrelation time:", sampler.get_autocorr_time())
print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
print "median, std 1: ", np.median(samples[:,0]), np.std(samples[:,0])
print "median, std 2: ", np.median(samples[:,1]), np.std(samples[:,1])
print "median, std 3: ", np.median(samples[:,2]), np.std(samples[:,2])
print "median, std 4: ", np.median(samples[:,3]), np.std(samples[:,3])
print "median, std 5: ", np.median(samples[:,4]), np.std(samples[:,4])
print "median, std 6: ", np.median(samples[:,5]), np.std(samples[:,5])
print "median, std 7: ", np.median(samples[:,6]), np.std(samples[:,6])
fig = corner.corner(samples, labels=["amplitude", "r_eff", "n", "x_0", "y_0", "ellip", "theta"],truths=[theta_true[0],theta_true[1],theta_true[2],theta_true[3],theta_true[4],theta_true[5],theta_true[6]])
fig.show()

# parallel tempering; takes longer to run because it runs emcee for each temperature
start_timept = time.time()
ntemps = 5
starting_guesses = np.zeros((ntemps,nwalkers,ndim))
for i in range(ntemps):
    starting_guesses[i] = np.c_[abs(np.random.normal(optimal.x[0],10,nwalkers)),abs(np.random.normal(optimal.x[1],2,nwalkers)),abs(np.random.normal(optimal.x[2],0.1,nwalkers)),abs(np.random.normal(optimal.x[3],2,nwalkers)),abs(np.random.normal(optimal.x[4],2,nwalkers)),abs(np.random.normal(optimal.x[5],0.1,nwalkers)),np.random.normal(optimal.x[6],0.1,nwalkers)]
from emcee import PTSampler
sampler = PTSampler(ntemps,nwalkers,ndim,log_like,log_prior,threads = 7)
print("Running MCMC...")
for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations*
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nsteps))
print("Time for PT with 7 threads: --- %s seconds ---" % (time.time() - start_timept))
print "Evidence: ", sampler.thermodynamic_integration_log_evidence()
# plot the time laps
plt.clf()
fig, axes = plt.subplots(7, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, :, 0].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(theta_true[0], color="#888888", lw=2)
axes[0].set_ylabel("amplitude")
axes[1].plot(sampler.chain[:, :, :, 1].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(theta_true[1], color="#888888", lw=2)
axes[1].set_ylabel("r_eff")
axes[2].plot(sampler.chain[:, :, :, 2].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(theta_true[2], color="#888888", lw=2)
axes[2].set_ylabel("n")
axes[3].plot(sampler.chain[:, :, :, 3].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[3].yaxis.set_major_locator(MaxNLocator(5))
axes[3].axhline(theta_true[3], color="#888888", lw=2)
axes[3].set_ylabel("x_0")
axes[4].plot(sampler.chain[:, :, :, 4].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[4].yaxis.set_major_locator(MaxNLocator(5))
axes[4].axhline(theta_true[4], color="#888888", lw=2)
axes[4].set_ylabel("y_0")
axes[5].plot(sampler.chain[:, :, :, 5].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[5].yaxis.set_major_locator(MaxNLocator(5))
axes[5].axhline(theta_true[5], color="#888888", lw=2)
axes[5].set_ylabel("ellip")
axes[6].plot(sampler.chain[:, :, :, 6].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
axes[6].yaxis.set_major_locator(MaxNLocator(5))
axes[6].axhline(theta_true[6], color="#888888", lw=2)
axes[6].set_ylabel("theta")
axes[2].set_xlabel("step number")
fig.tight_layout(h_pad=0.0)
fig.show()
plt.clf()
samples = sampler.chain[:, :, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
fig = corner.corner(samples, labels=["amplitude", "r_eff", "n", "x_0", "y_0", "ellip", "theta"],truths=[theta_true[0],theta_true[1],theta_true[2],theta_true[3],theta_true[4],theta_true[5],theta_true[6]])
fig.show()
emcee.autocorr.integrated_time(samples,axis=0)
emcee.autocorr.function(samples,axis=0)

#Time for scipy optimization: --- 52.5281159878 seconds ---
#Time for emcee with 1 thread: --- 323.612744093 seconds ---
#Time for emcee with 4 threads: --- 132.817188978 seconds ---
#Time for emcee with 7 threads: --- 118.459464073 seconds ---
#Time for PT with 7 threads: --- 452.805958033 seconds ---

# FUTURE WORKS
'''
1) save progress, run on multiple computers
1) PT and emcee give different distributions even after setting more stringent prior on ellipticity
1) why does it say AutocorrError: The chain is too short to reliably estimate the autocorrelation time ?
1) careful because right now all images with noise have integer pixel values
1) convolution central pixels: http://docs.astropy.org/en/stable/api/astropy.convolution.discretize_model.html#astropy.convolution.discretize_model oversample for sub-pixel; pyprofit seems the better option
3) make physical units: exptime, gain, scale, mag
4) lens it
5) mask central pixels
'''
