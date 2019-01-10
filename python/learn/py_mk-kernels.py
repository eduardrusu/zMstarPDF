# Code from Adam Tomczak, used to create kernels in order to broaden a narrow PSF in order to match a broader one. Approximates the PSF analytically.
# The actual convolution of images using the kernel is done by the convolution.py code.

import pylab
from scipy import signal
from scipy import ndimage
from scipy import optimize
from astropy.io import fits
from astropy import modeling
from skimage import restoration
from photutils import CircularAperture
from photutils import aperture_photometry

moffat2d = modeling.functional_models.Moffat2D()

class MoffatTomczak2D(modeling.Fittable2DModel):
    amplitude = modeling.Parameter(default=1)
    x_0 = modeling.Parameter(default=0)
    y_0 = modeling.Parameter(default=0)
    scale = modeling.Parameter(default=1)
    rotation = modeling.Parameter(default=0)
    gamma = modeling.Parameter(default=1)
    alpha = modeling.Parameter(default=1)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, scale, rotation, gamma, alpha):
		rr_gg = (((x - x_0)*scale) ** 2 + (y - y_0) ** 2) / gamma ** 2
		rr_gg_alpha = amplitude * (1 + rr_gg) ** (-alpha)
		rr_gg_alpha_rotate = ndimage.rotate(rr_gg_alpha, rotation, reshape=0)
		return rr_gg_alpha_rotate

moffatTomczak = MoffatTomczak2D()

model_fitter = modeling.fitting.LevMarLSQFitter()

lens_name = 'WFI2033'

#imnames = ['u', 'i']

#psfs = [fits.getdata('u_psfnorm.fits'),
#		fits.getdata('g_psfnorm.fits')]

imnames = ['r', 'ch1']

psfs = [fits.getdata('modifiedforTPHOT/ir_PSFfix.fits'),
        fits.getdata('ch1/ch1psffix_oversampled.fits')]



cogs = []
cogs_moffs = []
fwhms = []
fluxtot = pylab.array([])
profiles = []
moffat_psfs = []
moffatTomczak_psfs = []

bgs = pylab.array([])
scatters = pylab.array([])

cogs_conv = []
fwhms_conv = []

dy, dx = psfs[0].shape
y0, x0 = dy/2, dx/2

px_scale = 0.2018
radii = pylab.linspace(0.5, dx/2-1, 20)
apers = [CircularAperture([x0, y0], r=ri) for ri in radii]



fig = pylab.figure(figsize=(12.4, 12.9))
sp2 = fig.add_subplot(222)
sp1 = fig.add_subplot(221)
sp4 = fig.add_subplot(224)
sp3 = fig.add_subplot(223)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp4.minorticks_on()

sp1.grid()
sp2.grid()
sp3.grid()
sp4.grid()

sp1.set_xlabel('radius ["]')
sp2.set_xlabel('radius ["]')	
sp3.set_xlabel('radius ["]')
sp4.set_xlabel('radius ["]')

sp1.set_ylabel('Flux(r)')
sp3.set_ylabel('$\Sigma$ Flux(< r) / model')

#sp3.axhline(1, color='k', lw=1.5, ls='--')
#sp4.axhline(1, color='k', lw=1.5, ls='--')

fig.subplots_adjust(wspace=0)


for i in range(len(psfs)):

	psf = psfs[i]


	###  Estimating FWHM
	xx, yy = pylab.meshgrid(range(-dx/2+1, dx/2+1), range(-dx/2+1, dx/2+1))
	rr = (xx**2 + yy**2)**0.5
	r_uniq = pylab.unique(rr)
	profile = pylab.array([])
	for ri in r_uniq:
		rinds = pylab.where(rr == ri)
		profile = pylab.append(profile, pylab.mean(psf[rinds]))
	fwhm = 2 * pylab.interp(profile.max()/2., profile[::-1], r_uniq[::-1])


	###  Estimating background, noise, and total flux of psf
	inds_outer = pylab.where(rr > 2*fwhm)
	bg = pylab.median(psf[inds_outer])
	scatter = pylab.std(psf[inds_outer])
	bgs = pylab.append(bgs, bg)
	scatters = pylab.append(scatters, scatter)
	inds_seg = pylab.where(psf > bg + 2*scatter)
	inds_nonseg = pylab.where(psf <= bg + 2*scatter)

#	psf -= bg
	aper_tot = CircularAperture([x0, y0], r=4./0.2)
	phot_tot = aperture_photometry(psf, aper_tot)
	fluxtot = pylab.append(fluxtot, phot_tot['aperture_sum'])
	psf /= fluxtot[-1]


	###  Re-estimating FWHM from bg-subtracted psf
	rr = (xx**2 + yy**2)**0.5
	r_uniq = pylab.unique(rr)
	profile = pylab.array([])
	for ri in r_uniq:
		rinds = pylab.where(rr == ri)
		profile = pylab.append(profile, pylab.mean(psf[rinds]))
	fwhm = 2 * pylab.interp(profile.max()/2., profile[::-1], r_uniq[::-1])
	fwhms.append(fwhm)


	###  Fitting a moffat profile
	fit = model_fitter(moffat2d, xx, yy, psf)
	fits.writeto('moffat_' + imnames[i] + '.fits', fit(xx, yy), clobber=1)
	fits.writeto('resid_' + imnames[i] + '-moffat.fits', psf-fit(xx, yy), clobber=1)

	###  Fitting a moffat+tomczak profile
	fit2 = model_fitter(moffatTomczak, xx, yy, psf)
	fits.writeto('moffatTomczak_' + imnames[i] + '.fits', fit2(xx, yy), clobber=1)
	fits.writeto('resid_' + imnames[i] + '-moffatTomczak.fits', psf-fit2(xx, yy), clobber=1)

	moffat_psfs.append(fit(xx, yy))
	moffatTomczak_psfs.append(fit(xx, yy))



	###  Plotting profiles
	color = pylab.cm.jet(i * 1. / (len(psfs) - 1))
	sp1.plot(r_uniq * px_scale, profile / profile.max(), color=color, lw=3, label=imnames[i] + '  %.2f' % (fwhm * px_scale) + '"')
	profiles.append(profile / profile.max())

	sp2.axhline(10, color=color, label=imnames[i] + '  %.2f' % (fwhm * px_scale) + '"')


	###  Measuring curve of growth
	cog = pylab.array([])
	cog_moff = pylab.array([])
	for aper in apers:
		phot = aperture_photometry(psf, aper)
		phot_moff = aperture_photometry(fit(xx, yy), aper)
		cog = pylab.append(cog, phot['aperture_sum'])
		cog_moff = pylab.append(cog_moff, phot_moff['aperture_sum'])
	cogs.append(cog)
	cogs_moffs.append(cog_moff)


	print imnames[i], '%.2f' % fwhm


###  PSF with largest FWHM
ind_max_fwhm = fwhms.index(max(fwhms))
color = pylab.cm.jet(ind_max_fwhm * 1. / (len(psfs) - 1))
sp2.plot(r_uniq * px_scale, profiles[ind_max_fwhm], color=color, lw=3)





sp2.legend(loc=1, title=lens_name)
sp1.axis([-0.1, 4.4, -0.1, 1.1])
sp2.axis([-0.1, 4.4, -0.1, 1.1])






###  Plotting curves of growth
cogs = pylab.array(cogs)
mean_cog = pylab.average(cogs, axis=0)

for i in range(len(psfs)):
	sp3.plot(radii * px_scale, cogs[i] / cogs_moffs[ind_max_fwhm], color=pylab.cm.jet(i * 1. / (len(psfs) - 1)))



sp3.axis([-0.1, 4.4, 0.87, 1.13])
sp4.axis([-0.1, 4.4, 0.87, 1.13])






for i in range(len(fwhms)):

	if i == ind_max_fwhm:
		target_cog = cogs[i]
		cogs_conv.append(cogs[i])
		continue

	psf1 = psfs[i]
	psf2 = psfs[ind_max_fwhm]

	moff1 = moffat_psfs[i]
	moff2 = moffat_psfs[ind_max_fwhm]

	moffTomcz1 = moffatTomczak_psfs[i]
	moffTomcz2 = moffatTomczak_psfs[ind_max_fwhm]

	###  Masking out corners
	inds_corners = pylab.where(rr > dx/2)
	moff1[inds_corners] *= 0
	moff2[inds_corners] *= 0
	moffTomcz1[inds_corners] *= 0
	moffTomcz2[inds_corners] *= 0

	###  Creating hybrid PSFs by replacing pixes <(bg+3sigma) with the model
	indslo1 = pylab.where(psf1 < bgs[i] + 3*scatters[i])
	indslo2 = pylab.where(psf1 < bgs[ind_max_fwhm] + 3*scatters[ind_max_fwhm])

	hybrid1, hybrid2 = psf1 * 1., psf2 * 1.
	hybrid1[indslo1] = moffTomcz1[indslo1]
	hybrid2[indslo2] = moffTomcz2[indslo2]



	niter = 100

	'''
	###  Data kernel
	kernel_data = restoration.richardson_lucy(psf2, psf1, iterations=niter)
	kernel_data /= kernel_data.sum()
	k0 = 'kernel_' + imnames[ind_max_fwhm] + '-' + imnames[i] + '_data.fits'
	fits.writeto(k0, kernel_data, clobber=1)

	###  Moffat kernel
	kernel_moff = restoration.richardson_lucy(moff2, moff1, iterations=niter)
	kernel_moff /= kernel_moff.sum()
	k1 = 'kernel_' + imnames[ind_max_fwhm] + '-' + imnames[i] + '_moff.fits'
	fits.writeto(k1, kernel_moff, clobber=1)

	###  Moffat-Tomczak kernel
	kernel_moffTomcz = restoration.richardson_lucy(moffTomcz2, moffTomcz1, iterations=niter)
	kernel_moffTomcz /= kernel_moffTomcz.sum()
	k2 = 'kernel_' + imnames[ind_max_fwhm] + '-' + imnames[i] + '_moffTomcz.fits'
	fits.writeto(k2, kernel_moffTomcz, clobber=1)
	'''

	###  Hybrid kernel
	kernel_hybrid = restoration.richardson_lucy(hybrid2, hybrid1, iterations=niter)
	kernel_hybrid /= kernel_hybrid.sum()
	k3 = 'kernel_' + imnames[ind_max_fwhm] + '-' + imnames[i] + '_hybrid.fits'
	fits.writeto(k3, kernel_hybrid, clobber=1)



	###  Adding info to kernel headers
	for kern in [k3]:
		f = fits.open(kern, mode='update')
		h = f[0].header
		h['NITER'] = (niter, 'number of Richardson Lucy iterations')
		f.close()



	###  Testing by convolving psf1 by kernel
	conv = signal.convolve2d(psf1, kernel_hybrid, mode='same')

	fits.writeto('conv_'+ imnames[ind_max_fwhm] + '-' + imnames[i] + '.fits', conv, clobber=1)
	fits.writeto('resid_'+ imnames[ind_max_fwhm] + '-' + imnames[i] + '.fits', psf2 - conv, clobber=1)

	print 'generated kernels from', imnames[i], 'to', imnames[ind_max_fwhm]




	###  Estimating convolved FWHM
	profile = pylab.array([])
	for ri in r_uniq:
		rinds = pylab.where(rr == ri)
		profile = pylab.append(profile, pylab.mean(conv[rinds]))
	fwhm = 2 * pylab.interp(profile.max()/2., profile[::-1], r_uniq[::-1])
	fwhms_conv.append(fwhm)


	###  Plotting profiles
	color = pylab.cm.jet(i * 1. / (len(psfs) - 1))
	sp2.plot(r_uniq * px_scale, profile / profile.max(), color=color, lw=3)


	###  Measuring curve of growth
	cog = pylab.array([])
	for aper in apers:
		phot = aperture_photometry(conv, aper)
		cog = pylab.append(cog, phot['aperture_sum'])
	cogs_conv.append(cog)







###  Plotting curves of growth
cogs_conv = pylab.array(cogs_conv)
mean_cog_conv = pylab.average(cogs_conv, axis=0)

for i in range(len(psfs)):
	sp4.plot(radii * px_scale, cogs_conv[i] / cogs[ind_max_fwhm], color=pylab.cm.jet(i * 1. / (len(psfs) - 1)))



###  Plotting estimate for catalog aperture as maximum FWHM +30%
sp4.axvline(max(fwhms) * px_scale, color='k', ls='--', lw=3, label='FWHM')
sp4.axvline(max(fwhms) * px_scale * 1.3, color='k', lw=3, label='"ideal"\naperture')

sp4.legend(loc=1, prop={'size':18})
fig.savefig('psfdiagnostics.png' , dpi=250)

