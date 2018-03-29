# Header manupulation

from astropy.io import fits

image0 = fits.open('VISTAmatch_Ks_small.fits')
data0 = image0[0].data
head0 = image0[0].header

image1 = fits.open('lensfixcompanion_correctPSFshowlens_subtractK.fits')

#h['NITER'] = (niter, 'number of Richardson Lucy iterations')

image1_corr = image1
image1_corr[0].header = head0
image1_corr.writeto('lensfixcompanion_correctPSFshowlens_subtractK.fits',overwrite=True)
