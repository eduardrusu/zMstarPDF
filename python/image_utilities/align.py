# Subpixel cross-correlation alignment from http://image-registration.readthedocs.io/en/latest/#quick-example

from astropy.io import fits
from image_registration import chi2_shift
from image_registration.fft_tools import shift

image0 = fits.open('VISTAmatch_Ks_small.fits')
data0 = image0[0].data

# repeat this block for each image
image1 = fits.open('VISTAmatch_J_small.fits')
data1 = image1[0].data
xoff, yoff, exoff, eyoff = chi2_shift(data0, data1,err=None,return_error=True,upsample_factor='auto')
data1_corr = shift.shiftnd(data1, (-yoff, -xoff))
image1_corr = image1
image1_corr[0].data = data1_corr
image1_corr.writeto('VISTAmatch_J_small_shift.fits',overwrite=True)

image1 = fits.open('VISTAmatch_Y_small.fits')
data1 = image1[0].data
xoff, yoff, exoff, eyoff = chi2_shift(data0, data1,err=None,return_error=True,upsample_factor='auto')
data1_corr = shift.shiftnd(data1, (-yoff, -xoff))
image1_corr = image1
image1_corr[0].data = data1_corr
image1_corr.writeto('VISTAmatch_Y_small_shift.fits',overwrite=True)
