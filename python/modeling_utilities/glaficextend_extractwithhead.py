# Header manupulation

from astropy.io import fits

image0 = fits.open('../../F125W_visit1_highres_drz_sci_cut.fits')
data0 = image0[0].data
head0 = image0[0].header

image1 = fits.open('glafic125_1_extendSIEgamma_image.fits')
data1 = image1[0].data[5]

image2 = image0
image2[0].header = head0
image2[0].data = data0 - data1
image2.writeto('glafic125_1_extendSIEgamma_image_subtract.fits',overwrite=True)
