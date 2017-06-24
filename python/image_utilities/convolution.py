# Code to convolve an image with a kernel, such as produced by py_mk-kernels.py

from astropy.io import fits
from scipy import signal

im=fits.open("FINALweighted_u_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-u_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALweighted_u_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALweighted_r_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-r_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALweighted_r_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALweighted_i_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-i_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALweighted_i_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALweighted_z_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-z_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALweighted_z_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALweighted_Y_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-Y_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALweighted_Y_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALHEADmedian_J_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-J_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALHEADmedian_J_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALHEADmedian_H_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-H_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALHEADmedian_H_covernolens_matchg.fits', clobber=0)

im=fits.open("FINALHEADmedian_Ks_covernolens.fits")
kernel_hybrid=fits.open("kernel_g-Ks_hybrid.fits")
conv = signal.convolve2d(im[0].data, kernel_hybrid[0].data, mode='same')
im[0].data=conv
im.writeto('FINALHEADmedian_Ks_covernolens_matchg.fits', clobber=0)

