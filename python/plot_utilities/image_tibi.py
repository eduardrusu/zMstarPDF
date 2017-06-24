# creates the lens montage in Berghea et al. 2017

from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

plt.clf()

image_g = fits.getdata("/Users/eduardrusu/GITHUB/tibi/g_cut.fits")
image_r = fits.getdata("/Users/eduardrusu/GITHUB/tibi/r_cut.fits")
image_i = fits.getdata("/Users/eduardrusu/GITHUB/tibi/i_cut.fits")
image_z = fits.getdata("/Users/eduardrusu/GITHUB/tibi/z_cut.fits")
image_y = fits.getdata("/Users/eduardrusu/GITHUB/tibi/y_cut.fits")
image_g_sub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/g_sub.fits")
image_r_sub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/r_sub.fits")
image_i_sub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/i_sub.fits")
image_z_sub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/z_sub.fits")
image_y_sub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/y_sub.fits")
image_i_lenssub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/i_lenssub.fits")
image_z_lenssub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/z_lenssub.fits")
image_y_lenssub = fits.getdata("/Users/eduardrusu/GITHUB/tibi/y_lenssub.fits")
fig = plt.figure(figsize=(5,3))
ax1 = fig.add_subplot(1,1,1)
ax1.set_aspect(1)

for i in range(15):
        ax1 = plt.subplot(3,5,i+1, sharex=ax1, sharey=ax1)
        ax1.axes.get_xaxis().set_visible(False)
        ax1.axes.get_yaxis().set_visible(False)
        if i == 0: plt.imshow(image_g, cmap='gray_r', origin='lower')
        if i == 1: plt.imshow(image_r, cmap='gray_r', origin='lower')
        if i == 2: plt.imshow(image_i, cmap='gray_r', origin='lower')
        if i == 3: plt.imshow(image_z, cmap='gray_r', origin='lower')
        if i == 4: plt.imshow(image_y, cmap='gray_r', origin='lower')
        if i == 5: plt.imshow(image_g_sub, cmap='gray_r', origin='lower')
        if i == 6: plt.imshow(image_r_sub, cmap='gray_r', origin='lower')
        if i == 7: plt.imshow(image_i_sub, cmap='gray_r', origin='lower')
        if i == 8: plt.imshow(image_z_sub, cmap='gray_r', origin='lower')
        if i == 9: plt.imshow(image_y_sub, cmap='gray_r', origin='lower')
        if i == 10: plt.imshow(image_g_sub, cmap='gray_r', origin='lower')
        if i == 11: plt.imshow(image_r_sub, cmap='gray_r', origin='lower')
        if i == 12: plt.imshow(image_i_lenssub, cmap='gray_r', origin='lower')
        if i == 13: plt.imshow(image_z_lenssub, cmap='gray_r', origin='lower')
        if i == 14: plt.imshow(image_y_lenssub, cmap='gray_r', origin='lower')

#if i == 14: plt.imshow(image_y_lenssub, cmap='gray_r', norm=LogNorm(), origin='lower', vmin=1, vmax=10000000)
plt.subplots_adjust(bottom=0, left=0, right=1, top=1, wspace=0, hspace=0)
plt.savefig('hostlens.eps', dpi=300, bbox_inches='tight')
