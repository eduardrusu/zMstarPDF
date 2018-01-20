# creates the lens montage in Rusu et al. 2018a

from astropy.io import fits
#from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from astropy.convolution import Gaussian2DKernel, convolve

plt.clf()

lim = 2 # sigma limits for plot range
kernel = Gaussian2DKernel(stddev=2)

# the simple way to center the images properly is to do a match > wcs in ds9 and record the matching central pixel coordinates in each image; python does not allow to non-integer margins, so do as below:
center_g=[71,35]
center_r=[71,35]
center_i=[71,35]
center_z=[71,35]
center_y=[71,35]
center_Y=[54,26]
center_J=[54,27]
center_K=[54,26]
side_g = 60 # 15 arcsec
side_r = 60
side_i = 60
side_z = 60
side_y = 60
side_Y = 44
side_J = 44
side_K = 44

image_g = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_g/nolens_analPSF_subtract.fits")
image_r = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_r/nolens_correctPSF_subtract.fits")
image_i = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_i/nolens_correctPSF_subtract.fits")
image_z = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_z/nolens_correctPSF_subtract.fits")
image_y = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_y/nolens_analPSF_subtract.fits")
image_Y = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_VISTAY/nolens_analPSF_subtract.fits")
image_J = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_J/lens_correctPSF_subtract.fits")
image_K = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/final_K/lens_correctPSF_subtract.fits")

noise_g = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/g_sigma_small.fits")
noise_r = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/r_sigma_small.fits")
noise_i = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/i_sigma_small.fits")
noise_z = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/z_sigma_small.fits")
noise_y = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/y_sigma_small.fits")
noise_Y = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/VISTAmatch_Y_sigma_small.fits")
noise_J = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/VISTAmatch_J_sigma_small.fits")
noise_K = fits.getdata("/Users/cerusu/OneDrive - Subaru Telescope/2M1134-2103/images/VISTAmatch_Ks_sigma_small.fits")

cat_g = np.loadtxt("astrom_g.cat",unpack=True)
cat_r = np.loadtxt("astrom_r.cat",unpack=True)
cat_i = np.loadtxt("astrom_i.cat",unpack=True)
cat_z = np.loadtxt("astrom_z.cat",unpack=True)
cat_y = np.loadtxt("astrom_y.cat",unpack=True)
cat_Y = np.loadtxt("astrom_YVISTA.cat",unpack=True)
cat_J = np.loadtxt("astrom_J.cat",unpack=True)
cat_K = np.loadtxt("astrom_K.cat",unpack=True)

# I'm adding an extra column, in which I will load the colormap. Otherwise either there will be too much blank space or the color map will overlap the images. An alternative is to use ImageGrid, but that only works for images with exactly the same pixel size, which is not the case here.
fig, axes = plt.subplots(nrows=2, ncols=5, gridspec_kw={"width_ratios":[1,1,1,1,0.05]},figsize=(4,2)) # setting the colorbar column to be narrow; each of the 2x5 subplots gets an axis
# Unfortunately I couldn't find a way to make the margins of the plots touch perfectly. That is an interplay between figsize and subplots_adjust parameters

def imageplt(image,noise,center,side,cat,band,ax):
    #print np.std(image)
    image = image / noise # rescale pixel values so that the standard deviation is 1; this is in order to share one colorbar across images with different noise properties, and also so that I can have the colorbar in units of sigma
    image = convolve(image, kernel)
    global im
    im = ax.imshow(image[center[1]-int(side/2):center[1]+int(side/2),center[0]-int(side/2):center[0]+int(side/2)], cmap='gray', origin='lower', aspect='equal', vmin=-lim, vmax=lim)
    ax.text(side * 4.5/6, side * 5.0/6, band, size='10')
    ax.scatter(cat[0]-center[0]+int(side/2),cat[1]-center[1]+int(side/2),marker='*',s=5,facecolors='none', edgecolors='k')
    ax.set_xticks([])
    ax.set_yticks([])

i = 0
for ax in axes.flat:
        if i == 0: imageplt(image_g,noise_g,center_g,side_g,cat_g,'g',ax)
        if i == 1: imageplt(image_r,noise_r,center_r,side_r,cat_r,'r',ax)
        if i == 2: imageplt(image_i,noise_i,center_i,side_i,cat_i,'i',ax)
        if i == 3: imageplt(image_z,noise_z,center_z,side_z,cat_z,'z',ax)
        if i == 4: ax.axis('off')
        if i == 5: imageplt(image_y,noise_y,center_y,side_y,cat_y,'y',ax)
        if i == 6: imageplt(image_Y,noise_Y,center_Y,side_Y,cat_Y,'Y',ax)
        if i == 7: imageplt(image_J,noise_J,center_J,side_J,cat_J,'J',ax)
        if i == 8: imageplt(image_K,noise_K,center_K,side_K,cat_K,'Ks',ax)
        i += 1

ip = InsetPosition(axes.flat[8], [1.05,0,0.05,2]) # allows me to set the position of the colorbar wherever I want, in this case right after the final subplot, extended across two rows
axes.flat[9].set_axes_locator(ip)

cbar=fig.colorbar(im, cax=axes.flat[9], ax=[axes.flat[0],axes.flat[1],axes.flat[2],axes.flat[3],axes.flat[5],axes.flat[6],axes.flat[7],axes.flat[8]],ticks=np.linspace(-lim,lim,2*lim+1))
cbar.ax.tick_params(labelsize=6)
#cbar.ax.set_yticklabels([str(x) for x in np.linspace(-lim,lim,2*lim+1)])

plt.subplots_adjust(bottom=0, left=0, right=1, top=1, wspace=0, hspace=-0.02)
plt.savefig('hostlens.eps', dpi=150, bbox_inches='tight')
