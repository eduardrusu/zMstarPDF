# used to plot Fig. 1 from Rusu et al. 2017

from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
plt.clf()
star = np.loadtxt("star.cat",unpack=True)
specz = np.loadtxt("specz.cat",unpack=True)
photz = np.loadtxt("photz.cat",unpack=True)
image = fits.getdata("0435_scam_mar14_r_blk.fits")
image[image<0]=0.0001
mask = fits.getdata("mskHE0435_asecrad120_no5arcsec_blk.fits")
mask[mask!=0]=1
plt.scatter(star[0]/2,star[1]/2,marker='*',facecolors='none', edgecolors='k')
plt.scatter(photz[0]/2,photz[1]/2,marker='o',facecolors='none', edgecolors='k')
plt.scatter(specz[0]/2,specz[1]/2,marker='s',facecolors='none', edgecolors='k')
plt.imshow(image, cmap='gray_r', norm=LogNorm(), origin='lower', vmin=0.001, vmax=100)
plt.imshow(mask, cmap='Oranges', origin='lower', alpha=0.2)
circle1 = plt.Circle((300,300),225/2.0,color='k',fill=False)
circle2 = plt.Circle((300,300),300,color='k',fill=False)
fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)
fig = plt.gca()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig('FOV_small.png', dpi=300, bbox_inches='tight')
