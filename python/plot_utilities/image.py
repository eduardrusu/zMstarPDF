# creates figure 1 from Rusu et al. 2017

from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
plt.clf()
pix = 915 # number of pixels in 4 arcmin
scale = 240.0 / pix
maglim = 23
x,y,i,classify = np.loadtxt("/Users/eduardrusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat",usecols=[0,1,4,97],unpack=True)
sep = np.sqrt((x - 457.5)**2 + (y - 457.5)**2) * scale

x = x[(sep <= 120) & (i <= maglim)]
y = y[(sep <= 120) & (i <= maglim)]
classify = classify[(sep <= 120) & (i <= maglim)]
sep = sep[(sep <= 120) & (i <= maglim)]

x_spec = x[(classify == 0) | (classify == 1)]
y_spec = y[(classify == 0) | (classify == 1)]

x_star = x[classify < 0]
y_star = y[classify < 0]

x_galnospec = x[classify == 2]
y_galnospec = y[classify == 2]

image = fits.getdata("/Users/eduardrusu/Desktop/WFI2033/WFI2033analysis/FINALweighted_ir.fits")
image[image<0]=0.0001

mask = fits.getdata("6_120arcsec.fits")

plt.scatter(x_star,y_star,marker='*',facecolors='none', edgecolors='k')
plt.scatter(x_galnospec,y_galnospec,marker='o',facecolors='none', edgecolors='k')
plt.scatter(x_spec,y_spec,marker='s',facecolors='none', edgecolors='k')
plt.imshow(image, cmap='gray_r', norm=LogNorm(), origin='lower', vmin=0.001, vmax=100)
plt.imshow(mask, cmap='Oranges', origin='lower', alpha=0.2)
circle1 = plt.Circle((pix/2.0,pix/2.0),45/120.0*pix/2.0,color='k',fill=False)
circle2 = plt.Circle((pix/2.0,pix/2.0),pix/2.0,color='k',fill=False)
fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)
fig = plt.gca()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig('FOV_WFI2033.png', dpi=300, bbox_inches='tight')
