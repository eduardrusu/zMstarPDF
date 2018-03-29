# Plots on top of a png image

import matplotlib.image as mpimg
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
plt.clf()
#plt.axes().set_aspect('equal')
img=mpimg.imread('J1206gri.png')
plt.imshow(img)
#circle1 = plt.Circle((1/2.0,1/2.0),1/2.0,color='k',fill=False)
circle1 = plt.Circle((738,736.6),736,color='k',fill=False)
gal=np.loadtxt("J1206_i23_120gal.cat",usecols=[0,1],unpack=True)
plt.scatter(1.005*gal[0],0.99*(1488-gal[1]),marker='s',facecolors='none', edgecolors='r',s=30)
fig = plt.gcf()
fig.gca().add_artist(circle1)
fig = plt.gca()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig('test.png', dpi=300, bbox_inches='tight')

