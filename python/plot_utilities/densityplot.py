# produces density plot for the MS catalogue Mhalo vs catalogue Mstar

from matplotlib.colors import LogNorm
import scipy.optimize as optimization
from pylab import *
import numpy as np

font = 10
ticksize = 10

plt.clf()
fig = plt.figure(figsize=(10,12))
#fig, axes = plt.subplots(nrows=2, ncols=2)

ax1 = fig.add_subplot(1,1,1)
ax1.set_aspect(1)
x, y = np.loadtxt("/Users/perseus/Desktop/GGL_los_8_0_0_0_0_N_4096_ang_4_Bower_galaxies_on_plane_27_to_63.imagesplot.txt", usecols=(0, 2), unpack=True)
zlim = 14.5
zlim_ = 10
#x = np.log10(x)
#y = np.log10(y)
x = x[(abs(y) <= zlim) & (abs(y) > zlim_)]
y = y[(abs(y) <= zlim) & (abs(y) > zlim_)]
y = y[(abs(x) <= zlim) & (abs(x) > zlim_)]
x = x[(abs(x) <= zlim) & (abs(x) > zlim_)]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
colorbar()
delta = (y-x)/(1+x)
plt.xlabel('Mhalo', fontsize=font)
plt.ylabel('Mhalo (Behroozi)', fontsize=font)
plt.xlim(zlim_, zlim)
plt.ylim(zlim_, zlim)
plt.savefig('/Users/perseus/Desktop/Bower.png' , dpi=250)
