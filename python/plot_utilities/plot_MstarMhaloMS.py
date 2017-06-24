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
#ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images_forNAOJ.txt" % (i,j,k), usecols=(4, 5), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 15
outlim = 0.15
x = np.log10(x)
y = np.log10(y)
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
colorbar()
delta = (y-x)/(1+x)
plt.xlabel('Mhalo', fontsize=font)
plt.ylabel('Mstar', fontsize=font)
plt.xlim(8, 16)
plt.ylim(5, 12)
plt.savefig('/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/MstarMhaloMS.png' , dpi=250)
