##########################
# Given a list of spectroscopic redshifts, plots the histograms and marks the groups previously identified
##########################

#from matplotlib.colors import LogNorm
#import scipy.optimize as optimization
from pylab import *
import numpy as np
#ax=plt.subplot(111)
zlim = 1.8
#maglimit = 21
#outlim = 0.15

z = np.loadtxt("/Users/cerusu/Dropbox/Davis_work/code/J1206/spec/J1206_specfromChrisandSDSSusethis.tab", usecols=[4], unpack=True)
plt.hist(z[z < zlim],bins = 72)

#x = x[abs(y) <= zlim]
#y = y[abs(y) <= zlim]
#plt.scatter(x,y, color='k')
#plt.scatter(x, y)
#stdoutobj = "%d objects" % len(x)
#stdout = "scatter = %.3f" % std
#ax.plot(x, m+x, 'b--')
#ax.plot(x, x, 'g')
#plt.plot(x, 0.85*x-0.15, 'r--')
#plt.plot(x, 1.15*x+0.15, 'r--')
#ax.text(0.05, 0.95, stdoutobj, fontsize=15, color='black',transform=ax.transAxes)
plt.xlabel('spectroscopic redshift')
#plt.ylabel('EaZy')
plt.xlim(0, zlim)
#plt.ylim(0, zlim)
#ax1.set_yticklabels(np.arange(0.0, zlim, 0.2))
#ax1.set_xticklabels(np.arange(0.0, zlim, 0.2))
#plt.subplots_adjust(bottom=0.1, left =0.2, right=0.9, top=0.90, wspace=0, hspace=0)
#plt.tight_layout()
#fig.text(0.05, 0.5, 'photo-z', ha='center', va='center', size='20', rotation='vertical')
#plt.title('HE0435 ugri specz-photz')
plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/J1206/speczhist.png')
