##########################
# Simple code for scatter plot with y error bars. Computes bias, scatter and fraction of outliers. Ignores the outliers when computing bias and scatter
##########################

from matplotlib.colors import LogNorm
import scipy.optimize as optimization
from pylab import *
import numpy as np
ax=plt.subplot(111)
zlim = 1.0
outlim = 0.15

fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)
#ax1.set(aspect=1)
ax1.set_aspect(1, adjustable='datalim')

ax = plt.subplot(1,1,1, sharex=ax1, sharey=ax1)
#x, y, yinf, ysup = np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz_specflag0_bpz.cat", usecols=(9, 1, 2, 3), unpack=True)
x, y, yinf, ysup = np.loadtxt("bpzspeczeazy.cat", usecols=(3, 0, 1, 2), unpack=True)
#x = x[abs(y) <= zlim]
#y = y[abs(y) <= zlim]
#yinf = yinf[abs(y) <= zlim]
#ysup = ysup[abs(y) <= zlim]
#y = y[abs(x) <= zlim]
#x = x[abs(x) <= zlim]
#yinf = yinf[abs(x) <= zlim]
#ysup = ysup[abs(x) <= zlim]
ax.tick_params(labelsize=14)
#plt.scatter(x,y, color='k')
plt.errorbar(x, y, yerr=[y-yinf, ysup-y], fmt='o')
delta = (y-x)/(1+x)
#delta = (y-(m+x))/(1+x)
firstpass_std = np.std(delta)
#std = np.std(delta[abs(delta)<(4*firstpass_std)])
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % len(x)
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
#ax.plot(x, m*x + b, '--')
ax.plot(x, m+x, 'b--')
ax.plot(x, x, 'g')
plt.plot(x, 0.85*x-0.15, 'r--')
plt.plot(x, 1.15*x+0.15, 'r--')
ax.text(0.05, 0.95, stdoutobj, fontsize=15, color='black',transform=ax.transAxes)
ax.text(0.05, 0.90, stdout, fontsize=15, color='black',transform=ax.transAxes)
ax.text(0.05, 0.85, outlier, fontsize=15, color='black',transform=ax.transAxes)
ax.text(0.05, 0.80, bias, fontsize=15, color='black',transform=ax.transAxes)
#ax.text(0.05, 0.80, bestfit, fontsize=15, color='black',transform=ax.transAxes)
plt.xlabel('spec-z')
#plt.ylabel('phot-z')
plt.ylabel('bpz')
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(0, zlim)
plt.ylim(0, zlim)
ax1.set_yticklabels(np.arange(0.0, 1.0, 0.2))
ax1.set_xticklabels(np.arange(0.0, 1.0, 0.2))
#plt.tick_params(labelbottom='off')
#plt.title('HE0435 ugri specz-photz')

plt.subplots_adjust(bottom=0.1, left =0.2, right=0.9, top=0.90, wspace=0, hspace=0)
#plt.tight_layout()
#fig.text(0.05, 0.5, 'photo-z', ha='center', va='center', size='20', rotation='vertical')
#plt.title('HE0435 ugri specz-photz')
plt.savefig('specz-bpz.png')
