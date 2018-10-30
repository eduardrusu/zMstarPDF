# The code uses the output statistics produced by plotkappacompletestatistics.py/plotkappabiascompletestatistics.py in order to plot bars. Run without arguments. Make sure the uncomment the appropriate ax.set_ylim, ylabel and savefig lines

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np
import sys
import os
from os import system

root = "/Users/eduardrusu/Dropbox/Davis_work/code/GOODCODE/WFI2033kappa/"

#kappastat_45 = np.loadtxt("%skappastatistics_WFI2033_5_23_45.lst" % root, unpack=True)
#kappastat_120 = np.loadtxt("%skappastatistics_WFI2033_5_23_120.lst" % root, unpack=True)
kappastat_45 = np.loadtxt("%skappacomputebias_WFI2033_5_23_45.lst" % root, unpack=True)
kappastat_120 = np.loadtxt("%skappacomputebias_WFI2033_5_23_120.lst" % root, unpack=True)

N = 18
ind = 2.5 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

ax = plt.subplot(2,1,1)

col1 = (kappastat_45[0])
rects1 = ax.bar(ind + width, col1, width, color='r')
col2 = (kappastat_120[0])
rects2 = ax.bar(ind + 2*width, col2, width, color='b')

#ax.set_ylim([0.00,0.05])
ax.set_ylim([-0.1,0.1])
#ax.set_ylabel('median$_\kappa$')
ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')

ax = plt.subplot(2,1,2)

col3 = (kappastat_45[2])
rects3 = ax.bar(ind + width, col3, width, color='r')
col4 = (kappastat_120[2])
rects4 = ax.bar(ind + 2*width, col4, width, color='b')

ax.set_ylim([0,0.05])
#ax.set_ylabel('$\sigma_\kappa$')
ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
ax.legend((rects3[0], rects4[0]), ('45 23 gal+1/r+', '120 23 gal+1/r+'), bbox_to_anchor=(0.6, 0.5), fontsize=10)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
#plt.savefig('%skappahistbar_noshear.png' % root, dpi=250)
plt.savefig('%skappabiashistbar_noshear.png' % root, dpi=250)

ax = plt.subplot(2,1,1)

col5 = (kappastat_45[1])
rects5 = ax.bar(ind + width, col5, width, color='r')
col6 = (kappastat_120[1])
rects6 = ax.bar(ind + 2*width, col6, width, color='b')

#ax.set_ylim([0.0,0.2])
ax.set_ylim([-0.1,0.1])
#ax.set_ylabel('median$_\kappa$')
ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')

ax = plt.subplot(2,1,2)

col7 = (kappastat_45[3])
rects7 = ax.bar(ind + width, col7, width, color='r')
col8 = (kappastat_120[3])
rects8 = ax.bar(ind + 2*width, col8, width, color='b')

ax.set_ylim([0,0.1])
#ax.set_ylabel('$\sigma_\kappa$')
ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
ax.legend((rects7[0], rects8[0]), ('45 23 gal+1/r+$\gamma$+', '120 23 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.6, 0.3), fontsize=10)
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
#plt.savefig('%skappahistbar_shear.png' % root, dpi=250)
plt.savefig('%skappabiashistbar_shear.png' % root, dpi=250)

plt.clf()
