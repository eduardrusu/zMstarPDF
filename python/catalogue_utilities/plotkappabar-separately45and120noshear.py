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

root = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappa/"
data = np.genfromtxt('%smedstd.dat' % root,dtype=['S1000','f8','f8','f8','f8'])

kappastat = np.array([])
for i in range(np.shape(data)[0]):
    if i == 0:
        kappastat = np.array([data[i][1],data[i][2],data[i][3],data[i][4]])
    else:
        x = np.array([data[i][1],data[i][2],data[i][3],data[i][4]])
        kappastat = np.c_[kappastat,x]

kappastat_120 = np.c_[kappastat[:,22-1], # 1-1/r
                        kappastat[:,17-1], # z
                        kappastat[:,14-1], # M_\star
                        kappastat[:,6-1], # M^2_\star
                        kappastat[:,10-1], # M^3_\star$
                        kappastat[:,19-1], # 1/r
                        kappastat[:,18-1], # z/r
                        kappastat[:,15-1], # M_\star/r
                        kappastat[:,7-1], # M^2_\star/r
                        kappastat[:,11-1], # M^3_\star/r
                        kappastat[:,9-1], # M^2_{\star\mathrm{rms}}
                        kappastat[:,13-1], # M^3_{\star\mathrm{rms}}
                        kappastat[:,8-1], # M^2_\star/r_\mathrm{,rms}
                        kappastat[:,12-1], # M^3_\star/r_\mathrm{,rms}
                        kappastat[:,5-1], # flexion
                        kappastat[:,16-1], # tidal
                        kappastat[:,3-1], # SIS
                        kappastat[:,4-1]] # SIShalo
kappastat_45 = np.c_[kappastat[:,25-1], # 1-1/r
                        kappastat[:,47-1], # z
                        kappastat[:,44-1], # M_\star
                        kappastat[:,36-1], # M^2_\star
                        kappastat[:,40-1], # M^3_\star$
                        kappastat[:,33-1], # 1/r
                        kappastat[:,48-1], # z/r
                        kappastat[:,45-1], # M_\star/r
                        kappastat[:,37-1], # M^2_\star/r
                        kappastat[:,41-1], # M^3_\star/r
                        kappastat[:,39-1], # M^2_{\star\mathrm{rms}}
                        kappastat[:,43-1], # M^3_{\star\mathrm{rms}}
                        kappastat[:,38-1], # M^2_\star/r_\mathrm{,rms}
                        kappastat[:,42-1], # M^3_\star/r_\mathrm{,rms}
                        kappastat[:,35-1], # flexion
                        kappastat[:,46-1], # tidal
                        kappastat[:,34-1], # SIS
                        kappastat[:,34-1]] # SIShalo

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
ax.set_ylabel('median$_\kappa$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')

ax = plt.subplot(2,1,2)

col3 = (kappastat_45[1])
rects3 = ax.bar(ind + width, col3, width, color='r')
col4 = (kappastat_120[1])
rects4 = ax.bar(ind + 2*width, col4, width, color='b')

ax.set_ylim([0,0.05])
ax.set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+', '120 22.5 gal+1/r+'), bbox_to_anchor=(0.65, 1.4), fontsize=10)
#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
plt.savefig('%skappashistbar-separately45and120noshear.png' % root, dpi=250)

plt.clf()
