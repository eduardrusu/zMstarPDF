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

kappastat_gamma = np.c_[kappastat[:,81-1], # 1-1/r
                        kappastat[:,79-1], # z
                        kappastat[:,75-1], # mass
                        kappastat[:,67-1], # mass2
                        kappastat[:,71-1], # mass3
                        kappastat[:,77-1], # 1/r
                        kappastat[:,80-1], # z/r
                        kappastat[:,76-1], # massoverr
                        kappastat[:,68-1], # mass2overr
                        kappastat[:,72-1], # mass3overr
                        kappastat[:,70-1], # mass2rms
                        kappastat[:,74-1], # mass3rms
                        kappastat[:,69-1], # mass2overrrms
                        kappastat[:,73-1], # mass3overrrms
                        kappastat[:,66-1], # flexion
                        kappastat[:,78-1], # tidal
                        kappastat[:,64-1], # SIS
                        kappastat[:,65-1]] # SIShalo
kappastat_nogamma = np.c_[kappastat[:137,-1], # 1-1/r
                        kappastat[:,96-1], # z
                        kappastat[:,93-1], # mass
                        kappastat[:,85-1], # mass2
                        kappastat[:,89-1], # mass3
                        kappastat[:,98-1], # 1/r
                        kappastat[:,97-1], # z/r
                        kappastat[:,94-1], # massoverr
                        kappastat[:,86-1], # mass2overr
                        kappastat[:,90-1], # mass3overr
                        kappastat[:,88-1], # mass2rms
                        kappastat[:,92-1], # mass3rms
                        kappastat[:,87-1], # mass2overrrms
                        kappastat[:,91-1], # mass3overrrms
                        kappastat[:,84-1], # flexion
                        kappastat[:,95-1], # tidal
                        kappastat[:,82-1], # SIS
                        kappastat[:,83-1]] # SIShalo

N = 18
ind = 2.5 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

ax = plt.subplot(2,1,1)

col1 = (kappastat_gamma[0])
rects1 = ax.bar(ind + width, col1, width, color='r')
col2 = (kappastat_nogamma[0])
rects2 = ax.bar(ind + 2*width, col2, width, color='b')

#ax.set_ylim([0.00,0.05])
ax.set_ylim([0.00,0.09])
ax.set_ylabel('median$_\kappa$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')

ax = plt.subplot(2,1,2)

col3 = (kappastat_gamma[1])
rects3 = ax.bar(ind + width, col3, width, color='r')
col4 = (kappastat_nogamma[1])
rects4 = ax.bar(ind + 2*width, col4, width, color='b')

ax.set_ylim([0,0.07])
ax.set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
ax.legend((rects1[0], rects2[0]), ('22.5 45: gal+1/r 120: gal+$\gamma$+', '22.5 45: gal+1/r 120: gal+1/r+'), bbox_to_anchor=(0.65, 1.4), fontsize=10)
#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
plt.savefig('%skappashistbar-conjoint45120gammanogamma.png' % root, dpi=250)

plt.clf()
