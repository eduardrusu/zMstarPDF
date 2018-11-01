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

kappastat= np.c_[#kappastat[:,1-1], # 1
                        #kappastat[:,2-1], # 2
                        kappastat[:,3-1], # 3
                        #kappastat[:,4-1], # 4
                        #kappastat[:,5-1], # 5
                        #kappastat[:,6-1], # 6
                        kappastat[:,7-1], # 7
                        kappastat[:,42-1], # 20
                        #kappastat[:,8-1], # 8
                        kappastat[:,117-1], # 23
                        kappastat[:,121-1], # 27
                        kappastat[:,11-1],
                        kappastat[:,47-1], # 21
                        kappastat[:,46-1],
                        #kappastat[:,100-1], # 22
                        #kappastat[:,99-1],

                        #kappastat[:,12-1], # 9
                        kappastat[:,14-1], # 10
                        #kappastat[:,32-1], # 11
                        kappastat[:,33-1], # 12
                        kappastat[:,34-1], # 13
                        #kappastat[:,35-1], # 14
                        kappastat[:,36-1], # 15
                        #kappastat[:,37-1], # 16
                        #kappastat[:,38-1], # 17
                        #kappastat[:,39-1], # 18
                        #kappastat[:,41-1], # 19
                        #kappastat[:,118-1], # 24
                        #kappastat[:,119-1], # 25
                        #kappastat[:,120-1], # 26
                        #kappastat[:,122-1], # 28
                        #kappastat[:,123-1], # 29
                        #kappastat[:,124-1], # 30
                        #kappastat[:,125-1], # 31
]

labels = (#'22.5cham120galgamma1/r', # 1
            #'22.5cham120galgamma', # 2
            '22.5cham45gal120gamma', # 3
            #'22.5cham45gal1/r120galgamma1/r', # 4
            #'22.5comp120galgamma1/r', # 5
            #'22.5comp120galgamma', # 6
            '22.5comp45gal120gamma', # 7
            '22.5fid45gal120galgamma', # 20
            #'22.5comp45gal1/r120galgamma1/r', # 8
            '22.5group120galgamma1/r', # 23
            '22.5twogroup120galgamma1/r', # 27
            '22.5fid120galgamma1/r', #
            '23fid45galgamma1/r', # 21
            '22.5fid45galgamma1/r',
            #'23fid45gal1/r', # 22
            #'22.5fid45gal1/r',

            #'23fid120galgamma1/r', # 9
            '22.5fid45gal120galgamma', # 10
            #'23fid120gal1/r', # 11
            '22.5fid45galgamma120gal1/r', # 12
            '22.5fid45galgammamass120gal1/r', # 13
            #'23fid45galgamma1/r120gal1/r', # 14
            '22.5fid45galgammaz120gal1/r', # 15
            #'22.5fid45galmass120gal1/r', # 16
            #'23fid45gal1/r120gal1/r', # 17
            #'22.5fid45galz120gal1/r', # 18
            #'22.5fid45gal120gal', # 19
            #'22.5group45galgamma1/r120gal1/r', # 24
            #'22.5group45gal1/r120gal1/r', # 25
            #'22.5group45galgamma1/r', # 26
            #'22.5twogroup120gal1/r', # 28
            #'22.5twogroup45galgamma1/r120gal1/r', # 29
            #'22.5twogroup45gal1/r120gal1/r', # 30
            #'22.5twogroup45galgamma1/r', # 31
)

N = 12
ind = 2 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

ax = plt.subplot(2,1,1)

col1 = (kappastat[0])
rects1 = ax.bar(ind + 2*width, col1, width, color='b')

#ax.set_ylim([0.00,0.05])
ax.set_ylim([-0.00,0.12])
ax.set_ylabel('median$_\kappa$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(labels, fontsize=6, rotation='vertical')

ax = plt.subplot(2,1,2)

col2 = (kappastat[1])
rects2 = ax.bar(ind + 2*width, col2, width, color='b')

ax.set_ylim([0,0.10])
ax.set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels([])
#ax.legend((rects1[0], rects2[0]), ('22.5 45: gal+1/r 120: gal+$\gamma$+', '22.5 45: gal+1/r 120: gal+1/r+$\gamma$+'), bbox_to_anchor=(0.65, 1.4), fontsize=10)
#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
plt.subplots_adjust(left=0.15, bottom=0.02, right=0.95, top=0.95, wspace=0.7, hspace=1.2)
plt.savefig('%skappashistbar-others.png' % root, dpi=250)

plt.clf()
