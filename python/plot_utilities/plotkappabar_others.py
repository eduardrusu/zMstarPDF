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
        kappastat = np.array([data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]])
    else:
        x = np.array([data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]])
        kappastat = np.c_[kappastat,x]

kappastat = np.c_[      kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_45_gal_45_gamma_45_oneoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_45_gal_45_z_45_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_gamma_45_gal_45_zoverr_45_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_45_gal_45_gamma_45_mass_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_45_gal_45_gamma_45_z_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_45_gal_45_mass_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_oneoverr_45_gal_45_z_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_zoverr_45_gal_45_gamma_45_mass_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_zoverr_45_gal_45_gamma_45_z_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_zoverr_45_gal_45_mass_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_120_gal_120_zoverr_45_gal_45_z_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_120_gal_120_gamma_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_120_gal_22.5_med_increments4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_120_gal_120_z_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_gamma_120_gal_120_zoverr_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_22.5_med_increments4_4_emptymsk.cat'][0][1:].astype(float)]#,  SIS]

labels = ( '120galgamma1/r45gal1/r', # 1
            '120galgamma45galgamma1/r', # 2
            '120galgamma45galzSIS', # 3
            '120galgamma45galz/rSIS', # 4
            '120gal1/r45galgammamass', # 5
            '120gal1/r45galgammaz', # 6
            '120gal1/r45galmass', # 7
            '120gal1/r45galz', # 20
            '120galz/r45galgammamass', # 8
            '120galz/r45galgammaz', # 23
            '120galz/r45galmass', # 27
            '120galz/r45galz', #
            '45gal120galgamma', # 21
            '45gal120gal',
            '45galgamma120galzSIS', # 22
            '45galgamma120galz/rSIS',
            '45galz/r'
)

N = 17
ind = 2 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

ax = plt.subplot(2,1,1)

col1 = (kappastat[0])
rects1 = ax.bar(ind + 2*width, col1, width, color='b')

#ax.set_ylim([0.00,0.05])
ax.set_ylim([-0.02,0.08])
ax.set_ylabel('median$_\kappa$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(labels, fontsize=6, rotation='vertical')

ax = plt.subplot(2,1,2)

col2 = (kappastat[1])
rects2 = ax.bar(ind + 2*width, col2, width, color='b')

ax.set_ylim([0,0.08])
ax.set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels([])
#ax.legend((rects1[0], rects2[0]), ('22.5 45: gal+1/r 120: gal+$\gamma$+', '22.5 45: gal+1/r 120: gal+1/r+$\gamma$+'), bbox_to_anchor=(0.65, 1.4), fontsize=10)
#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
plt.subplots_adjust(left=0.15, bottom=0.02, right=0.95, top=0.95, wspace=0.7, hspace=1.2)
plt.savefig('%skappashistbar-others.png' % root, dpi=250)

plt.clf()
