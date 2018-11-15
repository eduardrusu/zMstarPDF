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

kappastat_45 = np.c_[   kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_45_gamma_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_z_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass2_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass3_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_oneoverr_120_gal_120_gamma_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_zoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_massoverr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass2overr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass3overr_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass2rms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass3rms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass2overrrms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_mass3overrrms_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_flexion_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_tidal_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_SIS_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_45_SIShalo_120_gal_120_gamma_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

kappastat_120 = np.c_[  kappastat.T[kappastat[0]=='fiducial_120_gal_120_zoverr_45_gal_45_gamma_22.5_med_increments4_4_4_4_emptymsk.cat'][0][1:].astype(float), # 1-1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_z_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass2_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass3_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_oneoverr_22.5_med_increments4_4_4_emptymsk.cat'][0][1:].astype(float), # 1/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_zoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # z/r
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_massoverr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # massoverr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass2overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass3overr_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overr
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass2rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass3rms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3rms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass2overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass2overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_mass3overrrms_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # mass3overrrms
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_flexion_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # flexion
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_tidal_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # tidal
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_SIS_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float), # SIS
                        kappastat.T[kappastat[0]=='fiducial_45_gal_45_zoverr_120_gal_120_gamma_120_SIShalo_22.5_med_increments4_4_4_4_4_emptymsk.cat'][0][1:].astype(float)] # SIShalo

N = 18
ind = 2.5 * np.arange(N)  # the x locations for the groups
width = 0.8       # the width of the bars

ax = plt.subplot(2,1,1)

col1 = (kappastat_45[0])
rects1 = ax.bar(ind + width, col1, width, color='r')
col2 = (kappastat_120[0])
rects2 = ax.bar(ind + 2*width, col2, width, color='b')

#ax.set_ylim([0.00,0.05])
ax.set_ylim([-0.02,0.08])
ax.set_ylabel('median$_\kappa$')
#ax.set_ylabel('$\mathrm{median}_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')

ax = plt.subplot(2,1,2)

col3 = (kappastat_45[1])
rects3 = ax.bar(ind + width, col3, width, color='r')
col4 = (kappastat_120[1])
rects4 = ax.bar(ind + 2*width, col4, width, color='b')

ax.set_ylim([0,0.08])
ax.set_ylabel('$\sigma_\kappa$')
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(('$1-1/r$', '$z$',  '$M_\star$', '$M^2_\star$', '$M^3_\star$', '$1/r$', '$z/r$', '$M_\star/r$', '$M^2_\star/r$', '$M^3_\star/r$', '$M^2_{\star\mathrm{rms}}$', '$M^3_{\star\mathrm{rms}}$', '$M^2_\star/r_\mathrm{,rms}$', '$M^3_\star/r_\mathrm{,rms}$', '$M_\star/r^3$', '$M_\star/r^2$', '$\sqrt{M_\star}/r$', '$\sqrt{M_h}/r$'), fontsize=10, rotation='vertical')
ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.65, 1.4), fontsize=10)
#ax.legend((rects1[0], rects2[0]), ('45 22.5 gal+1/r+$\gamma$+', '120 22.5 gal+1/r+$\gamma$+'), bbox_to_anchor=(0.3, 0.97), fontsize=10)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.7, hspace=0.7)
plt.savefig('%skappashistbar-conjointgalzoverrgalgamma.png' % root, dpi=250)

plt.clf()
