
# run from the Guo_galaxies folder as: python /Users/perseus/Dropbox/Davis_work/code/kappaplotcustom.py


import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
from scipy import stats
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
from numpy import linspace

file1="GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_orig_size45_i24_ratioquick_kappa.dat"
file2="GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size45_i24_ratioquick_gal_2.0_0.2_kappa.dat"
file3="GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size45_i24_ratioquick_gal_2.0_0.2_oneoverr_2.2_0.2_kappa.dat"
file4="GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size45_i24_ratioquick_gal_2.0_0.2_oneoverr_2.2_0.2_massoverr_4.0_1.0_kappa.dat"

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
data3 = np.loadtxt(file3)
data4 = np.loadtxt(file4)
kappa1=data1.T
kappa2=data2.T
kappa3=data3.T
kappa4=data4.T



output=str(file1)[0:len(str(file1))-21]+"_custom_gal_oneoverr_mass.eps"
BINS=100
plt.suptitle(r'%s' % output[68:len(str(file1))-21], fontsize=15, y=0.998)
#x = linspace(0,2,500)
#plt.subplot(451)
#n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.subplot(1,1,1)
if file1[101:len(file1)-10]!="":
    plt.hist(kappa1, histtype='step', color='b', label="No constraints", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
else:
    plt.hist(kappa1, histtype='step', color='b', label="No constraints", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
if file3[101:len(file3)-10]!="":
    plt.hist(kappa3, histtype='step', color='g', label="ratios: count 2.0$\pm$0.2 $1/r$ 2.2$\pm$0.2", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
else:
    plt.hist(kappa3, histtype='step', color='g', label="ratios: count 2.0$\pm$0.2 $1/r$ 2.2$\pm$0.2", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
if file2[101:len(file2)-10]!="":
    plt.hist(kappa2, histtype='step', color='r', label="ratios: count 2.0$\pm$0.2", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
else:
    plt.hist(kappa2, histtype='step', color='r', label="ratios: count 2.0$\pm$0.2", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
if file4[101:len(file4)-10]!="":
    plt.hist(kappa4, histtype='step', color='m', label="ratios: count 2.0$\pm$0.2 $1/r$ 2.2$\pm$0.2 $M/r$ 4.1$\pm$1.0", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
else:
    plt.hist(kappa4, histtype='step', color='m', label="ratios: count 2.0$\pm$0.2 $1/r$ 2.2$\pm$0.2 $M/r$ 4.1$\pm$1.0", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])

#if file1[101:len(file1)-10]!="":
#    plt.hist(kappa1, histtype='step', color='b', label=file1[101:len(file1)-10], linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#else:
#    plt.hist(kappa1, histtype='step', color='b', label="no", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#if file3[101:len(file3)-10]!="":
#    plt.hist(kappa3, histtype='step', color='g', label=file3[101:len(file3)-10], linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#else:
#    plt.hist(kappa3, histtype='step', color='g', label="no", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#if file4[101:len(file4)-10]!="":
#    plt.hist(kappa4, histtype='step', color='m', label=file4[101:len(file4)-10], linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#else:
#    plt.hist(kappa4, histtype='step', color='m', label="no", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#if file2[101:len(file2)-10]!="":
#    plt.hist(kappa2, histtype='step', color='r', label=file2[101:len(file2)-10], linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#else:
#    plt.hist(kappa2, histtype='step', color='r', label="no", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
#plt.hist(kappa4, histtype='step', color='m', label="no", linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
ax=plt.subplot(111)
s = "$\kappa_\mathrm{med}$=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa1),517500)
ax.text(0.3, 0.75, s, fontsize=15, color='b',transform=ax.transAxes)
s = "$\kappa_\mathrm{med}$=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa2),len(kappa2))
ax.text(0.3, 0.70, s, fontsize=15, color='r',transform=ax.transAxes)
s = "$\kappa_\mathrm{med}$=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa3),len(kappa3))
ax.text(0.3, 0.65, s, fontsize=15, color='g',transform=ax.transAxes)
s = "$\kappa_\mathrm{med}$=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa4),len(kappa4))
ax.text(0.3, 0.60, s, fontsize=15, color='m',transform=ax.transAxes)
#s = "med=%.3f, std=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa1),np.std(kappa1),len(kappa1))
#ax.text(0.3, 0.75, s, fontsize=15, color='b',transform=ax.transAxes)
#s = "med=%.3f, std=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa2),np.std(kappa2),len(kappa2))
#ax.text(0.3, 0.70, s, fontsize=15, color='r',transform=ax.transAxes)
#s = "med=%.3f, std=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa3),np.std(kappa3),len(kappa3))
#ax.text(0.3, 0.65, s, fontsize=15, color='g',transform=ax.transAxes)
#s = "med=%.3f, std=%.3f, $N\_\mathrm{LOS}$=%d" % (np.average(kappa4),np.std(kappa4),len(kappa4))
#ax.text(0.3, 0.60, s, fontsize=15, color='m',transform=ax.transAxes)
#text(0.5, 0.5,'matplotlib',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
#plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=20)
plt.xlabel(r'$\kappa_\mathrm{ext}$', fontsize=20)
plt.tick_params(axis='x', labelsize=13)
plt.tick_params(axis='y', labelsize=13)
plt.setp(plt.xticks()[1], rotation=90)
plt.legend(bbox_to_anchor=(0.30, 0.9), loc='center left', borderaxespad=0., fontsize=10)
#plt.subplots_adjust(top=0.6)
#plt.tight_layout()
plt.savefig('%s' % output, dpi=500)
#plt.show()
print 'Done!'


