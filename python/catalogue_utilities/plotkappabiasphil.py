# Use as python plot plotkappabiasphil.py WFI2033 5 23 45

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np
import sys
import os
from os import system

lens = str(sys.argv[1])
inner = str(sys.argv[2])
mag = str(sys.argv[3])
rad = str(sys.argv[4])

bins = np.linspace(-0.20,0.20,300)
halfwidth = (bins[1]-bins[0])/2

root = "/Users/eduardrusu/Dropbox/Davis_work/code/GOODCODE/WFI2033kappa/"

col1  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col2  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_z_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col3  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col4  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass2_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col5  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass3_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col6  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col7  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_zoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col8  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_massoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col9  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass2overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col10  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass3overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col11  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass2rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col12  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass3rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col13  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass2overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col14  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_mass3overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col15  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_flexion_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col16  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_tidal_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col17  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_SIS_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col18  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_oneoverr_SIShalo_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
'''
col1  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col2  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_z_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col3  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col4  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col5  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col6  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col7  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_zoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col8  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_massoverr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col9  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col10  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overr_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col11  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col12  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3rms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col13  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass2overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col14  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_mass3overrrms_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col15  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_flexion_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col16  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_tidal_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col17  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIS_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
col18  = np.loadtxt("%skappasimbiasphil_%s_%sinnermask_nobeta_gal_gamma_oneoverr_SIShalo_%s_%s_meds_8of64fields_increments.cat" % (root,lens,inner,mag,rad), unpack=True)
'''
min_kappa_plot = -0.04
max_kappa_plot = 0.04
nRows = 6
nCols = 3
nPlots = 18
fontlabel = 18

plt.clf()
#plt.axis([min_kappa_plot, max_kappa_plot, 0, 800])
fig = plt.figure()
ax1 = fig.add_subplot(6,3,1)
ax = plt.subplot(6,3,1, sharex=ax1, sharey=ax1)
#ax.set_aspect(1, adjustable='datalim')

ax=plt.subplot(6,3,1, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col1/np.sum(col1  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$1$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,2, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col2/np.sum(col2  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$z$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,3, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col3/np.sum(col3  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M_\star$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,4, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col4/np.sum(col4  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^2_\star$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,5, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col5/np.sum(col5  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^3_\star$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,6, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col6/np.sum(col6  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$1/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,7, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col7/np.sum(col7  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$z/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,8, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col8/np.sum(col8  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M_\star/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,9, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col9/np.sum(col9  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^2_\star/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,10, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col10/np.sum(col10  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^3_\star/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,11, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col11/np.sum(col11  *  (bins[1]-bins[0])),'k-',label="xxx")
#ax.text(0.8, 0.85, "$M/r^3$", fontsize=fontlabel, color='k',transform=ax.transAxes)
ax.text(0.6, 0.6, "$M^2_\mathrm{\star,rms}$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,12, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col12/np.sum(col12  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^3_\mathrm{\star,rms}$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,13, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col13/np.sum(col13  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^2_\star/r_\mathrm{,rms}$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,14, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col14/np.sum(col14  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M^3_\star/r_\mathrm{,rms}$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,15, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.plot(bins[:-1],col15/np.sum(col15  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M_\star/r^3$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,16, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
ax.set_xticklabels([-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03])
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col16/np.sum(col16  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$M_\star/r^2$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,17, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
ax.set_xticklabels([-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03])
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.xticks(rotation='vertical')
plt.plot(bins[:-1],col17/np.sum(col17  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.6, 0.6, "$\sqrt{M_\star}/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

ax=plt.subplot(6,3,18, sharex=ax1, sharey=ax1)
ax.tick_params(labelsize=14)
ax.set_xticklabels([-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04])
plt.xlim(min_kappa_plot, max_kappa_plot)

plt.plot(bins[:-1],col18/np.sum(col18  *  (bins[1]-bins[0])),'k-',label="xxx")
ax.text(0.5, 0.6, "$\sqrt{M_h}/r$", fontsize=fontlabel, color='k',transform=ax.transAxes)

# hide the plots with no data in the grid
ax=plt.subplot(6,3,16, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.xticks(rotation='vertical')
ax=plt.subplot(6,3,18, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.xticks(rotation='vertical')

index = 1
for r in range(1, nRows +1):
    for c in range(1, nCols + 1):
        ax = plt.subplot(nRows, nCols, index, sharex=ax1, sharey=ax1)
        index += 1
        # Turn off y tick labels for all but the first column.
        if ((c != 1) and (index <= nPlots)):
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
        # Turn off x tick lables for all but the bottom plot in each column.
        if ((nPlots - index) >= nCols):
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
        if index == 18:
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

fig.text(0.5, 0.025, '$\kappa-\kappa_\mathrm{true}$', ha='center', va='center', size='22')
plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.15)

plt.savefig('%sbiasphil_%s_%s_%s_%s.png' % (root,lens,inner,mag,rad), dpi=250)
#plt.savefig('%sbiasphil_%s_%s_%s_%s_gamma.png' % (root,lens,inner,mag,rad), dpi=250)
print 'Done!'
