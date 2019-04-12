# produces density plot

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.optimize as optimization
from pylab import *
import numpy as np
import corner

font = 10
ticksize = 10

import sys
file = str(sys.argv[1])
plt.clf()
fig = plt.figure(figsize=(10,12))
#fig, axes = plt.subplots(nrows=2, ncols=2)

file1 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_22.5_med_increments4_1000_4_emptymsk_shearwithoutprior_test_kappagamma.cat"
file2 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_22.5_med_increments4_1000_4_emptymsk_shearwithoutprior_test_kappagamma_fiducialzeta.cat"
file3 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments4_1000_4_4_4_emptymsk_shearwithoutprior_test_kappagamma.cat"
file4 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments4_1000_4_4_4_emptymsk_shearwithoutprior_test_kappagamma_fiducialzeta.cat"
file5 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_gamma_45_oneoverr_22.5_med_increments4_1000_4_emptymsk_shearwithoutprior_test_kappagamma.cat"
file6 = "../WFI2033/MSkapparesults/kappahist_WFI2033_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_gamma_45_oneoverr_22.5_med_increments4_1000_4_emptymsk_shearwithoutprior_test_kappagamma_fiducialzeta.cat"
ax1 = fig.add_subplot(1,1,1)
ax1.set_aspect(1)
if file == "1": x, y = np.loadtxt(file1, usecols=(0, 1), unpack=True)
if file == "2": x, y = np.loadtxt(file2, usecols=(0, 1), unpack=True)
if file == "3": x, y = np.loadtxt(file3, usecols=(0, 1), unpack=True)
if file == "4": x, y = np.loadtxt(file4, usecols=(0, 1), unpack=True)
if file == "5": x, y = np.loadtxt(file5, usecols=(0, 1), unpack=True)
if file == "6": x, y = np.loadtxt(file6, usecols=(0, 1), unpack=True)
xlim = 0.15
xlim_ = 0.04
ylim = 0.15
ylim_ = 0.08
x_ = x[(x >= xlim_) & (x <= xlim) & (y >= ylim_) & (y <= ylim)]
y_ = y[(x >= xlim_) & (x <= xlim) & (y >= ylim_) & (y <= ylim)]
x = x_
y = y_
#hist2d(x, y, bins=[100, 100], norm=LogNorm())
#plt.xticks(rotation='vertical',size = ticksize)
#plt.yticks(size = ticksize)
#colorbar()

figure = corner.corner(np.c_[x,y], labels=np.linspace(1,5,5).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84, 0.955, 0.997],show_titles=True, title_kwargs={"fontsize": 12})

#plt.xlabel('$\kappa$', fontsize=font)
#plt.ylabel('$\gamma$', fontsize=font)
#plt.xlim(xlim_, xlim)
#plt.ylim(ylim_, ylim)
if file == "1": figure.savefig('../WFI2033/MSkapparesults/1.png' , dpi=250)
if file == "2": figure.savefig('../WFI2033/MSkapparesults/2.png' , dpi=250)
if file == "3": figure.savefig('../WFI2033/MSkapparesults/3.png' , dpi=250)
if file == "4": figure.savefig('../WFI2033/MSkapparesults/4.png' , dpi=250)
if file == "5": figure.savefig('../WFI2033/MSkapparesults/5.png' , dpi=250)
if file == "6": figure.savefig('../WFI2033/MSkapparesults/6.png' , dpi=250)
