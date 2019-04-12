# produces density plot

from matplotlib.colors import LogNorm
import scipy.optimize as optimization
from pylab import *
import numpy as np

font = 10
ticksize = 10

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
x, y = np.loadtxt(file1, usecols=(0, 1), unpack=True)
xlim = 0.2
xlim_ = -0.2
ylim = 0.2
ylim_ = 0.0
x_ = x[(x >= xlim_) & (x <= xlim) & (y >= ylim_) & (y <= ylim)]
y_ = y[(x >= xlim_) & (x <= xlim) & (y >= ylim_) & (y <= ylim)]
x = x_
y = y_
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
colorbar()
#delta = (y-x)/(1+x)
plt.xlabel('$\kappa$', fontsize=font)
plt.ylabel('$\gamma$', fontsize=font)
plt.xlim(xlim_, xlim)
plt.ylim(ylim_, ylim)
plt.savefig('../WFI2033/MSkapparesults/1.png' , dpi=250)
