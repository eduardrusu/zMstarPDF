#

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np

min_kappa = -0.20
max_kappa = 1
#min_kappa_plot = -0.1
#max_kappa_plot = 0.15
bin_stat = 2000
halfwidth = (max_kappa - min_kappa) / (bin_stat * 2.0)

root = "/Users/cerusu/Dropbox/Davis_work/code/0408/test/"

def statistics(kappa_all_,bin_stat_,min_kappa_,max_kappa_):
    a, kappa_values = np.histogram([0], bins = bin_stat_, range=(min_kappa_,max_kappa_)) # create an empty histogram of the correct shape

    sum = np.sum(kappa_all_)
    #meanX = np.sum(kappa_counts * (kappa_values[:-1] + halfwidth)) / sum
    #meanX2 = np.sum(kappa_counts * (kappa_values[:-1] + halfwidth) ** 2) / sum
    #std = np.sqrt(meanX2 - meanX**2)

    med = 0
    i = 0
    ok = False
    while (med <= sum/2.0) and (ok == False):
        med = med + kappa_all_[i]
        if med > sum/2.0:
            median = kappa_values[i] + halfwidth
            ok = True
        if med == sum/2.0:
            median = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1

    std = 0
    ok = False
    i = 0
    while (std <= sum * 0.16) and (ok == False):
        std = std + kappa_all_[i]
        if std > sum * 0.16:
            std1_ = kappa_values[i] + halfwidth
            ok = True
        if med == sum*0.16:
            std1_ = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1

    std = 0
    ok = False
    i = 0
    while (std <= sum*0.84) and (ok == False):
        std = std + kappa_all_[i]
        if std > sum*0.84:
            std1 = kappa_values[i] + halfwidth
            ok = True
        if med == sum*0.84:
            std1 = kappa_values[i] + 2 * halfwidth
            ok = True
        i = i + 1

    stddev = (std1 - std1_) / 2

    return median,stddev,kappa_values

kappa_45_1_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_z1,stddev0,kappa_values = statistics(kappa_45_1_z1,bin_stat,min_kappa,max_kappa)
kappa_45_1_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_z2,stddev0,kappa_values = statistics(kappa_45_1_z2,bin_stat,min_kappa,max_kappa)
kappa_45_1_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_z05,stddev0,kappa_values = statistics(kappa_45_1_z05,bin_stat,min_kappa,max_kappa)

kappa_45_1_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_z1_H,stddev0,kappa_values = statistics(kappa_45_1_z1_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_z2_H,stddev0,kappa_values = statistics(kappa_45_1_z2_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_22.5_med_increments2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_z05_H,stddev0,kappa_values = statistics(kappa_45_1_z05_H,bin_stat,min_kappa,max_kappa)

kappa_45_1_1r_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z1,stddev0,kappa_values = statistics(kappa_45_1_1r_z1,bin_stat,min_kappa,max_kappa)
kappa_45_1_1r_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z2,stddev0,kappa_values = statistics(kappa_45_1_1r_z2,bin_stat,min_kappa,max_kappa)
kappa_45_1_1r_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z05,stddev0,kappa_values = statistics(kappa_45_1_1r_z05,bin_stat,min_kappa,max_kappa)

kappa_45_1_1r_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z1_H,stddev0,kappa_values = statistics(kappa_45_1_1r_z1_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_1r_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z2_H,stddev0,kappa_values = statistics(kappa_45_1_1r_z2_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_1r_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_oneoverr_22.5_med_increments2_2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_1r_z05_H,stddev0,kappa_values = statistics(kappa_45_1_1r_z05_H,bin_stat,min_kappa,max_kappa)

kappa_45_1_zr_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z1,stddev0,kappa_values = statistics(kappa_45_1_zr_z1,bin_stat,min_kappa,max_kappa)
kappa_45_1_zr_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z2,stddev0,kappa_values = statistics(kappa_45_1_zr_z2,bin_stat,min_kappa,max_kappa)
kappa_45_1_zr_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z05,stddev0,kappa_values = statistics(kappa_45_1_zr_z05,bin_stat,min_kappa,max_kappa)

kappa_45_1_zr_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z1_H,stddev0,kappa_values = statistics(kappa_45_1_zr_z1_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_zr_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z2_H,stddev0,kappa_values = statistics(kappa_45_1_zr_z2_H,bin_stat,min_kappa,max_kappa)
kappa_45_1_zr_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_45_gal_45_zoverr_22.5_med_increments2_2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_45_1_zr_z05_H,stddev0,kappa_values = statistics(kappa_45_1_zr_z05_H,bin_stat,min_kappa,max_kappa)

kappa_120_1_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_z1,stddev0,kappa_values = statistics(kappa_120_1_z1,bin_stat,min_kappa,max_kappa)
kappa_120_1_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_z2,stddev0,kappa_values = statistics(kappa_120_1_z2,bin_stat,min_kappa,max_kappa)
kappa_120_1_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_z05,stddev0,kappa_values = statistics(kappa_120_1_z05,bin_stat,min_kappa,max_kappa)

kappa_120_1_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_z1_H,stddev0,kappa_values = statistics(kappa_120_1_z1_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_z2_H,stddev0,kappa_values = statistics(kappa_120_1_z2_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_22.5_med_increments2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_z05_H,stddev0,kappa_values = statistics(kappa_120_1_z05_H,bin_stat,min_kappa,max_kappa)

kappa_120_1_1r_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z1,stddev0,kappa_values = statistics(kappa_120_1_1r_z1,bin_stat,min_kappa,max_kappa)
kappa_120_1_1r_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z2,stddev0,kappa_values = statistics(kappa_120_1_1r_z2,bin_stat,min_kappa,max_kappa)
kappa_120_1_1r_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z05,stddev0,kappa_values = statistics(kappa_120_1_1r_z05,bin_stat,min_kappa,max_kappa)

kappa_120_1_1r_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z1_H,stddev0,kappa_values = statistics(kappa_120_1_1r_z1_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_1r_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z2_H,stddev0,kappa_values = statistics(kappa_120_1_1r_z2_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_1r_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_oneoverr_22.5_med_increments2_2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_1r_z05_H,stddev0,kappa_values = statistics(kappa_120_1_1r_z05_H,bin_stat,min_kappa,max_kappa)

kappa_120_1_zr_z1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z1,stddev0,kappa_values = statistics(kappa_120_1_zr_z1,bin_stat,min_kappa,max_kappa)
kappa_120_1_zr_z2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z2,stddev0,kappa_values = statistics(kappa_120_1_zr_z2,bin_stat,min_kappa,max_kappa)
kappa_120_1_zr_z05  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z05,stddev0,kappa_values = statistics(kappa_120_1_zr_z05,bin_stat,min_kappa,max_kappa)

kappa_120_1_zr_z1_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2Henriques_zeta1.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z1_H,stddev0,kappa_values = statistics(kappa_120_1_zr_z1_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_zr_z2_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2Henriques_zeta2.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z2_H,stddev0,kappa_values = statistics(kappa_120_1_zr_z2_H,bin_stat,min_kappa,max_kappa)
kappa_120_1_zr_z05_H  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_22.5_med_increments2_2Henriques_zeta05.cat" % root, usecols=[0], unpack=True)
median_120_1_zr_z05_H,stddev0,kappa_values = statistics(kappa_120_1_zr_z05_H,bin_stat,min_kappa,max_kappa)

plt.subplot(1,1,1)
ax = plt.subplot(1,1,1)
ax.tick_params(labelsize=15)
plt.clf()
#winlen = 12

plt.plot([0.5,1,2],[median_45_1_z05,median_45_1_z1,median_45_1_z2],color='k', linewidth=1, linestyle='-', label ='SA $45: 1$')
plt.plot([0.5,1,2],[median_45_1_z05_H,median_45_1_z1_H,median_45_1_z2_H],color='k', linewidth=1, linestyle='--', label ='Henriques $45: 1$')
plt.plot([0.5,1,2],[median_45_1_1r_z05,median_45_1_1r_z1,median_45_1_1r_z2],color='b', linewidth=1, linestyle='-', label ='SA $45: 1+1/r$')
plt.plot([0.5,1,2],[median_45_1_1r_z05_H,median_45_1_1r_z1_H,median_45_1_1r_z2_H],color='b', linewidth=1, linestyle='--', label ='Henriques $45: 1+1/r$')
plt.plot([0.5,1,2],[median_45_1_zr_z05,median_45_1_zr_z1,median_45_1_zr_z2],color='r', linewidth=1, linestyle='-', label ='SA $45: 1+z/r$')
plt.plot([0.5,1,2],[median_45_1_zr_z05_H,median_45_1_zr_z1_H,median_45_1_zr_z2_H],color='r', linewidth=1, linestyle='--', label ='Henriques $45: 1+z/r$')

plt.plot([0.5,1,2],[median_120_1_z05,median_120_1_z1,median_120_1_z2],color='k', linewidth=3, linestyle='-', label ='SA $120: 1$')
plt.plot([0.5,1,2],[median_120_1_z05_H,median_120_1_z1_H,median_120_1_z2_H],color='k', linewidth=3, linestyle='--', label ='Henriques $120: 1$')
plt.plot([0.5,1,2],[median_120_1_1r_z05,median_120_1_1r_z1,median_120_1_1r_z2],color='b', linewidth=3, linestyle='-', label ='SA $120: 1+1/r$')
plt.plot([0.5,1,2],[median_120_1_1r_z05_H,median_120_1_1r_z1_H,median_120_1_1r_z2_H],color='b', linewidth=3, linestyle='--', label ='Henriques $120: 1+1/r$')
plt.plot([0.5,1,2],[median_120_1_zr_z05,median_120_1_zr_z1,median_120_1_zr_z2],color='r', linewidth=3, linestyle='-', label ='SA $120: 1+z/r$')
plt.plot([0.5,1,2],[median_120_1_zr_z05_H,median_120_1_zr_z1_H,median_120_1_zr_z2_H],color='r', linewidth=3, linestyle='--', label ='Henriques $120: 1+z/r$')

plt.xlabel(r'$\zeta$', fontsize=20)
plt.ylabel(r'median ($\kappa)$', fontsize=20)
plt.xticks([0.5,1,2])
plt.legend(loc="upper left",fontsize=7)
plt.title("SA - Henriques comparison $z_s=2.4$",fontsize=14)
plt.savefig('%sSAHenriquescomparison0408.png' % root, dpi=250, bbox_inches='tight')
