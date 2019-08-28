# Plots the equivalent of Figure 13 in Rusu et al. 2017

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np

min_kappa = -0.20
max_kappa = 1
min_kappa_plot = -0.1
max_kappa_plot = 0.15
bin_stat = 2000
halfwidth = (max_kappa - min_kappa) / (bin_stat * 2.0)

#root = "/Users/cerusu/Dropbox/"
#root = "/Volumes/LaCieSubaru/kapparesults/"
root = "/Users/cerusu/Dropbox/Davis_work/code/0408/"

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

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.

        input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.

        output: the smoothed signal

        see also:

        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter

        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

plt.clf()

kappa_0  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_45_gal_45_gamma_45_oneoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median0,stddev0,kappa_values = statistics(kappa_0,bin_stat,min_kappa,max_kappa)
kappa_0 = kappa_0 / np.sum(kappa_0 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_1  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_45_gal_45_gamma_45_z_120_gal_120_z_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
#median1,stddev1,kappa_values = statistics(kappa_1,bin_stat,min_kappa,max_kappa)
#kappa_1 = kappa_1 / np.sum(kappa_1 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_2  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_45_gal_45_gamma_45_zoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median2,stddev2,kappa_values = statistics(kappa_2,bin_stat,min_kappa,max_kappa)
kappa_2 = kappa_2 / np.sum(kappa_2 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_3  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_120_gal_120_gamma_120_oneoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median3,stddev3,kappa_values = statistics(kappa_3,bin_stat,min_kappa,max_kappa)
kappa_3 = kappa_3 / np.sum(kappa_3 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_4  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median4,stddev4,kappa_values = statistics(kappa_4,bin_stat,min_kappa,max_kappa)
kappa_4 = kappa_4 / np.sum(kappa_4 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_5  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_120_gal_120_gamma_120_oneoverr_45_gal_45_zoverr_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median5,stddev5,kappa_values = statistics(kappa_5,bin_stat,min_kappa,max_kappa)
kappa_5 = kappa_5 / np.sum(kappa_5 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_6  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_composite_120_gal_120_gamma_120_zoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median6,stddev6,kappa_values = statistics(kappa_6,bin_stat,min_kappa,max_kappa)
kappa_6 = kappa_6 / np.sum(kappa_6 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_7  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_gamma_45_oneoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median7,stddev7,kappa_values = statistics(kappa_7,bin_stat,min_kappa,max_kappa)
kappa_7 = kappa_7 / np.sum(kappa_7 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_8  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_gamma_45_z_120_gal_120_z_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
#median8,stddev8,kappa_values = statistics(kappa_8,bin_stat,min_kappa,max_kappa)
#kappa_8 = kappa_8 / np.sum(kappa_8 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_9  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_gamma_45_zoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median9,stddev9,kappa_values = statistics(kappa_9,bin_stat,min_kappa,max_kappa)
kappa_9 = kappa_9 / np.sum(kappa_9 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_10  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_gamma_45_zoverr_120_gal_120_zoverr_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median10,stddev10,kappa_values = statistics(kappa_10,bin_stat,min_kappa,max_kappa)
kappa_10 = kappa_10 / np.sum(kappa_10 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_11  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_oneoverr_22.5_med_increments1_1Henriques.cat" % root, usecols=[0], unpack=True)
median11,stddev11,kappa_values = statistics(kappa_11,bin_stat,min_kappa,max_kappa)
kappa_11 = kappa_11 / np.sum(kappa_11 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_12  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_45_gal_45_zoverr_22.5_med_increments1_1Henriques.cat" % root, usecols=[0], unpack=True)
median12,stddev12,kappa_values = statistics(kappa_12,bin_stat,min_kappa,max_kappa)
kappa_12 = kappa_12 / np.sum(kappa_12 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_13  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_gamma_120_oneoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median13,stddev13,kappa_values = statistics(kappa_13,bin_stat,min_kappa,max_kappa)
kappa_13 = kappa_13 / np.sum(kappa_13 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_14  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median14,stddev14,kappa_values = statistics(kappa_14,bin_stat,min_kappa,max_kappa)
kappa_14 = kappa_14 / np.sum(kappa_14 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_15  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_gamma_120_oneoverr_45_gal_45_zoverr_22.5_med_increments1_1_1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median15,stddev15,kappa_values = statistics(kappa_15,bin_stat,min_kappa,max_kappa)
kappa_15 = kappa_15 / np.sum(kappa_15 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_16  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_gamma_120_zoverr_22.5_med_increments1_1_1_shearwithoutpriorHenriques.cat" % root, usecols=[0], unpack=True)
median16,stddev16,kappa_values = statistics(kappa_16,bin_stat,min_kappa,max_kappa)
kappa_16 = kappa_16 / np.sum(kappa_16 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_17  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_oneoverr_22.5_med_increments1_1Henriques.cat" % root, usecols=[0], unpack=True)
median17,stddev17,kappa_values = statistics(kappa_17,bin_stat,min_kappa,max_kappa)
kappa_17 = kappa_17 / np.sum(kappa_17 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_18  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_oneoverr_45_gal_45_oneoverr_22.5_med_increments1_1_1_1Henriques.cat" % root, usecols=[0], unpack=True)
median18,stddev18,kappa_values = statistics(kappa_18,bin_stat,min_kappa,max_kappa)
kappa_18 = kappa_18 / np.sum(kappa_18 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_19  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_oneoverr_45_gal_45_zoverr_22.5_med_increments1_1_1_1Henriques.cat" % root, usecols=[0], unpack=True)
median19,stddev19,kappa_values = statistics(kappa_19,bin_stat,min_kappa,max_kappa)
kappa_19 = kappa_19 / np.sum(kappa_19 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_20  = np.loadtxt("%skappahist_0408_measured_5innermask_nobeta_removehandpicked_zgap-1.0_-1.0_powerlaw_120_gal_120_zoverr_22.5_med_increments1_1Henriques.cat" % root, usecols=[0], unpack=True)
median20,stddev20,kappa_values = statistics(kappa_20,bin_stat,min_kappa,max_kappa)
kappa_20 = kappa_20 / np.sum(kappa_20 * np.abs((kappa_values[:-1]+halfwidth)))

#s = "%.3f %.3f LOS=%d" % (median,std1,LOS)
s0 = "%.3f %.3f" % (median0,stddev0)
#s1 = "%.3f %.3f" % (median1,stddev1)
s2 = "%.3f %.3f" % (median2,stddev2)
s3 = "%.3f %.3f" % (median3,stddev3)
s4 = "%.3f %.3f" % (median4,stddev4)
s5 = "%.3f %.3f" % (median5,stddev5)
s6 = "%.3f %.3f" % (median6,stddev6)
s7 = "%.3f %.3f" % (median7,stddev7)
#s8 = "%.3f %.3f" % (median8,stddev8)
s9 = "%.3f %.3f" % (median9,stddev9)
s10 = "%.3f %.3f" % (median10,stddev10)
s11 = "%.3f %.3f" % (median11,stddev11)
s12 = "%.3f %.3f" % (median12,stddev12)
s13 = "%.3f %.3f" % (median13,stddev13)
s14 = "%.3f %.3f" % (median14,stddev14)
s15 = "%.3f %.3f" % (median15,stddev15)
s16 = "%.3f %.3f" % (median16,stddev16)
s17 = "%.3f %.3f" % (median17,stddev17)
s18 = "%.3f %.3f" % (median18,stddev18)
s19 = "%.3f %.3f" % (median19,stddev19)
s20 = "%.3f %.3f" % (median20,stddev20)

plt.subplot(1,1,1)
ax = plt.subplot(1,1,1)
ax.tick_params(labelsize=15)
plt.xlim(min_kappa_plot, max_kappa_plot)
plt.ylim(0, 0.25)

winlen = 12
#smooth(kappa_3,winlen,'flat')
#smooth(kappa_3,winlen,'hanning')
#smooth(kappa_3,winlen,'hamming')
#smooth(kappa_3,winlen,'bartlett')
#smooth(kappa_3,winlen,'blackman')

plt.plot(kappa_values[:-1],smooth(kappa_11,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='b', linewidth=1, linestyle=':', label ='%s; $45: 1,1/r$' %s11)
plt.plot(kappa_values[:-1],smooth(kappa_12,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='g', linewidth=1, linestyle=':', label ='%s; $45: 1,z/r$' %s12)
plt.plot(kappa_values[:-1],smooth(kappa_17,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='r', linewidth=1, linestyle=':', label ='%s; $120: 1,1/r$' %s17)
plt.plot(kappa_values[:-1],smooth(kappa_20,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='c', linewidth=1, linestyle=':', label ='%s; $120: 1,z/r$' %s20)
plt.plot(kappa_values[:-1],smooth(kappa_18,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='m', linewidth=1, linestyle=':', label ='%s; $120: 1,1/r; 45:1,1/r$' %s18)
plt.plot(kappa_values[:-1],smooth(kappa_19,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='k', linewidth=1, linestyle=':', label ='%s; $120: 1,1/r; 45:1,z/r$' %s19)

plt.plot(kappa_values[:-1],smooth(kappa_0,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='b', linewidth=1, linestyle='--', label ='%s; $45: 1,\gamma_c,1/r$' %s0)
plt.plot(kappa_values[:-1],smooth(kappa_2,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='g', linewidth=1, linestyle='--', label ='%s; $45: 1,\gamma_c,z/r$' %s2)
plt.plot(kappa_values[:-1],smooth(kappa_3,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='r', linewidth=1, linestyle='--', label ='%s; $120: 1,\gamma_c,1/r$' %s3)
plt.plot(kappa_values[:-1],smooth(kappa_6,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='c', linewidth=1, linestyle='--', label ='%s; $120: 1,\gamma_c,z/r$' %s6)
#plt.plot(kappa_values[:-1],smooth(kappa_1,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='m', linewidth=1, linestyle='--', label ='%s; $45: 1,\gamma_c,z; 120:1,z$' %s1)
plt.plot(kappa_values[:-1],smooth(kappa_4,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='m', linewidth=1, linestyle='--', label ='%s; $120: 1,\gamma_c,1/r; 45:1,1/r$' %s4)
plt.plot(kappa_values[:-1],smooth(kappa_5,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='k', linewidth=1, linestyle='--', label ='%s; $120: 1,\gamma_c,1/r; 45:1,z/r$' %s5)

plt.plot(kappa_values[:-1],smooth(kappa_7,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='b', linewidth=1, linestyle='-', label ='%s; $45: 1,\gamma_p,1/r$' %s7)
plt.plot(kappa_values[:-1],smooth(kappa_9,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='g', linewidth=1, linestyle='-', label ='%s; $45: 1,\gamma_p,z/r$' %s9)
plt.plot(kappa_values[:-1],smooth(kappa_13,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='r', linewidth=1, linestyle='-', label ='%s; $120: 1,\gamma_p,1/r$' %s13)
plt.plot(kappa_values[:-1],smooth(kappa_16,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='c', linewidth=1, linestyle='-', label ='%s; $120: 1,\gamma_p,z/r$' %s16)
#plt.plot(kappa_values[:-1],smooth(kappa_8,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='m', linewidth=1, linestyle='-', label ='%s; $45: 1,\gamma_p,z; 120:1,z$' %s8)
plt.plot(kappa_values[:-1],smooth(kappa_10,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='y', linewidth=1, linestyle='-', label ='%s; $45: 1,\gamma_p,z/r; 120:1,z/r$' %s10)
plt.plot(kappa_values[:-1],smooth(kappa_14,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='m', linewidth=1, linestyle='-', label ='%s; $120: 1,\gamma_p,1/r; 45:1,1/r$' %s14)
plt.plot(kappa_values[:-1],smooth(kappa_15,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='k', linewidth=1, linestyle='-', label ='%s; $120: 1,\gamma_p,1/r; 45:1,z/r$' %s15)

plt.xlabel(r'$\kappa$', fontsize=20)
plt.ylabel(r'normalized counts', fontsize=20)
plt.legend(loc="lower right",fontsize=7)
plt.title("Convergence distributions using Henriques galaxies",fontsize=14)
plt.savefig('%skappahist_handpickedHenriques.png' % root, dpi=250, bbox_inches='tight')
