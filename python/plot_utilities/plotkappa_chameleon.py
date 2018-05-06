# Plots the equivalent of Figure 13 in Rusu et al. 2017

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.stats import norm
import numpy as np

min_kappa = -0.10
max_kappa = 1
min_kappa_plot = -0.1
max_kappa_plot = 0.30
bin_stat = 2000
halfwidth = (max_kappa - min_kappa) / (bin_stat * 2.0)

#root = "/Users/cerusu/Dropbox/"
root = "/Volumes/LaCieSubaru/kapparesults/"

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

kappa_0  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_chameleon_120_gal_120_gamma_45_gal_23_meds_increments2_2_2.cat" % root, usecols=[0], unpack=True)
median0,stddev0,kappa_values = statistics(kappa_0,bin_stat,min_kappa,max_kappa)
kappa_0 = kappa_0 / np.sum(kappa_0 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_1  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_chameleon_120_gal_120_gamma_120_oneoverr_45_gal_23_meds_increments2_2_2_2.cat" % root, usecols=[0], unpack=True)
median1,stddev1,kappa_values = statistics(kappa_1,bin_stat,min_kappa,max_kappa)
kappa_1 = kappa_1 / np.sum(kappa_1 * np.abs((kappa_values[:-1]+halfwidth)))

kappa_2  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_chameleon_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_23_meds_increments2_2_2_2_2.cat" % root, usecols=[0], unpack=True)
median2,stddev2,kappa_values = statistics(kappa_2,bin_stat,min_kappa,max_kappa)
kappa_2 = kappa_2 / np.sum(kappa_2 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_3  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_handpicked_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_23_meds_increments2_16_8_2_2.cat" % root, usecols=[0], unpack=True)
#median3,stddev3,kappa_values = statistics(kappa_3,bin_stat,min_kappa,max_kappa)
#kappa_3 = kappa_3 / np.sum(kappa_3 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_4  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_23_meds_increments2_2_2_2_2.cat" % root, usecols=[0], unpack=True)
#median4,stddev4,kappa_values = statistics(kappa_4,bin_stat,min_kappa,max_kappa)
#kappa_4 = kappa_4 / np.sum(kappa_4 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_5  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_fiducial_45_gal_45_gamma_120_gal_23_meds_increments2_2_2.cat" % root, usecols=[0], unpack=True)
#median5,stddev5,kappa_values = statistics(kappa_5,bin_stat,min_kappa,max_kappa)
#kappa_5 = kappa_5 / np.sum(kappa_5 * np.abs((kappa_values[:-1]+halfwidth)))

#kappa_6  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_fiducial_45_gal_45_gamma_45_oneoverr_120_gal_23_meds_increments2_2_2_2.cat" % root, usecols=[0], unpack=True)
#median6,stddev6,kappa_values = statistics(kappa_6,bin_stat,min_kappa,max_kappa)
#kappa_6 = kappa_6 / np.sum(kappa_6 * np.abs((kappa_values[:-1]+halfwidth)))

#s = "med=%.3f std=%.3f LOS=%d" % (median,std1,LOS)
s0 = "med=%.3f std=%.3f" % (median0,stddev0)
s1 = "med=%.3f std=%.3f" % (median1,stddev1)
s2 = "med=%.3f std=%.3f" % (median2,stddev2)
#s3 = "med=%.3f std=%.3f" % (median3,stddev3)
#s4 = "med=%.3f std=%.3f" % (median4,stddev4)
#s5 = "med=%.3f std=%.3f" % (median5,stddev5)
#s6 = "med=%.3f std=%.3f" % (median6,stddev6)
plt.subplot(1,1,1)
ax = plt.subplot(1,1,1)
ax.tick_params(labelsize=15)
plt.xlim(min_kappa_plot, max_kappa_plot)
plt.ylim(0, 0.02)

plt.plot(kappa_values[:-1][::1],kappa_0[::1],linewidth=2, label ='$120: 1 + \gamma$; 45: 1', linestyle=':') # every 1th point
ax.text(0.6, 0.90, s0, fontsize=10, transform=ax.transAxes)
plt.plot(kappa_values[:-1][::1],kappa_1[::1], linewidth=2, label ='$120: 1 + 1/r + \gamma$; 45: 1', linestyle='-.') # every 1th point
ax.text(0.6, 0.85, s1, fontsize=10, transform=ax.transAxes)
plt.plot(kappa_values[:-1][::1],kappa_2[::1], linewidth=2, label ='$120: 1 + 1/r + \gamma$; 45: 1 + 1/r', linestyle='--') # every 1th point
ax.text(0.6, 0.80, s2, fontsize=10, transform=ax.transAxes)
#plt.plot(kappa_values[:-1][::1],kappa_3[::1], linewidth=2, label ='$120: 1 + 1/r + \gamma$; 45: 1 + 1/r') # every 1th point
#ax.text(0.6, 0.75, s3, fontsize=10,transform=ax.transAxes)
#plt.plot(kappa_values[:-1][::1],kappa_4[::1],linewidth=2, label ='$120: 1 + 1/r + \gamma$; 45: 1 + 1/r') # every 1th point
#ax.text(0.6, 0.70, s4, fontsize=10,transform=ax.transAxes)
#plt.plot(kappa_values[:-1][::1],kappa_5[::1],linewidth=2, label ='$45: 1 + \gamma$; 120: 1') # every 1th point
#ax.text(0.6, 0.65, s5, fontsize=10,transform=ax.transAxes)
#plt.plot(kappa_values[:-1][::1],kappa_6[::1],linewidth=2, label ='$45: 1 + 1/r + \gamma$; 120: 1') # every 1th point
#ax.text(0.6, 0.60, s6, fontsize=10,transform=ax.transAxes)

winlen = 12
#smooth(kappa_3,winlen,'flat')
#smooth(kappa_3,winlen,'hanning')
#smooth(kappa_3,winlen,'hamming')
#smooth(kappa_3,winlen,'bartlett')
#smooth(kappa_3,winlen,'blackman')
#plt.plot(kappa_values[:-1],smooth(kappa_3,winlen,'flat')[(winlen/2-1):-(winlen/2)],color='k', linewidth=2, label ='$1 + 1/r + \gamma$')
plt.xlabel(r'$\kappa$', fontsize=20)
plt.ylabel(r'normalized counts', fontsize=20)
plt.legend(loc="lower right")
plt.savefig('%skappahist_chameleon.png' % root, dpi=250, bbox_inches='tight')
