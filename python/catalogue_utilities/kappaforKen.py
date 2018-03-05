# Produces a normalized distribution with the kappa axis attached. the third axis is the smooth distribution

import scipy as sp
from scipy.stats import norm
import numpy as np

min_kappa = -0.10
max_kappa = 1
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

kappa  = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_chameleon_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_23_meds_increments2_2_2_2_2.cat" % root, usecols=[0], unpack=True)
median,stddev,kappa_values = statistics(kappa,bin_stat,min_kappa,max_kappa)
kappa = kappa / np.sum(kappa * np.abs((kappa_values[:-1]+halfwidth)))

winlen = 12
#smooth(kappa_3,winlen,'flat')
#smooth(kappa_3,winlen,'hanning')
#smooth(kappa_3,winlen,'hamming')
#smooth(kappa_3,winlen,'bartlett')
#smooth(kappa_3,winlen,'blackman')
data = np.c_[kappa_values[:-1],kappa,smooth(kappa,winlen,'flat')[(winlen/2-1):-(winlen/2)]]
np.savetxt('%skappaforKen_WFI2033_5innermask_nobeta_zgap-1.0_-1.0_chameleon_120_gal_120_gamma_120_oneoverr_45_gal_45_oneoverr_23_meds_increments2_2_2_2_2.cat' % root, data, fmt='%.4f %.3e %.3e')
