import numpy
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
    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

from numpy import *
from pylab import *

'''
def smooth_demo():
    
    t=linspace(-4,4,100)
    x=sin(t)
    xn=x+randn(len(t))*0.1
    y=smooth(x)
    ws=31
    subplot(211)
    plot(ones(ws))
    windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    hold(True)
    for w in windows[1:]:
        eval('plot('+w+'(ws) )')
    axis([0,30,0,1.1])
    legend(windows)
    title("The smoothing windows")
    subplot(212)
    plot(x)
    plot(xn)
    for w in windows:
        plot(smooth(xn,10,w))
    l=['original signal', 'signal with noise']
    l.extend(windows)
    legend(l)
    title("Smoothing a noisy signal")
    show()

if __name__=='__main__':
    smooth_demo()
'''

#plt.clf()
#x = np.linspace(-0.1,1,2000)
#y = np.loadtxt("kappahist_WFI2033_5innermask_nobeta_gal_gamma_oneoverr_23_45_meds_increments2_2_2.cat",unpack=True)
#winlen = 36
#z = smooth(y,winlen,'flat')
#z = smooth(y,winlen,'hanning')
#z = smooth(y,winlen,'hamming')
#z = smooth(y,winlen,'bartlett')
#z = smooth(y,winlen,'blackman')
#plt.plot(x,y)
#plt.plot(x,z[(winlen/2-1):-(winlen/2)])
#plt.show()
