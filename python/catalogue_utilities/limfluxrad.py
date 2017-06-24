# The code computes the histogram of Sextractor fluxradius parameter and fits it with a Gaussian. THis is useful for the star-galaxy classification employed by CFHTLenS

import numpy as np
import pylab as plt
from scipy.optimize import leastsq
plt.clf()
mag,rad=np.loadtxt("i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazy.cat",usecols=(4,6),unpack=True)
rad = rad[(mag<21) & (rad<6)]
fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
init = [20.0, 2.5, 1]
bins = 50
ydata,xdatabin = np.histogram(rad,bins)
bin = xdatabin[1]-xdatabin[0]
xdata = np.linspace(xdatabin[0]+bin/2.0,xdatabin[len(xdatabin)-1]-bin/2.0,len(xdatabin)-1)
out = leastsq( errfunc, init, args=(xdata, ydata))
c = out[0]
c[2] = np.abs(c[2]) # because by definition, c[2] may be negative
plt.plot(xdata, fitfunc(c, xdata))
plt.hist(rad, bins)
plt.title(r'$A = %.3f\  \mu = %.3f\  \sigma = %.3f\ $' %(c[0],c[1],c[2]));
plt.show()
