##########################
# Simple code for scatter plot without error bars. Computes bias, scatter and fraction of outliers
##########################

#from matplotlib.colors import LogNorm
#import scipy.optimize as optimization
from pylab import *
import numpy as np
plt.clf()

#A_125_1 = 1.476526e+06
#B_125_1 = 2.233584e+05
#C_125_1 = 1.344782e+05
#D_125_1 = 7.739332e+04
#A_125_2 = 1.436465e+06
#B_125_2 = 2.133077e+05
#C_125_2 = 1.474414e+05
#D_125_2 = 6.260539e+04
#A_160_1 = 3.499114e+06
#B_160_1 = 3.805276e+05
#C_160_1 = 2.701676e+05
#D_160_1 = 1.401271e+05
#A_160_2 = 3.157662e+06
#B_160_2 = 3.534201e+05
#C_160_2 = 2.970849e+05
#D_160_2 = 1.389452e+05

A_125_1 = 10 ** ((-18.18 + 26.23)/2.5)
B_125_1 = 10 ** ((-20.28 + 26.23)/2.5)
C_125_1 = 10 ** ((-20.86 + 26.23)/2.5)
D_125_1 = 10 ** ((-21.45 + 26.23)/2.5)
A_125_2 = 10 ** ((-18.26 + 26.23)/2.5)
B_125_2 = 10 ** ((-20.31 + 26.23)/2.5)
C_125_2 = 10 ** ((-20.69 + 26.23)/2.5)
D_125_2 = 10 ** ((-21.74 + 26.23)/2.5)
A_160_1 = 10 ** ((-17.53 + 25.95)/2.5)
B_160_1 = 10 ** ((-20.05 + 25.95)/2.5)
C_160_1 = 10 ** ((-20.44 + 25.95)/2.5)
D_160_1 = 10 ** ((-21.09 + 25.95)/2.5)
A_160_2 = 10 ** ((-17.75 + 25.95)/2.5)
B_160_2 = 10 ** ((-20.09 + 25.95)/2.5)
C_160_2 = 10 ** ((-20.26 + 25.95)/2.5)
D_160_2 = 10 ** ((-21.11 + 25.95)/2.5)

muAnopert,muBnopert,muCnopert,muDnopert = np.loadtxt("pointSIEgamma_einstmagniftime_out_mcmc.dat",usecols=[14,10,6,2],unpack=True)
muApert,muBpert,muCpert,muDpert = np.loadtxt("point1pertSIEgamma_einstmagniftime_out_mcmc.dat",usecols=[14,10,6,2],unpack=True)

plt.hist(np.abs(muCnopert/muAnopert),bins=100,normed=True,label='predicted C/A w/o GX',color='r',histtype='stepfilled')
plt.hist(np.abs(muCnopert/muBnopert),bins=100,normed=True,label='predicted C/B w/o GX',color='b',histtype='stepfilled')
plt.hist(np.abs(muCnopert/muDnopert),bins=100,normed=True,label='predicted C/D w/o GX',color='g',histtype='stepfilled')
plt.hist(np.abs(muCpert/muApert),bins=100,normed=True,label='predicted C/A w/ GX',color='r',histtype='step')
plt.hist(np.abs(muCpert/muBpert),bins=100,normed=True,label='predicted C/B w/ GX',color='b',histtype='step')
plt.hist(np.abs(muCpert/muDpert),bins=100,normed=True,label='predicted C/D w/ GX',color='g',histtype='step')

plt.axvline(C_125_1/A_125_1, 0, 1, linestyle='--', linewidth=3, color='r', label='F125 C/A visit 1')
plt.axvline(C_125_1/B_125_1, 0, 1, linestyle='--', linewidth=3, color='b', label='F125 C/B visit 1')
plt.axvline(C_125_1/D_125_1, 0, 1, linestyle='--', linewidth=3, color='g', label='F125 C/D visit 1')
plt.axvline(C_160_1/A_160_1, 0, 1, linestyle='--', linewidth=1, color='r', label='F160 C/A visit 1')
plt.axvline(C_160_1/B_160_1, 0, 1, linestyle='--', linewidth=1, color='b', label='F160 C/B visit 1')
plt.axvline(C_160_1/D_160_1, 0, 1, linestyle='--', linewidth=1, color='g', label='F160 C/D visit 1')
plt.axvline(C_125_2/A_125_2, 0, 1, linestyle='-', linewidth=3, color='r', label='F125 C/A visit 2')
plt.axvline(C_125_2/B_125_2, 0, 1, linestyle='-', linewidth=3, color='b', label='F125 C/B visit 2')
plt.axvline(C_125_2/D_125_2, 0, 1, linestyle='-', linewidth=3, color='g', label='F125 C/D visit 2')
plt.axvline(C_160_2/A_160_2, 0, 1, linestyle='-', linewidth=1, color='r', label='F160 C/A visit 2')
plt.axvline(C_160_2/B_160_2, 0, 1, linestyle='-', linewidth=1, color='b', label='F160 C/B visit 2')
plt.axvline(C_160_2/D_160_2, 0, 1, linestyle='-', linewidth=1, color='g', label='F160 C/D visit 2')

plt.legend(prop={'size': 8}) # font size
plt.tick_params(labelsize=14)
plt.xlabel('Flux ratio', fontsize=16)
plt.ylabel('Normalized distribution', fontsize=16)
plt.savefig('fluxratiosC.png')

plt.clf()

plt.hist(np.abs(muBnopert/muAnopert),bins=100,normed=True,label='predicted B/A w/o GX',color='r',histtype='stepfilled')
plt.hist(np.abs(muBnopert/muCnopert),bins=100,normed=True,label='predicted B/C w/o GX',color='b',histtype='stepfilled')
plt.hist(np.abs(muBnopert/muDnopert),bins=100,normed=True,label='predicted B/D w/o GX',color='g',histtype='stepfilled')
plt.hist(np.abs(muBpert/muApert),bins=100,normed=True,label='predicted B/A w/ GX',color='r',histtype='step')
plt.hist(np.abs(muBpert/muCpert),bins=100,normed=True,label='predicted B/C w/ GX',color='b',histtype='step')
plt.hist(np.abs(muBpert/muDpert),bins=100,normed=True,label='predicted B/D w/ GX',color='g',histtype='step')

plt.axvline(B_125_1/A_125_1, 0, 1, linestyle='--', linewidth=3, color='r', label='F125 B/A visit 1')
plt.axvline(B_125_1/C_125_1, 0, 1, linestyle='--', linewidth=3, color='b', label='F125 B/C visit 1')
plt.axvline(B_125_1/D_125_1, 0, 1, linestyle='--', linewidth=3, color='g', label='F125 B/D visit 1')
plt.axvline(B_160_1/A_160_1, 0, 1, linestyle='--', linewidth=1, color='r', label='F160 B/A visit 1')
plt.axvline(B_160_1/C_160_1, 0, 1, linestyle='--', linewidth=1, color='b', label='F160 B/C visit 1')
plt.axvline(B_160_1/D_160_1, 0, 1, linestyle='--', linewidth=1, color='g', label='F160 B/D visit 1')
plt.axvline(B_125_2/A_125_2, 0, 1, linestyle='-', linewidth=3, color='r', label='F125 B/A visit 2')
plt.axvline(B_125_2/C_125_2, 0, 1, linestyle='-', linewidth=3, color='b', label='F125 B/C visit 2')
plt.axvline(B_125_2/D_125_2, 0, 1, linestyle='-', linewidth=3, color='g', label='F125 B/D visit 2')
plt.axvline(B_160_2/A_160_2, 0, 1, linestyle='-', linewidth=1, color='r', label='F160 B/A visit 2')
plt.axvline(B_160_2/C_160_2, 0, 1, linestyle='-', linewidth=1, color='b', label='F160 B/C visit 2')
plt.axvline(B_160_2/D_160_2, 0, 1, linestyle='-', linewidth=1, color='g', label='F160 B/D visit 2')

plt.legend(prop={'size': 8}) # font size
plt.tick_params(labelsize=14)
plt.xlabel('Flux ratio', fontsize=16)
plt.ylabel('Normalized distribution', fontsize=16)
plt.savefig('fluxratiosB.png')

plt.clf()

plt.hist(np.abs(muDnopert/muAnopert),bins=100,normed=True,label='predicted D/A w/o GX',color='r',histtype='stepfilled')
plt.hist(np.abs(muDnopert/muBnopert),bins=100,normed=True,label='predicted D/B w/o GX',color='b',histtype='stepfilled')
plt.hist(np.abs(muDnopert/muCnopert),bins=100,normed=True,label='predicted D/D w/o GX',color='g',histtype='stepfilled')
plt.hist(np.abs(muDpert/muApert),bins=100,normed=True,label='predicted D/A w/ GX',color='r',histtype='step')
plt.hist(np.abs(muDpert/muBpert),bins=100,normed=True,label='predicted D/B w/ GX',color='b',histtype='step')
plt.hist(np.abs(muDpert/muCpert),bins=100,normed=True,label='predicted D/D w/ GX',color='g',histtype='step')

plt.axvline(D_125_1/A_125_1, 0, 1, linestyle='--', linewidth=3, color='r', label='F125 D/A visit 1')
plt.axvline(D_125_1/B_125_1, 0, 1, linestyle='--', linewidth=3, color='b', label='F125 D/B visit 1')
plt.axvline(D_125_1/C_125_1, 0, 1, linestyle='--', linewidth=3, color='g', label='F125 D/C visit 1')
plt.axvline(D_160_1/A_160_1, 0, 1, linestyle='--', linewidth=1, color='r', label='F160 D/A visit 1')
plt.axvline(D_160_1/B_160_1, 0, 1, linestyle='--', linewidth=1, color='b', label='F160 D/B visit 1')
plt.axvline(D_160_1/C_160_1, 0, 1, linestyle='--', linewidth=1, color='g', label='F160 D/C visit 1')
plt.axvline(D_125_2/A_125_2, 0, 1, linestyle='-', linewidth=3, color='r', label='F125 D/A visit 2')
plt.axvline(D_125_2/B_125_2, 0, 1, linestyle='-', linewidth=3, color='b', label='F125 D/B visit 2')
plt.axvline(D_125_2/C_125_2, 0, 1, linestyle='-', linewidth=3, color='g', label='F125 D/C visit 2')
plt.axvline(D_160_2/A_160_2, 0, 1, linestyle='-', linewidth=1, color='r', label='F160 D/A visit 2')
plt.axvline(D_160_2/B_160_2, 0, 1, linestyle='-', linewidth=1, color='b', label='F160 D/B visit 2')
plt.axvline(D_160_2/C_160_2, 0, 1, linestyle='-', linewidth=1, color='g', label='F160 D/C visit 2')

plt.legend(prop={'size': 8}) # font size
plt.tick_params(labelsize=14)
plt.xlabel('Flux ratio', fontsize=16)
plt.ylabel('Normalized distribution', fontsize=16)
plt.savefig('fluxratiosD.png')


