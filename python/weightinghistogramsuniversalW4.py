# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid
# usage: use one of the following arguments: lens name, followed by orig or samp, followed by number of bins, followed by radius (45,60,90 or 120) and by maglimit

import numpy as np
import sys
import os
from os import system
#import scipy
#from scipy import special
#from astropy.io import fits
#from astropy.wcs import WCS
#from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.io import ascii
#from astropy.table import Table, Column
import time
import matplotlib.pyplot as plt
#from numpy.random import normal
from scipy.stats.kde import gaussian_kde
from numpy import linspace

print("Arguments: \n Lens field: %s \n Original values or samples drawn from P(z) and P(Mstar): %s \n Number of bins: %s \n Radius of each cell: %s \n Limiting i-band magnitude: %s \n Classification: %s" % (str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6])))

if (str(sys.argv[2]) == "samp") or (str(sys.argv[2]) == "tab"):
    print "This process is both processor and memory intensive and will take a couple of hours for a sampling of 1000..."
    start_time = time.time()

with open('fieldsforhist50try_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

cols=1

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.suptitle(r'HE0435 weight histogram test W1-W4', fontsize=10, y=0.998)


gauss_q_W4_50 = gaussian_kde(q_W4_50)

gauss_q_W4_75 = gaussian_kde(q_W4_75)

x = linspace(0,2,500)

plt.subplot(451)
rangemax=4

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(451)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 1
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=3

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,2,500)


gauss_q_W4_50 = gaussian_kde(q_W4_50)

gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(452)
rangemax=4

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(452)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 2
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=5

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(453)
rangemax=4

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(453)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 3
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=7

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(454)
rangemax=4

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(454)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 4
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=9

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(456)
rangemax=7

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(456)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 5
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=11

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(457)
rangemax=6

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(457)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 6
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=13

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(458)
rangemax=6

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(458)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 7
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=15

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(459)
rangemax=7

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(459)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 8
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=17

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,11)
rangemax=4

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,11)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 9
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=19

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,12)
rangemax=7

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,12)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 10
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=21

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,13)
rangemax=6

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,13)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 11
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=23

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,14)
rangemax=6

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,14)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 12
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=25

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,3,500)


gauss_q_W4_50 = gaussian_kde(q_W4_50)

gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(4,5,16)
rangemax=8

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(4,5,16)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 13
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=27

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,3,500)


gauss_q_W4_50 = gaussian_kde(q_W4_50)

gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(4,5,17)
rangemax=8

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,17)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 14
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=29

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,18)
rangemax=7

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,18)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 15
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=31

q_W4_50read = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]

q_W4_75read = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_%s.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,19)
rangemax=10

n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,19)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 16
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()

plt.savefig('%s_overdensities_%s_size%s_i%s_%s.png' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])), dpi=1000)

#plt.show()

os.system("rm fieldshistW4_50_%s_%s_size%s_i%s_%s.lst" % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])))
os.system("rm fieldshistW4_75_%s_%s_size%s_i%s_%s.lst" % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6])))

if str(sys.argv[2]) == "samp":
    print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
