# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid
# usage: use one of the following arguments: orig or samp, followed by number of bins

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

if str(sys.argv[1]) == "samp":
    print "This process is both processor and memory intensive and will take a couple of hours for a sampling of 1000..."
    start_time = time.time()

with open('fieldsforhist50try_%s_i23.lst' % str(sys.argv[1])) as f:
    listfields = f.readlines()

with open('fieldshistW1_50_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_50_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_50_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_50_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_i23.lst' % str(sys.argv[1])) as f:
    listfields = f.readlines()

with open('fieldshistW1_75_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_75_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_75_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_75_i23.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

cols=1
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.suptitle(r'HE0435 weight histogram test W1-W4', fontsize=10, y=0.998)

gauss_q_W1_50 = gaussian_kde(q_W1_50)
gauss_q_W2_50 = gaussian_kde(q_W2_50)
gauss_q_W3_50 = gaussian_kde(q_W3_50)
gauss_q_W4_50 = gaussian_kde(q_W4_50)
gauss_q_W1_75 = gaussian_kde(q_W1_75)
gauss_q_W2_75 = gaussian_kde(q_W2_75)
gauss_q_W3_75 = gaussian_kde(q_W3_75)
gauss_q_W4_75 = gaussian_kde(q_W4_75)

x = linspace(0,2,500)

plt.subplot(451)
rangemax=2
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(451)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_50(x))],np.average(q_W1_50),np.median(q_W1_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_50(x))],np.average(q_W2_50),np.median(q_W2_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_50(x))],np.average(q_W3_50),np.median(q_W3_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_50(x))],np.average(q_W4_50),np.median(q_W4_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_75(x))],np.average(q_W2_75),np.median(q_W2_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_75(x))],np.average(q_W3_75),np.median(q_W3_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_75(x))],np.average(q_W4_75),np.median(q_W4_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 1
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=3
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,2,500)

gauss_q_W1_50 = gaussian_kde(q_W1_50)
gauss_q_W2_50 = gaussian_kde(q_W2_50)
gauss_q_W3_50 = gaussian_kde(q_W3_50)
gauss_q_W4_50 = gaussian_kde(q_W4_50)
gauss_q_W1_75 = gaussian_kde(q_W1_75)
gauss_q_W2_75 = gaussian_kde(q_W2_75)
gauss_q_W3_75 = gaussian_kde(q_W3_75)
gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(452)
rangemax=2
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(452)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_50(x))],np.average(q_W1_50),np.median(q_W1_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_50(x))],np.average(q_W2_50),np.median(q_W2_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_50(x))],np.average(q_W3_50),np.median(q_W3_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_50(x))],np.average(q_W4_50),np.median(q_W4_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_75(x))],np.average(q_W2_75),np.median(q_W2_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_75(x))],np.average(q_W3_75),np.median(q_W3_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_75(x))],np.average(q_W4_75),np.median(q_W4_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 2
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=5
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(453)
rangemax=2
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(453)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 3
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=7
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(454)
rangemax=3
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(454)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 4
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=9
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(456)
rangemax=6
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(456)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 5
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=11
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(457)
rangemax=4
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(457)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 6
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=13
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(458)
rangemax=6
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(458)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 7
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=15
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(459)
rangemax=6
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(459)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 8
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=17
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,11)
rangemax=2
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,11)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 9
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=19
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,12)
rangemax=3
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,12)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 10
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=21
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,13)
rangemax=6
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,13)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 11
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=23
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,14)
rangemax=6
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,14)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 12
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=25
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,3,500)

gauss_q_W1_50 = gaussian_kde(q_W1_50)
gauss_q_W2_50 = gaussian_kde(q_W2_50)
gauss_q_W3_50 = gaussian_kde(q_W3_50)
gauss_q_W4_50 = gaussian_kde(q_W4_50)
gauss_q_W1_75 = gaussian_kde(q_W1_75)
gauss_q_W2_75 = gaussian_kde(q_W2_75)
gauss_q_W3_75 = gaussian_kde(q_W3_75)
gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(4,5,16)
rangemax=4
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(4,5,16)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_50(x))],np.average(q_W1_50),np.median(q_W1_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_50(x))],np.average(q_W2_50),np.median(q_W2_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_50(x))],np.average(q_W3_50),np.median(q_W3_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_50(x))],np.average(q_W4_50),np.median(q_W4_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_75(x))],np.average(q_W2_75),np.median(q_W2_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_75(x))],np.average(q_W3_75),np.median(q_W3_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_75(x))],np.average(q_W4_75),np.median(q_W4_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 13
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=27
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

x = linspace(0,3,500)

gauss_q_W1_50 = gaussian_kde(q_W1_50)
gauss_q_W2_50 = gaussian_kde(q_W2_50)
gauss_q_W3_50 = gaussian_kde(q_W3_50)
gauss_q_W4_50 = gaussian_kde(q_W4_50)
gauss_q_W1_75 = gaussian_kde(q_W1_75)
gauss_q_W2_75 = gaussian_kde(q_W2_75)
gauss_q_W3_75 = gaussian_kde(q_W3_75)
gauss_q_W4_75 = gaussian_kde(q_W4_75)

plt.subplot(4,5,17)
rangemax=5
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,17)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_50(x))],np.average(q_W1_50),np.median(q_W1_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_50(x))],np.average(q_W2_50),np.median(q_W2_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_50(x))],np.average(q_W3_50),np.median(q_W3_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "50: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_50(x))],np.average(q_W4_50),np.median(q_W4_50))
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W2_75(x))],np.average(q_W2_75),np.median(q_W2_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W3_75(x))],np.average(q_W3_75),np.median(q_W3_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W4_75(x))],np.average(q_W4_75),np.median(q_W4_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='r',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 14
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=29
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,18)
rangemax=4
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,18)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 15
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

cols=31
q_W1_50read = np.loadtxt('fieldshistW1_50_i23.lst', usecols=[cols], unpack=True)
q_W1_50 = q_W1_50read[q_W1_50read < 10]
q_W2_50read = np.loadtxt('fieldshistW2_50_i23.lst', usecols=[cols], unpack=True)
q_W2_50 = q_W2_50read[q_W2_50read < 10]
q_W3_50read = np.loadtxt('fieldshistW3_50_i23.lst', usecols=[cols], unpack=True)
q_W3_50 = q_W3_50read[q_W3_50read < 10]
q_W4_50read = np.loadtxt('fieldshistW4_50_i23.lst', usecols=[cols], unpack=True)
q_W4_50 = q_W4_50read[q_W4_50read < 10]
q_W1_75read = np.loadtxt('fieldshistW1_75_i23.lst', usecols=[cols], unpack=True)
q_W1_75 = q_W1_75read[q_W1_75read < 10]
q_W2_75read = np.loadtxt('fieldshistW2_75_i23.lst', usecols=[cols], unpack=True)
q_W2_75 = q_W2_75read[q_W2_75read < 10]
q_W3_75read = np.loadtxt('fieldshistW3_75_i23.lst', usecols=[cols], unpack=True)
q_W3_75 = q_W3_75read[q_W3_75read < 10]
q_W4_75read = np.loadtxt('fieldshistW4_75_i23.lst', usecols=[cols], unpack=True)
q_W4_75 = q_W4_75read[q_W4_75read < 10]

plt.subplot(4,5,19)
rangemax=7
n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax])
plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[2]), range=[0, rangemax], linestyle='dotted')
ax=plt.subplot(4,5,19)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
ax.text(0.15, 0.6, s, fontsize=5, color='g',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
ax.text(0.15, 0.7, s, fontsize=5, color='b',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
ax.text(0.15, 0.5, s, fontsize=5, color='g',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 16
print "finished subplot %d/16; fraction of points inside the q < 10 cut: W1_50 %.3f W1_75 %.3f W2_50 %.3f W2_75 %.3f W3_50 %.3f W3_75 %.3f W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W1_50.size)/q_W1_50read.size, float(q_W1_75.size)/q_W1_75read.size, float(q_W2_50.size)/q_W2_50read.size, float(q_W2_75.size)/q_W2_75read.size, float(q_W3_50.size)/q_W3_50read.size, float(q_W3_75.size)/q_W3_75read.size, float(q_W4_50.size)/q_W4_50read.size, float(q_W4_75.size)/q_W4_75read.size)

plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()


plt.savefig('HE0435overdensities_%s_i23.png' % str(sys.argv[1]), dpi=1000)

#plt.show()

os.system("rm fieldshistW1_50_i23.lst")
os.system("rm fieldshistW2_50_i23.lst")
os.system("rm fieldshistW3_50_i23.lst")
os.system("rm fieldshistW4_50_i23.lst")
os.system("rm fieldshistW1_75_i23.lst")
os.system("rm fieldshistW2_75_i23.lst")
os.system("rm fieldshistW3_75_i23.lst")
os.system("rm fieldshistW4_75_i23.lst")

if str(sys.argv[1]) == "samp":
    print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
