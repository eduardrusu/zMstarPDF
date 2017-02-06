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

print("Arguments: \n Lens field: %s \n Original values or samples drawn from P(z) and P(Mstar): %s \n Number of bins: %s \n Radius of each cell: %s \n Limiting i-band magnitude: %s \n" % (str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5])))

if (str(sys.argv[2]) == "samp") or (str(sys.argv[2]) == "tab"):
    print "This process is both processor and memory intensive and will take a couple of hours for a sampling of 1000..."
    start_time = time.time()

with open('fieldsforhist50try_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist50try_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5]))) as f:
    listfields = f.readlines()

with open('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

cols=1

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.suptitle(r'HE0435 weight histogram test W1-W4', fontsize=10, y=0.998)


gauss_q_W4_50old = gaussian_kde(q_W4_50old)

gauss_q_W4_75old = gaussian_kde(q_W4_75old)

gauss_q_W4_50new = gaussian_kde(q_W4_50new)

gauss_q_W4_75new = gaussian_kde(q_W4_75new)

x = linspace(0,2,500)

plt.subplot(451)
rangemax=4

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(451)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75old[np.argmax(n_q_W4_75old)],np.average(q_W4_75old),np.median(q_W4_75old))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50new[np.argmax(n_q_W4_50new)],np.average(q_W4_50new),np.median(q_W4_50new))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 1
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50new.size)/q_W4_50readnew.size, float(q_W4_75new.size)/q_W4_75readnew.size)

cols=3

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

x = linspace(0,2,500)


gauss_q_W4_50old = gaussian_kde(q_W4_50old)

gauss_q_W4_75old = gaussian_kde(q_W4_75old)

gauss_q_W4_50new = gaussian_kde(q_W4_50new)

gauss_q_W4_75new = gaussian_kde(q_W4_75new)

plt.subplot(452)
rangemax=4

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='k', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='k', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='k', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='k', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(452)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75old[np.argmax(n_q_W4_75old)],np.average(q_W4_75old),np.median(q_W4_75old))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50new[np.argmax(n_q_W4_50new)],np.average(q_W4_50new),np.median(q_W4_50new))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 2
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=5

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(453)
rangemax=4

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(453)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 3
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=7

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(454)
rangemax=4

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(454)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 4
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=9

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(456)
rangemax=7

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(456)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 5
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=11

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(457)
rangemax=6

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(457)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 6
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=13

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(458)
rangemax=6

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(458)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 7
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=15

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(459)
rangemax=7

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(459)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 8
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=17

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,11)
rangemax=4

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,11)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 9
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=19

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,12)
rangemax=7

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,12)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 10
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=21

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,13)
rangemax=6

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,13)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 11
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=23

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,14)
rangemax=6

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,14)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 12
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=25

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

x = linspace(0,3,500)


gauss_q_W4_50old = gaussian_kde(q_W4_50old)

gauss_q_W4_75old = gaussian_kde(q_W4_75old)

gauss_q_W4_50new = gaussian_kde(q_W4_50new)

gauss_q_W4_75new = gaussian_kde(q_W4_75new)

plt.subplot(4,5,16)
rangemax=8

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

#plt.plot(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(4,5,16)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75old[np.argmax(n_q_W4_75old)],np.average(q_W4_75old),np.median(q_W4_75old))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50new[np.argmax(n_q_W4_50new)],np.average(q_W4_50new),np.median(q_W4_50new))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 13
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=27

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

x = linspace(0,3,500)


gauss_q_W4_50old = gaussian_kde(q_W4_50old)

gauss_q_W4_75old = gaussian_kde(q_W4_75old)

gauss_q_W4_50new = gaussian_kde(q_W4_50new)

gauss_q_W4_75new = gaussian_kde(q_W4_75new)

plt.subplot(4,5,17)
rangemax=8

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,17)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75old[np.argmax(n_q_W4_75old)],np.average(q_W4_75old),np.median(q_W4_75old))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50new[np.argmax(n_q_W4_50new)],np.average(q_W4_50new),np.median(q_W4_50new))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)
#s = "75: %.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_W1_75(x))],np.average(q_W1_75),np.median(q_W1_75))
s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 14
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=29

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,18)
rangemax=7

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,18)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 15
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

cols=31

q_W4_50readnew = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50new = q_W4_50readnew[q_W4_50readnew < 10]

q_W4_75readnew = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_new.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75new = q_W4_75readnew[q_W4_75readnew < 10]

q_W4_50readold = np.loadtxt('fieldshistW4_50_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_50old = q_W4_50readold[q_W4_50readold < 10]

q_W4_75readold = np.loadtxt('fieldshistW4_75_%s_%s_size%s_i%s_old.lst' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), usecols=[cols], unpack=True)
q_W4_75old = q_W4_75readold[q_W4_75readold < 10]

plt.subplot(4,5,19)
rangemax=10

n_q_W4_50old, bins_q_W4_50old, patches = plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_50new, bins_q_W4_50new, patches = plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

n_q_W4_75old, bins_q_W4_75old, patches = plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

n_q_W4_75new, bins_q_W4_75new, patches = plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_50old, histtype='step', color='k', label='W4_50old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_50new, histtype='step', color='r', label='W4_50new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])

plt.hist(q_W4_75old, histtype='step', color='k', label='W4_75old', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

plt.hist(q_W4_75new, histtype='step', color='r', label='W4_75new', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')

ax=plt.subplot(4,5,19)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.4, s, fontsize=5, color='r',transform=ax.transAxes)

s = "50: %.3f,%.3f,%.3f" % (bins_q_W4_50old[np.argmax(n_q_W4_50old)],np.average(q_W4_50old),np.median(q_W4_50old))
ax.text(0.15, 0.3, s, fontsize=5, color='r',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.2, s, fontsize=5, color='k',transform=ax.transAxes)

s = "75: %.3f,%.3f,%.3f" % (bins_q_W4_75new[np.argmax(n_q_W4_75new)],np.average(q_W4_75new),np.median(q_W4_75new))
ax.text(0.15, 0.1, s, fontsize=5, color='k',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 16
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n  W4_50 %.3f W4_75 %.3f " % (subplot, float(q_W4_50old.size)/q_W4_50readold.size, float(q_W4_75old.size)/q_W4_75readold.size)

plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()

plt.savefig('%s_overdensities_%s_size%s_i%s_oldnew.png' % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])), dpi=1000)

#plt.show()

os.system("rm fieldshistW4_50_%s_%s_size%s_i%s_new.lst" % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])))
os.system("rm fieldshistW4_75_%s_%s_size%s_i%s_new.lst" % (str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[4]),str(sys.argv[5])))


if str(sys.argv[2]) == "samp":
    print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
