# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid
# usage: use one of the following arguments: orig or samp, followed by maglimit

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
#from scipy.stats.kde import gaussian_kde
from numpy import linspace

print("Arguments: \n Original values or samples drawn from P(z) and P(Mstar): %s \n Limiting i-band magnitude: %s" % (str(sys.argv[1]), str(sys.argv[2])))

if str(sys.argv[1]) == "samp":
    print "This process is both processor and memory intensive and will take a couple of hours for a sampling of 1000..."
    start_time = time.time()

with open('fieldsforhist50try_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist50try_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist50try_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
    listfields = f.readlines()

with open('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

#with open('fieldsforhist50try_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
#    listfields = f.readlines()

#with open('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#           with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#               outfile.write(infile.read())

#with open('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldsforhist75try_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2]))) as f:
#    listfields = f.readlines()

#with open('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

#with open('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), 'w') as outfile:
#    for i in range(len(listfields)):
#        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
#            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
#                outfile.write(infile.read())

cols=1
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.suptitle(r'HE0435 weight histogram test W1-W4', fontsize=10, y=0.998)



x = linspace(0,2,500)

plt.subplot(451)
rangemax=2
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(451)
#s = "50: %.3f,%.3f,%.3f" % (np.median(q_W1_50_120),np.average(q_W1_50_120),np.median(q_W1_50_120)) # only the peak depends on the binning
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 1

cols=3
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

x = linspace(0,2,500)



plt.subplot(452)
rangemax=2
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(452)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 2

cols=5
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(453)
rangemax=2
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(453)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 3

cols=7
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(454)
rangemax=3
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(454)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 4

cols=9
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(456)
rangemax=6
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(456)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 5

cols=11
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(457)
rangemax=4
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(457)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 6

cols=13
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(458)
rangemax=6
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(458)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 7

cols=15
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(459)
rangemax=6
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(459)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 8

cols=17
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,11)
rangemax=2
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,11)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 9

cols=19
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,12)
rangemax=3
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,12)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 10

cols=21
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,13)
rangemax=6
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,13)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 11

cols=23
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,14)
rangemax=6
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,14)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 12

cols=25
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

x = linspace(0,3,500)



plt.subplot(4,5,16)
rangemax=4
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(x,gauss_q_W1_50(x),'b', linewidth=0.5)
ax=plt.subplot(4,5,16)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 13

cols=27
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

x = linspace(0,3,500)



plt.subplot(4,5,17)
rangemax=5
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,17)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 14

cols=29
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,18)
rangemax=4
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,18)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 15

cols=31
q_read = np.loadtxt('fieldshistW1_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size120_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_120 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size90_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_90 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_50_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_50_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW1_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W1_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW2_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W2_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW3_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W3_75_60 = q_read[q_read < 10]
q_read = np.loadtxt('fieldshistW4_75_%s_size60_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
q_W4_75_60 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_50_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_50_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW1_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W1_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW2_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W2_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW3_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W3_75_45 = q_read[q_read < 10]
#q_read = np.loadtxt('fieldshistW4_75_%s_size45_i%s.lst' % (str(sys.argv[1]),str(sys.argv[2])), usecols=[cols], unpack=True)
#q_W4_75_45 = q_read[q_read < 10]

plt.subplot(4,5,19)
rangemax=7
plt.xlim(25,140)
#plt.scatter(45, np.median(q_W1_50_45), color='b', label='W1_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W2_50_45), color='g', label='W2_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W3_50_45), color='r', label='W3_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W4_50_45), color='k', label='W4_50', linewidth=0.5)
#plt.scatter(45, np.median(q_W1_75_45), color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W2_75_45), color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W3_75_45), color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
#plt.scatter(45, np.median(q_W4_75_45), color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W1_50_60),np.median(q_W1_50_90),np.median(q_W1_50_120)], color='b', label='W1_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W2_50_60),np.median(q_W2_50_90),np.median(q_W2_50_120)], color='g', label='W2_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W3_50_60),np.median(q_W3_50_90),np.median(q_W3_50_120)], color='r', label='W3_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W4_50_60),np.median(q_W4_50_90),np.median(q_W4_50_120)], color='k', label='W4_50', linewidth=0.5)
plt.plot([60,90,120], [np.median(q_W1_75_60),np.median(q_W1_75_90),np.median(q_W1_75_120)], color='b', label='W1_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W2_75_60),np.median(q_W2_75_90),np.median(q_W2_75_120)], color='g', label='W2_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W3_75_60),np.median(q_W3_75_90),np.median(q_W3_75_120)], color='r', label='W3_75', linewidth=0.5, linestyle='dotted')
plt.plot([60,90,120], [np.median(q_W4_75_60),np.median(q_W4_75_90),np.median(q_W4_75_120)], color='k', label='W4_75', linewidth=0.5, linestyle='dotted')
ax=plt.subplot(4,5,19)
#s = "50: %.3f,%.3f,%.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50))
#ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Median", fontsize=7)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 16

plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()

plt.savefig('HE0435sizeplots_%s_i%s.png' % (str(sys.argv[1]),str(sys.argv[2])), dpi=1000)

#plt.show()

os.system("rm fieldshistW1_50_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_50_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_50_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_50_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW1_75_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_75_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_75_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_75_%s_size120_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW1_50_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_50_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_50_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_50_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW1_75_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_75_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_75_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_75_%s_size90_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW1_50_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_50_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_50_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_50_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW1_75_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW2_75_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW3_75_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
os.system("rm fieldshistW4_75_%s_size60_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW1_50_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW2_50_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW3_50_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW4_50_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW1_75_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW2_75_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW3_75_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))
#os.system("rm fieldshistW4_75_%s_size45_i%s.lst" % (str(sys.argv[1]),str(sys.argv[2])))

if str(sys.argv[1]) == "samp":
    print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
