# The code uses the weighted count ratios derived by weightinguniversal_overlap_sampling_nobeta_WFI2033rethought.py to produce histograms and compute statistics
# run as python weightinguniversal_histograms_samples_WFI2033.py WFI2033 45 23 meds bpz deti IRAC 5 100

import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt

plt.clf()

fontabsciss = 10
fontlabel = 4
pltrange = 3
samples = 10
limit = 10**30
root = "/Volumes/perseus_1/CFHTLens_galphotmstar/"

lens = str(sys.argv[1])
radius = str(sys.argv[2])
mag = str(sys.argv[3])
mode = str(sys.argv[4])
photz = str(sys.argv[5])
detect = str(sys.argv[6])
irac = str(sys.argv[7])
inner = str(sys.argv[8])
bin = int(str(sys.argv[9]))

start_time = time.time()

print "Working on samples:"

medsum50W1 = np.zeros((18,samples))
medsum75W1 = np.zeros((18,samples))
medsum50W2 = np.zeros((18,samples))
medsum75W2 = np.zeros((18,samples))
medsum50W3 = np.zeros((18,samples))
medsum75W3 = np.zeros((18,samples))
medsum50W4 = np.zeros((18,samples))
medsum75W4 = np.zeros((18,samples))

for nr in range(samples):
    print '%s...' %nr
    lstW1_50 = [x for x in os.listdir(root) if ('W1' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)] # select from the files in the root directory
    lstW1_75 = [x for x in os.listdir(root) if ('W1' in x) and ('_24galphotmstar_75_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW2_50 = [x for x in os.listdir(root) if ('W2' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW2_75 = [x for x in os.listdir(root) if ('W2' in x) and ('_24galphotmstar_75_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW3_50 = [x for x in os.listdir(root) if ('W3' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW3_75 = [x for x in os.listdir(root) if ('W3' in x) and ('_24galphotmstar_75_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW4_50 = [x for x in os.listdir(root) if ('W4' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW4_75 = [x for x in os.listdir(root) if ('W4' in x) and ('_24galphotmstar_75_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]

    if mag == "24" and photz == "bpz": cols=[4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38]
    if mag == "24" and photz == "eazy": cols=[40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74]
    if mag == "23" and photz == "bpz": cols=[5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]
    if mag == "23" and photz == "eazy": cols=[41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75]
                
    print "W1..."
    for i in range(len(lstW1_50)):
        if i == 0:
            q_W1_50read = np.loadtxt(root+lstW1_50[i], usecols=cols, unpack=True)
        else:
            q_W1_50read = np.r_['1',q_W1_50read,np.loadtxt(root+ lstW1_50[i], usecols=cols, unpack=True)]
    for i in range(len(lstW1_75)):
        if i == 0:
            q_W1_75read = np.loadtxt(root+lstW1_75[i], usecols=cols, unpack=True)
        else:
            q_W1_75read = np.r_['1',q_W1_75read,np.loadtxt(root+ lstW1_75[i], usecols=cols, unpack=True)]
            
    print "W2..."
    for i in range(len(lstW2_50)):
        if i == 0:
            q_W2_50read = np.loadtxt(root+lstW2_50[i], usecols=cols, unpack=True)
        else:
            q_W2_50read = np.r_['1',q_W2_50read,np.loadtxt(root+ lstW2_50[i], usecols=cols, unpack=True)]
    for i in range(len(lstW2_75)):
        if i == 0:
            q_W2_75read = np.loadtxt(root+lstW2_75[i], usecols=cols, unpack=True)
        else:
            q_W2_75read = np.r_['1',q_W2_75read,np.loadtxt(root+ lstW2_75[i], usecols=cols, unpack=True)]
          
    print "W3..."
    for i in range(len(lstW3_50)):
        if i == 0:
            q_W3_50read = np.loadtxt(root+lstW3_50[i], usecols=cols, unpack=True)
        else:
            q_W3_50read = np.r_['1',q_W3_50read,np.loadtxt(root+ lstW3_50[i], usecols=cols, unpack=True)]
    for i in range(len(lstW3_75)):
        if i == 0:
            q_W3_75read = np.loadtxt(root+lstW3_75[i], usecols=cols, unpack=True)
        else:
            q_W3_75read = np.r_['1',q_W3_75read,np.loadtxt(root+ lstW3_75[i], usecols=cols, unpack=True)]
    
    print "W4..."
    for i in range(len(lstW4_50)):
        if i == 0:
            q_W4_50read = np.loadtxt(root+lstW4_50[i], usecols=cols, unpack=True)
        else:
            q_W4_50read = np.r_['1',q_W4_50read,np.loadtxt(root+ lstW4_50[i], usecols=cols, unpack=True)]
    for i in range(len(lstW4_75)):
        if i == 0:
            q_W4_75read = np.loadtxt(root+lstW4_75[i], usecols=cols, unpack=True)
        else:
            q_W4_75read = np.r_['1',q_W4_75read,np.loadtxt(root+ lstW4_75[i], usecols=cols, unpack=True)]
    
    for j in range(18):
        q_W1_50 = q_W1_50read[j][q_W1_50read[j] < limit]
        if mode == "sum": q_W1_50 = abs(q_W1_50) # fix the negative halo convergence
        q_W1_75 = q_W1_75read[j][q_W1_75read[j] < limit]
        if mode == "sum": q_W1_75 = abs(q_W1_75)
        q_W2_50 = q_W2_50read[j][q_W2_50read[j] < limit]
        if mode == "sum": q_W2_50 = abs(q_W2_50) # fix the negative halo convergence
        q_W2_75 = q_W2_75read[j][q_W2_75read[j] < limit]
        if mode == "sum": q_W2_75 = abs(q_W2_75)
        q_W3_50 = q_W3_50read[j][q_W3_50read[j] < limit]
        if mode == "sum": q_W3_50 = abs(q_W3_50) # fix the negative halo convergence
        q_W3_75 = q_W3_75read[j][q_W3_75read[j] < limit]
        if mode == "sum": q_W3_75 = abs(q_W3_75)
        q_W4_50 = q_W4_50read[j][q_W4_50read[j] < limit]
        if mode == "sum": q_W4_50 = abs(q_W4_50) # fix the negative halo convergence
        q_W4_75 = q_W4_75read[j][q_W4_75read[j] < limit]
        if mode == "sum": q_W4_75 = abs(q_W4_75)

        if mode == "sum":
            medsum50W1[j][nr] = np.average(q_W1_50)
            medsum75W1[j][nr] = np.average(q_W1_75)
            medsum50W2[j][nr] = np.average(q_W2_50)
            medsum75W2[j][nr] = np.average(q_W2_75)
            medsum50W3[j][nr] = np.average(q_W3_50)
            medsum75W3[j][nr] = np.average(q_W3_75)
            medsum50W4[j][nr] = np.average(q_W4_50)
            medsum75W4[j][nr] = np.average(q_W4_75)
        if mode == "meds":
            medsum50W1[j][nr] = np.median(q_W1_50)
            medsum75W1[j][nr] = np.median(q_W1_75)
            medsum50W2[j][nr] = np.median(q_W2_50)
            medsum75W2[j][nr] = np.median(q_W2_75)
            medsum50W3[j][nr] = np.median(q_W3_50)
            medsum75W3[j][nr] = np.median(q_W3_75)
            medsum50W4[j][nr] = np.median(q_W4_50)
            medsum75W4[j][nr] = np.median(q_W4_75)
                
std50W1 = np.zeros(18)
std75W1 = np.zeros(18)
std50W2 = np.zeros(18)
std75W2 = np.zeros(18)
std50W3 = np.zeros(18)
std75W3 = np.zeros(18)
std50W4 = np.zeros(18)
std75W4 = np.zeros(18)

for i in range(18):
    std50W1[i] = np.std(medsum50W1[i],ddof=1) # ddof=1 stands for sample standard deviation, not population standard deviation
    std75W1[i] = np.std(medsum75W1[i],ddof=1)
    std50W2[i] = np.std(medsum50W2[i],ddof=1)
    std75W2[i] = np.std(medsum75W2[i],ddof=1)
    std50W3[i] = np.std(medsum50W3[i],ddof=1)
    std75W3[i] = np.std(medsum75W3[i],ddof=1)
    std50W4[i] = np.std(medsum50W4[i],ddof=1)
    std75W4[i] = np.std(medsum75W4[i],ddof=1)

print "Plotting..."

plt.suptitle(r'%s weighted counts histogram W1-W4 %s %s %s %s %s arcsec %s %s' % (lens, radius, mag, mode, photz, irac, inner, detect), fontsize=10, y=0.998)

for i in range(18):

    #n_q_W1_50, bins_q_W1_50, patches = plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax]) # in case I want to compute the peak
    #n_q_W2_50, bins_q_W2_50, patches = plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])
    #n_q_W3_50, bins_q_W3_50, patches = plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])
    #n_q_W4_50, bins_q_W4_50, patches = plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax])
    #n_q_W1_75, bins_q_W1_75, patches = plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
    #n_q_W2_75, bins_q_W2_75, patches = plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
    #n_q_W3_75, bins_q_W3_75, patches = plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
    #n_q_W4_75, bins_q_W4_75, patches = plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=int(sys.argv[3]), range=[0, rangemax], linestyle='dotted')
    
    if i == 0: ax=plt.subplot(5,4,1)
    if i == 1: ax=plt.subplot(5,4,2)
    if i == 2: ax=plt.subplot(5,4,3)
    if i == 3: ax=plt.subplot(5,4,4)
    if i == 4: ax=plt.subplot(5,4,5)
    if i == 5: ax=plt.subplot(5,4,6)
    if i == 6: ax=plt.subplot(5,4,7)
    if i == 7: ax=plt.subplot(5,4,8)
    if i == 8: ax=plt.subplot(5,4,9)
    if i == 9: ax=plt.subplot(5,4,10)
    if i == 10: ax=plt.subplot(5,4,11)
    if i == 11: ax=plt.subplot(5,4,12)
    if i == 12: ax=plt.subplot(5,4,13)
    if i == 13: ax=plt.subplot(5,4,14)
    if i == 14: ax=plt.subplot(5,4,15)
    if i == 15: ax=plt.subplot(5,4,17)
    if i == 16: ax=plt.subplot(5,4,18)
    if i == 17: ax=plt.subplot(5,4,19)
    
    q_W1_50 = q_W1_50read[i][q_W1_50read[i] < limit]
    if mode == "sum": q_W1_50 = abs(q_W1_50) # fix the negative halo convergence
    q_W1_75 = q_W1_75read[i][q_W1_75read[i] < limit]
    if mode == "sum": q_W1_75 = abs(q_W1_75)
    q_W2_50 = q_W2_50read[i][q_W2_50read[i] < limit]
    if mode == "sum": q_W2_50 = abs(q_W2_50) # fix the negative halo convergence
    q_W2_75 = q_W2_75read[i][q_W2_75read[i] < limit]
    if mode == "sum": q_W2_75 = abs(q_W2_75)
    q_W3_50 = q_W3_50read[i][q_W3_50read[i] < limit]
    if mode == "sum": q_W3_50 = abs(q_W3_50) # fix the negative halo convergence
    q_W3_75 = q_W3_75read[i][q_W3_75read[i] < limit]
    if mode == "sum": q_W3_75 = abs(q_W3_75)
    q_W4_50 = q_W4_50read[i][q_W4_50read[i] < limit]
    if mode == "sum": q_W4_50 = abs(q_W4_50) # fix the negative halo convergence
    q_W4_75 = q_W4_75read[i][q_W4_75read[i] < limit]
    if mode == "sum": q_W4_75 = abs(q_W4_75)
    
    plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
                
    #s = "50: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W1_50[np.argmax(n_q_W1_50)],np.average(q_W1_50),np.median(q_W1_50)) # only the peak depends on the binning
    ##ax.text(0.05, 0.8, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
    #s = "50: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W2_50[np.argmax(n_q_W2_50)],np.average(q_W2_50),np.median(q_W2_50))
    ##ax.text(0.05, 0.6, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
    #s = "50: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W3_50[np.argmax(n_q_W3_50)],np.average(q_W3_50),np.median(q_W3_50))
    ##ax.text(0.05, 0.4, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
    #s = "50: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W4_50[np.argmax(n_q_W4_50)],np.average(q_W4_50),np.median(q_W4_50))
    ##ax.text(0.05, 0.2, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
    #s = "75: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W1_75[np.argmax(n_q_W1_75)],np.average(q_W1_75),np.median(q_W1_75))
    ##ax.text(0.05, 0.7, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
    #s = "75: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W2_75[np.argmax(n_q_W2_75)],np.average(q_W2_75),np.median(q_W2_75))
    ##ax.text(0.05, 0.5, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
    #s = "75: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W3_75[np.argmax(n_q_W3_75)],np.average(q_W3_75),np.median(q_W3_75))
    ##ax.text(0.05, 0.3, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
    #s = "75: peak = %.3f ave = %.3f med = %.3f" % (bins_q_W4_75[np.argmax(n_q_W4_75)],np.average(q_W4_75),np.median(q_W4_75))
    ##ax.text(0.05, 0.1, s, fontsize=fontlabel, color='k',transform=ax.transAxes)

    if mode == "sum":
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.average(medsum50W1[i]),std50W1[i],np.average(medsum75W1[i]),std75W1[i])
        ax.text(0.05, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.average(medsum50W2[i]),std50W2[i],np.average(medsum75W2[i]),std75W2[i])
        ax.text(0.05, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.average(medsum50W3[i]),std50W3[i],np.average(medsum75W3[i]),std75W3[i])
        ax.text(0.05, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.average(medsum50W4[i]),std50W4[i],np.average(medsum75W4[i]),std75W4[i])
        ax.text(0.05, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 50: %.3f +/- %.3f" % (np.average([np.average(medsum50W1[i]),np.average(medsum50W2[i]),np.average(medsum50W3[i]),np.average(medsum50W4[i])]),np.sqrt(np.std([np.average(medsum50W1[i]),np.average(medsum50W2[i]),np.average(medsum50W3[i]),np.average(medsum50W4[i])],ddof=1)**2 + np.average([std50W1[i],std50W2[i],std50W3[i],std50W4[i]])**2))  # the errors include in quadrature the scatter from 10 samplings, and from W1-4, for 50%
        ax.text(0.05, 0.1, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 75: %.3f +/- %.3f" % (np.average([np.average(medsum75W1[i]),np.average(medsum75W2[i]),np.average(medsum75W3[i]),np.average(medsum75W4[i])]),np.sqrt(np.std([np.average(medsum75W1[i]),np.average(medsum75W2[i]),np.average(medsum75W3[i]),np.average(medsum75W4[i])],ddof=1)**2 + np.average([std75W1[i],std75W2[i],std75W3[i],std75W4[i]])**2))
        ax.text(0.05, 0.0, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
                                                
    if mode == "meds":
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.median(medsum50W1[i]),std50W1[i],np.median(medsum75W1[i]),std75W1[i])
        ax.text(0.05, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.median(medsum50W2[i]),std50W2[i],np.median(medsum75W2[i]),std75W2[i])
        ax.text(0.05, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.median(medsum50W3[i]),std50W3[i],np.median(medsum75W3[i]),std75W3[i])
        ax.text(0.05, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "50: %.3f +/- %.3f 75: %.3f +/- %.3f" % (np.median(medsum50W4[i]),std50W4[i],np.median(medsum75W4[i]),std75W4[i])
        ax.text(0.05, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 50: %.3f +/- %.3f" % (np.average([np.median(medsum50W1[i]),np.median(medsum50W2[i]),np.median(medsum50W3[i]),np.median(medsum50W4[i])]),np.sqrt(np.std([np.median(medsum50W1[i]),np.median(medsum50W2[i]),np.median(medsum50W3[i]),np.median(medsum50W4[i])],ddof=1)**2 + np.average([std50W1[i],std50W2[i],std50W3[i],std50W4[i]])**2))  # the errors include in quadrature the scatter from 10 samplings, and from W1-4, for 50%
        ax.text(0.05, 0.1, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 75: %.3f +/- %.3f" % (np.average([np.median(medsum75W1[i]),np.median(medsum75W2[i]),np.median(medsum75W3[i]),np.median(medsum75W4[i])]),np.sqrt(np.std([np.median(medsum75W1[i]),np.median(medsum75W2[i]),np.median(medsum75W3[i]),np.median(medsum75W4[i])],ddof=1)**2 + np.average([std75W1[i],std75W2[i],std75W3[i],std75W4[i]])**2))
        ax.text(0.05, 0.0, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
                
    if i == 0: plt.xlabel(r'$\zeta_{gal}$', fontsize=fontabsciss)
    if i == 1: plt.xlabel(r'$\zeta_{z}$', fontsize=fontabsciss)
    if i == 2: plt.xlabel(r'$\zeta_{M_\star}$', fontsize=fontabsciss)
    if i == 3: plt.xlabel(r'$\zeta_{M^2_\star}$', fontsize=fontabsciss)
    if i == 4: plt.xlabel(r'$\zeta_{M^3_\star}$', fontsize=fontabsciss)
    if i == 5: plt.xlabel(r'$\zeta_{1/r}$', fontsize=fontabsciss)
    if i == 6: plt.xlabel(r'$\zeta_{z/r}$', fontsize=fontabsciss)
    if i == 7: plt.xlabel(r'$\zeta_{M_\star/r}$', fontsize=fontabsciss)
    if i == 8: plt.xlabel(r'$\zeta_{M^2_\star/r}$', fontsize=fontabsciss)
    if i == 9: plt.xlabel(r'$\zeta_{M^3_\star/r}$', fontsize=fontabsciss)
    if i == 10: plt.xlabel(r'$\zeta_{M^2_{\star\mathrm{,rms}}}$', fontsize=fontabsciss)
    if i == 11: plt.xlabel(r'$\zeta_{M^3_{\star\mathrm{,rms}}}$', fontsize=fontabsciss)
    if i == 12: plt.xlabel(r'$\zeta_{M^2_{\star\mathrm{,rms}}/r}$', fontsize=fontabsciss)
    if i == 13: plt.xlabel(r'$\zeta_{M^3_{\star\mathrm{,rms}}/r}$', fontsize=fontabsciss)
    if i == 14: plt.xlabel(r'$\zeta_\mathrm{flexion}$', fontsize=fontabsciss)
    if i == 15: plt.xlabel(r'$\zeta_\mathrm{tidal}$', fontsize=fontabsciss)
    if i == 16: plt.xlabel(r'$\zeta_\mathrm{SIS}$', fontsize=fontabsciss)
    if i == 17: plt.xlabel(r'$\zeta_\mathrm{SIShalo}$', fontsize=fontabsciss)
    if i in [0,4,8,12,15]:
        plt.ylabel("Normalized counts", fontsize=5)
    plt.tick_params(axis='x', labelsize=4)
    plt.tick_params(axis='y', labelsize=4)
    plt.setp(plt.xticks()[1], rotation=90)
    subplot = i+1
    print "finished subplot %d/18; fraction of points inside the < %s cut: \n W1_50 %.3f W1_75 %.3f \n W2_50 %.3f W2_75 %.3f \n W3_50 %.3f W3_75 %.3f \n W4_50 %.3f W4_75 %.3f " % (subplot, limit, float(q_W1_50.size)/q_W1_50read[0].size, float(q_W1_75.size)/q_W1_75read[0].size, float(q_W2_50.size)/q_W2_50read[0].size, float(q_W2_75.size)/q_W2_75read[0].size, float(q_W3_50.size)/q_W3_50read[0].size, float(q_W3_75.size)/q_W3_75read[0].size, float(q_W4_50.size)/q_W4_50read[0].size, float(q_W4_75.size)/q_W4_75read[0].size)

plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.6)
plt.subplot(5,4,5)
plt.legend(bbox_to_anchor=(5, -5), loc='lower right', borderaxespad=0., fontsize=10)
plt.savefig('%s%s_weightedcountshist_%s_%s_%s_%s_%s_%s_%sarcsec.png' % (root, lens, radius, mag, mode, photz, detect, irac, inner), dpi=500)
                
print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
