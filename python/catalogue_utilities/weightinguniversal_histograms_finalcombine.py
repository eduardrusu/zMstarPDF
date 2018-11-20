# CE Rusu, Feb 13 2018
# Combines the results produced by weightinguniversal_histograms_samples.py into a final text file with the final distributions and widths
# "local" only outputs the median for the fiducial, but the std is still the one of global
# run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_histograms_finalcombine.py WFI2033 meds 5 -1.0 -1.0 global 22.5 removehandpicked
# It can be used to compute not only standard deviation around the global mean ("global"), but around a given ("local") point as well

import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt

lens = str(sys.argv[1])
mode = str(sys.argv[2])
inner = str(sys.argv[3])
zinf = str(sys.argv[4])
zsup = str(sys.argv[5])
local = str(sys.argv[6])
maglim = str(sys.argv[7])
try: handpicked = '_'+str(sys.argv[8])
except: handpicked = ''
rootin = "/Volumes/LaCieSubaru/weightedcounts/%s/" % lens
rootout = "/Users/cerusu/Dropbox/Davis_work/code/%s" % lens

# select the desired files
# edit the conditions as desired, to restrict the included criteria:
if handpicked == '':
    #lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('4samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' not in x)) and (('IRAC' in x) | ('noIRAC' in x)) (('_detir_' in x) | ('_deti_' in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('21samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' not in x)) and (('IRAC' in x) | ('noIRAC' in x)) (('_detir_' in x) | ('_deti_' in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    #lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('4samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
    lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('21samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
    #lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('1samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) & ('eazy' not in x)) and (('IRAC' in x) & ('noIRAC' not in x)) (('_detir_' in x) & ('_deti_' not in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    #lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('1samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('handpicked' not in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) & ('eazy' not in x)) and (('IRAC' in x) & ('noIRAC' not in x)) and (('_detir_' in x) & ('_deti_' not in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
else:
    #lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('4samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('21samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    #lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('4samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
    lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('21samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) | ('eazy' in x)) and (('IRAC' in x) | ('noIRAC' in x)) and (('_detir_' in x) | ('_deti_' in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
    #lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('1samples' in x) and ('45arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) & ('eazy' not in x)) and (('IRAC' in x) & ('noIRAC' not in x)) and (('_detir_' in x) & ('_deti_' not in x)) ] #  # CHOOSE WHAT YOU WANT HERE
    #lst120 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('nolim' in x) and ('1samples' in x) and ('120arcsec' in x) and ('_%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and (handpicked in x) and (maglim in x) and (('W1' in x) | ('W2' in x) | ('W3' in x) | ('W4' in x)) and (('50' in x) | ('75' in x)) and (('bpz' in x) & ('eazy' not in x)) and (('IRAC' in x) & ('noIRAC' not in x)) and (('_detir_' in x) & ('_deti_' not in x)) ] # and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE

# read the samples and classify by photoz, type, and detection
bpz_deti_irac45 = 0
#bpz_deti_noirac45 = 0
bpz_detir_irac45 = 0
#bpz_detir_noirac45 = 0
eazy_deti_irac45 = 0
#eazy_deti_noirac45 = 0
eazy_detir_irac45 = 0
#eazy_detir_noirac45 = 0
bpz_deti_irac120 = 0
#bpz_deti_noirac120 = 0
bpz_detir_irac120 = 0
#bpz_detir_noirac120 = 0
eazy_deti_irac120 = 0
#eazy_deti_noirac120 = 0
eazy_detir_irac120 = 0
#eazy_detir_noirac120 = 0
for i in range(len(lst45)):
    print lst45[i]
    if i == 0: x45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    else: x45 = np.c_[x45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)] # I checked using np.shape that if I use unpack=True I need to use np.c_
    #if ('bpz' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(bpz_deti_irac45) == int): bpz_deti_irac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('bpz' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' in lst45[i]) and (type(bpz_deti_noirac45) == int): bpz_deti_noirac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    if ('bpz' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(bpz_detir_irac45) == int): bpz_detir_irac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('bpz' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' in lst45[i]) and (type(bpz_detir_noirac45) == int): bpz_detir_noirac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('eazy' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(eazy_deti_irac45) == int): eazy_deti_irac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('eazy' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' in lst45[i]) and (type(eazy_deti_noirac45) == int): eazy_deti_noirac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('eazy' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(eazy_detir_irac45) == int): eazy_detir_irac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('eazy' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' in lst45[i]) and (type(eazy_detir_noirac45) == int): eazy_detir_noirac45 = np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)
    #if ('bpz' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(bpz_deti_irac45) != int): bpz_deti_irac45 = np.c_[bpz_deti_irac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('bpz' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' in lst45[i]) and (type(bpz_deti_noirac45) != int): bpz_deti_noirac45 = np.c_[bpz_deti_noirac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    if ('bpz' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(bpz_detir_irac45) != int): bpz_detir_irac45 = np.c_[bpz_detir_irac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('bpz' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' in lst45[i]) and (type(bpz_detir_noirac45) != int): bpz_detir_noirac45 = np.c_[bpz_detir_noirac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('eazy' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(eazy_deti_irac45) != int): eazy_deti_irac45 = np.c_[eazy_deti_irac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('eazy' in lst45[i]) and ('detir' not in lst45[i]) and ('noIRAC' in lst45[i]) and (type(eazy_deti_noirac45) != int): eazy_deti_noirac45 = np.c_[eazy_deti_noirac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('eazy' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' not in lst45[i]) and (type(eazy_detir_irac45) != int): eazy_detir_irac45 = np.c_[eazy_detir_irac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
    #if ('eazy' in lst45[i]) and ('detir' in lst45[i]) and ('noIRAC' in lst45[i]) and (type(eazy_detir_noirac45) != int): eazy_detir_noirac45 = np.c_[eazy_detir_noirac45,np.loadtxt('%s%s' %(rootin,lst45[i]), unpack=True)]
for i in range(len(lst120)):
    if i == 0: x120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    else: x120 = np.c_[x120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('bpz' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(bpz_deti_irac120) == int): bpz_deti_irac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('bpz' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' in lst120[i]) and (type(bpz_deti_noirac120) == int): bpz_deti_noirac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    if ('bpz' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(bpz_detir_irac120) == int): bpz_detir_irac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('bpz' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' in lst120[i]) and (type(bpz_detir_noirac120) == int): bpz_detir_noirac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('eazy' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(eazy_deti_irac120) == int): eazy_deti_irac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('eazy' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' in lst120[i]) and (type(eazy_deti_noirac120) == int): eazy_deti_noirac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('eazy' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(eazy_detir_irac120) == int): eazy_detir_irac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    #if ('eazy' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' in lst120[i]) and (type(eazy_detir_noirac120) == int): eazy_detir_noirac120 = np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)
    if ('bpz' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(bpz_deti_irac120) != int): bpz_deti_irac120 = np.c_[bpz_deti_irac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('bpz' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' in lst120[i]) and (type(bpz_deti_noirac120) != int): bpz_deti_noirac120 = np.c_[bpz_deti_noirac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('bpz' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(bpz_detir_irac120) != int): bpz_detir_irac120 = np.c_[bpz_detir_irac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('bpz' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' in lst120[i]) and (type(bpz_detir_noirac120) != int): bpz_detir_noirac120 = np.c_[bpz_detir_noirac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('eazy' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(eazy_deti_irac120) != int): eazy_deti_irac120 = np.c_[eazy_deti_irac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('eazy' in lst120[i]) and ('detir' not in lst120[i]) and ('noIRAC' in lst120[i]) and (type(eazy_deti_noirac120) != int): eazy_deti_noirac120 = np.c_[eazy_deti_noirac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('eazy' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' not in lst120[i]) and (type(eazy_detir_irac120) != int): eazy_detir_irac120 = np.c_[eazy_detir_irac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]
    #if ('eazy' in lst120[i]) and ('detir' in lst120[i]) and ('noIRAC' in lst120[i]) and (type(eazy_detir_noirac120) != int): eazy_detir_noirac120 = np.c_[eazy_detir_noirac120,np.loadtxt('%s%s' %(rootin,lst120[i]), unpack=True)]

def percentile(a,b,c):
    str_a = '%.2f' %a
    str_b = '%.2f' %b
    str_c = '%.2f' %c
    if str_b >= str_a: b = a - 0.01
    if str_c <= str_a: c = a + 0.01
    return b,a,c

# plot for each classification
fontlegend = 8
fontsize = 6
fontordonate = 6
fontabsciss = 6
fontlegend = 7
fontlabel = 2
linewidth = 0.5

def plot(mag,radius):
    plt.clf()
    plt.suptitle(r'50 and 75 W1-W4 %s %s arcsec %s inner %s %s %s zgap %s %s' % (lens, radius, inner, mag, mode, handpicked, zinf, zsup), fontsize=fontsize, y=0.998)
    #for i in range(18):
    for i in range(2):
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

        #z=np.linspace(0,1,1000)
        #w=np.ones(1000) * percentile(np.median(bpz_deti_irac[i]),np.percentile(bpz_deti_irac[i], 16),np.percentile(bpz_deti_irac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(bpz_deti_irac[i]),np.percentile(bpz_deti_irac[i], 16),np.percentile(bpz_deti_irac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(bpz_deti_irac[i]),np.percentile(bpz_deti_irac[i], 16),np.percentile(bpz_deti_irac[i], 84))[2]
        #plt.plot(w,z,'b-',label='bpz_deti_irac', linewidth=linewidth)
        #plt.plot(winf,z,'b--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'b--',label=None, linewidth=linewidth)
        #z=np.linspace(1,2,1000)
        #w=np.ones(1000) * percentile(np.median(bpz_deti_noirac[i]),np.percentile(bpz_deti_noirac[i], 16),np.percentile(bpz_deti_noirac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(bpz_deti_noirac[i]),np.percentile(bpz_deti_noirac[i], 16),np.percentile(bpz_deti_noirac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(bpz_deti_noirac[i]),np.percentile(bpz_deti_noirac[i], 16),np.percentile(bpz_deti_noirac[i], 84))[2]
        #plt.plot(w,z,'g-',label='bpz_deti_noirac', linewidth=linewidth)
        #plt.plot(winf,z,'g--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'g--',label=None, linewidth=linewidth)
        z=np.linspace(2,3,1000)
        w=np.ones(1000) * percentile(np.median(bpz_detir_irac[i]),np.percentile(bpz_detir_irac[i], 16),np.percentile(bpz_detir_irac[i], 84))[1]
        winf=np.ones(1000) * percentile(np.median(bpz_detir_irac[i]),np.percentile(bpz_detir_irac[i], 16),np.percentile(bpz_detir_irac[i], 84))[0]
        wsup=np.ones(1000) * percentile(np.median(bpz_detir_irac[i]),np.percentile(bpz_detir_irac[i], 16),np.percentile(bpz_detir_irac[i], 84))[2]
        plt.plot(w,z,'k-',label='bpz_detir_irac', linewidth=linewidth)
        plt.plot(winf,z,'k--',label=None, linewidth=linewidth)
        plt.plot(wsup,z,'k--',label=None, linewidth=linewidth)
        #z=np.linspace(3,4,1000)
        #w=np.ones(1000) * percentile(np.median(bpz_detir_noirac[i]),np.percentile(bpz_detir_noirac[i], 16),np.percentile(bpz_detir_noirac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(bpz_detir_noirac[i]),np.percentile(bpz_detir_noirac[i], 16),np.percentile(bpz_detir_noirac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(bpz_detir_noirac[i]),np.percentile(bpz_detir_noirac[i], 16),np.percentile(bpz_detir_noirac[i], 84))[2]
        #plt.plot(w,z,'r-',label='bpz_detir_noirac', linewidth=linewidth)
        #plt.plot(winf,z,'r--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'r--',label=None, linewidth=linewidth)
        #z=np.linspace(4,5,1000)
        #w=np.ones(1000) * percentile(np.median(eazy_deti_irac[i]),np.percentile(eazy_deti_irac[i], 16),np.percentile(eazy_deti_irac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(eazy_deti_irac[i]),np.percentile(eazy_deti_irac[i], 16),np.percentile(eazy_deti_irac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(eazy_deti_irac[i]),np.percentile(eazy_deti_irac[i], 16),np.percentile(eazy_deti_irac[i], 84))[2]
        #plt.plot(w,z,'m-',label='eazy_deti_irac', linewidth=linewidth)
        #plt.plot(winf,z,'m--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'m--',label=None, linewidth=linewidth)
        #z=np.linspace(5,6,1000)
        #w=np.ones(1000) * percentile(np.median(eazy_deti_noirac[i]),np.percentile(eazy_deti_noirac[i], 16),np.percentile(eazy_deti_noirac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(eazy_deti_noirac[i]),np.percentile(eazy_deti_noirac[i], 16),np.percentile(eazy_deti_noirac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(eazy_deti_noirac[i]),np.percentile(eazy_deti_noirac[i], 16),np.percentile(eazy_deti_noirac[i], 84))[2]
        #plt.plot(w,z,'y-',label='eazy_deti_noirac', linewidth=linewidth)
        #plt.plot(winf,z,'y--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'y--',label=None, linewidth=linewidth)
        #z=np.linspace(6,7,1000)
        #w=np.ones(1000) * percentile(np.median(eazy_detir_irac[i]),np.percentile(eazy_detir_irac[i], 16),np.percentile(eazy_detir_irac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(eazy_detir_irac[i]),np.percentile(eazy_detir_irac[i], 16),np.percentile(eazy_detir_irac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(eazy_detir_irac[i]),np.percentile(eazy_detir_irac[i], 16),np.percentile(eazy_detir_irac[i], 84))[2]
        #plt.plot(w,z,'c-',label='eazy_detir_irac', linewidth=linewidth)
        #plt.plot(winf,z,'c--',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,'c--',label=None, linewidth=linewidth)
        #z=np.linspace(7,8,1000)
        #w=np.ones(1000) * percentile(np.median(eazy_detir_noirac[i]),np.percentile(eazy_detir_noirac[i], 16),np.percentile(eazy_detir_noirac[i], 84))[1]
        #winf=np.ones(1000) * percentile(np.median(eazy_detir_noirac[i]),np.percentile(eazy_detir_noirac[i], 16),np.percentile(eazy_detir_noirac[i], 84))[0]
        #wsup=np.ones(1000) * percentile(np.median(eazy_detir_noirac[i]),np.percentile(eazy_detir_noirac[i], 16),np.percentile(eazy_detir_noirac[i], 84))[2]
        #plt.plot(w,z,linestyle='-',color='grey',label='eazy_detir_noirac', linewidth=linewidth)
        #plt.plot(winf,z,linestyle='--',color='grey',label=None, linewidth=linewidth)
        #plt.plot(wsup,z,linestyle='--',color='grey',label=None, linewidth=linewidth)
        #plt.xlim(0, 2.5)
        #plt.ylim(0, vertlimit)

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
        plt.tick_params(axis='x', labelsize=2)
        plt.tick_params(axis='y', labelsize=2)
        plt.setp(plt.xticks()[1], rotation=90)
        subplot = i+1
    plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.6)
    plt.subplot(5,4,5)
    plt.legend(bbox_to_anchor=(5, -5), loc='lower right', borderaxespad=0., fontsize=fontlegend)
    plt.savefig('%s/weightedcountshist_%sarcsec_%sinner_%s_%s%s_zgap%s_%s.png' % (rootout, radius, inner, mag, mode, handpicked, zinf, zsup), dpi=500)

#bpz_deti_irac = bpz_deti_irac45
#bpz_deti_noirac = bpz_deti_noirac45
bpz_detir_irac = bpz_detir_irac45
#bpz_detir_noirac = bpz_detir_noirac45
#eazy_deti_irac = eazy_deti_irac45
#eazy_deti_noirac = eazy_deti_noirac45
#eazy_detir_irac = eazy_detir_irac45
#eazy_detir_noirac = eazy_detir_noirac45
plot(maglim,'45')
#bpz_deti_irac = bpz_deti_irac120
#bpz_deti_noirac = bpz_deti_noirac120
bpz_detir_irac = bpz_detir_irac120
#bpz_detir_noirac = bpz_detir_noirac120
#eazy_deti_irac = eazy_deti_irac120
#eazy_deti_noirac = eazy_deti_noirac120
#eazy_detir_irac = eazy_detir_irac120
#eazy_detir_noirac = eazy_detir_noirac120
plot(maglim,'120')

# output the final summary file
if local == 'global':
    f = open('%s/weightedcounts_%s_%s_%s_%sinner%s_zgap%s_%s_%s.cat' %(rootout,lens,mode,maglim,inner,handpicked,zinf,zsup,local),'w')
    str = '# weight      45med    45inf  45sup 120med 120inf 120sup \n'
    str += 'gal           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[1],percentile(np.median(x45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[0],percentile(np.median(x45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[2],percentile(np.median(x120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[1],percentile(np.median(x120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[0],percentile(np.median(x120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[2])
    #str += 'z             %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[1],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[0],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[2],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[1],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[0],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[2])
    #str += 'mass          %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[1],percentile(np.median(x45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[0],percentile(np.median(x45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[2],percentile(np.median(x120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[1],percentile(np.median(x120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[0],percentile(np.median(x120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[2])
    #str += 'mass2         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[1],percentile(np.median(x45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[0],percentile(np.median(x45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[2],percentile(np.median(x120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[1],percentile(np.median(x120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[0],percentile(np.median(x120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[2])
    #str += 'mass3         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[1],percentile(np.median(x45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[0],percentile(np.median(x45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[2],percentile(np.median(x120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[1],percentile(np.median(x120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[0],percentile(np.median(x120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[2])
    #str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[1],percentile(np.median(x45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[0],percentile(np.median(x45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[2],percentile(np.median(x120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[1],percentile(np.median(x120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[0],percentile(np.median(x120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[2])
    str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[1],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[0],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[2],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[1],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[0],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[2])
    #str += 'zoverr        %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[1],percentile(np.median(x45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[0],percentile(np.median(x45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[2],percentile(np.median(x120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[1],percentile(np.median(x120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[0],percentile(np.median(x120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[2])
    #str += 'massoverr     %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[1],percentile(np.median(x45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[0],percentile(np.median(x45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[2],percentile(np.median(x120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[1],percentile(np.median(x120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[0],percentile(np.median(x120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[2])
    #str += 'mass2overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[1],percentile(np.median(x45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[0],percentile(np.median(x45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[2],percentile(np.median(x120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[1],percentile(np.median(x120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[0],percentile(np.median(x120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[2])
    #str += 'mass3overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[1],percentile(np.median(x45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[0],percentile(np.median(x45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[2],percentile(np.median(x120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[1],percentile(np.median(x120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[0],percentile(np.median(x120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[2])
    #str += 'mass2rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[1],percentile(np.median(x45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[0],percentile(np.median(x45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[2],percentile(np.median(x120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[1],percentile(np.median(x120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[0],percentile(np.median(x120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[2])
    #str += 'mass3rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[1],percentile(np.median(x45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[0],percentile(np.median(x45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[2],percentile(np.median(x120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[1],percentile(np.median(x120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[0],percentile(np.median(x120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[2])
    #str += 'mass2overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[1],percentile(np.median(x45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[0],percentile(np.median(x45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[2],percentile(np.median(x120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[1],percentile(np.median(x120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[0],percentile(np.median(x120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[2])
    #str += 'mass3overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[1],percentile(np.median(x45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[0],percentile(np.median(x45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[2],percentile(np.median(x120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[1],percentile(np.median(x120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[0],percentile(np.median(x120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[2])
    #str += 'flexion       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[1],percentile(np.median(x45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[0],percentile(np.median(x45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[2],percentile(np.median(x120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[1],percentile(np.median(x120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[0],percentile(np.median(x120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[2])
    #str += 'tidal         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[1],percentile(np.median(x45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[0],percentile(np.median(x45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[2],percentile(np.median(x120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[1],percentile(np.median(x120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[0],percentile(np.median(x120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[2])
    #str += 'SIS           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[1],percentile(np.median(x45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[0],percentile(np.median(x45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[2],percentile(np.median(x120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[1],percentile(np.median(x120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[0],percentile(np.median(x120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[2])
    #str += 'SIShalo       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[1],percentile(np.median(x45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[0],percentile(np.median(x45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[2],percentile(np.median(x120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[1],percentile(np.median(x120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[0],percentile(np.median(x120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[2])
    f.write(str)
    f.close()

else:
    print "fiducial:"
    fiducial45 = [x for x in os.listdir(rootin) if ('%s_weightedcountshist_45arcsec_%sinner_%s_%s_bpz_detir_IRAC%s_zgap%s_%s_21samples' % (lens,inner,maglim,mode,handpicked,zinf,zsup) in x) and ('.lst' in x)] # this will contain W1-W4 and 50/75; CHOOSE WHAT YOU WANT HERE
    fiducial120 = [x for x in os.listdir(rootin) if ('%s_weightedcountshist_120arcsec_%sinner_%s_%s_bpz_detir_IRAC%s_zgap%s_%s_21samples' % (lens,inner,maglim,mode,handpicked,zinf,zsup) in x) and ('.lst' in x)]
    for i in range(len(fiducial45)):
        print fiducial45[i]
        if i == 0: f45 = np.loadtxt('%s%s' %(rootin,fiducial45[0]), unpack=True)
        else: f45 = np.c_[f45,np.loadtxt('%s%s' %(rootin,fiducial45[i]), unpack=True)]
    for i in range(len(fiducial120)):
        if i == 0: f120 = np.loadtxt('%s%s' %(rootin,fiducial120[0]), unpack=True)
        else: f120 = np.c_[f120,np.loadtxt('%s%s' %(rootin,fiducial120[i]), unpack=True)]

    f = open('%s/weightedcounts_%s_%s_%s_%sinner%s_zgap%s_%s_%s.cat' %(rootout,lens,mode,maglim,inner,handpicked,zinf,zsup,local),'w')
    str = '# weight      45med    45inf  45sup 120med 120inf 120sup \n'
    str += 'gal           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[1],percentile(np.median(f45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[0],percentile(np.median(f45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84))[2],percentile(np.median(f120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[1],percentile(np.median(f120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[0],percentile(np.median(f120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))[2])
    #str += 'z             %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[1],percentile(np.median(f45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[0],percentile(np.median(f45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[2],percentile(np.median(f120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[1],percentile(np.median(f120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[0],percentile(np.median(f120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[2])
    #str += 'mass          %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[1],percentile(np.median(f45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[0],percentile(np.median(f45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84))[2],percentile(np.median(f120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[1],percentile(np.median(f120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[0],percentile(np.median(f120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))[2])
    #str += 'mass2         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[1],percentile(np.median(f45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[0],percentile(np.median(f45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84))[2],percentile(np.median(f120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[1],percentile(np.median(f120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[0],percentile(np.median(f120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))[2])
    #str += 'mass3         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[1],percentile(np.median(f45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[0],percentile(np.median(f45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84))[2],percentile(np.median(f120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[1],percentile(np.median(f120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[0],percentile(np.median(f120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))[2])
    #str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[1],percentile(np.median(f45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[0],percentile(np.median(f45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84))[2],percentile(np.median(f120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[1],percentile(np.median(f120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[0],percentile(np.median(f120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))[2])
    str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[1],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[0],percentile(np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84))[2],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[1],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[0],percentile(np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))[2])
    #str += 'zoverr        %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[1],percentile(np.median(f45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[0],percentile(np.median(f45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84))[2],percentile(np.median(f120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[1],percentile(np.median(f120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[0],percentile(np.median(f120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))[2])
    #str += 'massoverr     %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[1],percentile(np.median(f45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[0],percentile(np.median(f45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84))[2],percentile(np.median(f120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[1],percentile(np.median(f120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[0],percentile(np.median(f120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))[2])
    #str += 'mass2overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[1],percentile(np.median(f45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[0],percentile(np.median(f45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84))[2],percentile(np.median(f120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[1],percentile(np.median(f120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[0],percentile(np.median(f120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))[2])
    #str += 'mass3overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[1],percentile(np.median(f45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[0],percentile(np.median(f45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84))[2],percentile(np.median(f120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[1],percentile(np.median(f120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[0],percentile(np.median(f120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))[2])
    #str += 'mass2rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[1],percentile(np.median(f45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[0],percentile(np.median(f45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84))[2],percentile(np.median(f120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[1],percentile(np.median(f120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[0],percentile(np.median(f120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))[2])
    #str += 'mass3rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[1],percentile(np.median(f45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[0],percentile(np.median(f45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84))[2],percentile(np.median(f120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[1],percentile(np.median(f120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[0],percentile(np.median(f120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))[2])
    #str += 'mass2overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[1],percentile(np.median(f45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[0],percentile(np.median(f45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84))[2],percentile(np.median(f120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[1],percentile(np.median(f120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[0],percentile(np.median(f120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))[2])
    #str += 'mass3overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[1],percentile(np.median(f45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[0],percentile(np.median(f45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84))[2],percentile(np.median(f120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[1],percentile(np.median(f120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[0],percentile(np.median(f120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))[2])
    #str += 'flexion       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[1],percentile(np.median(f45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[0],percentile(np.median(f45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84))[2],percentile(np.median(f120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[1],percentile(np.median(f120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[0],percentile(np.median(f120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))[2])
    #str += 'tidal         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[1],percentile(np.median(f45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[0],percentile(np.median(f45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84))[2],percentile(np.median(f120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[1],percentile(np.median(f120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[0],percentile(np.median(f120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))[2])
    #str += 'SIS           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[1],percentile(np.median(f45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[0],percentile(np.median(f45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84))[2],percentile(np.median(f120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[1],percentile(np.median(f120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[0],percentile(np.median(f120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))[2])
    #str += 'SIShalo       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (percentile(np.median(f45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[1],percentile(np.median(f45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[0],percentile(np.median(f45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84))[2],percentile(np.median(f120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[1],percentile(np.median(f120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[0],percentile(np.median(f120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))[2])
    f.write(str)
    f.close()
