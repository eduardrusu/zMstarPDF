# CE Rusu, Feb 13 2018
# The code uses the weighted count ratios derived by weightinguniversal_overlap_sampling_nobeta_WFI2033rethought.py to produce histograms and compute the 16th, 50th and 84th percentiles, using the 10 samples
# What is actually plotted is the distribution for the final sample (of the various samples), but the plotted statistics refer to all samples. If I want to plot the first (fiducial) sample, I need to run the code with samples=1
# run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_histograms_samplesFITSinput.py WFI2033 45 5 22.5 meds bpz detir IRAC -1.0 -1.0 100 removelensgrouphandpicked testduplicatesamples/testothersamples
# After running this code, you need to combine the results into a final text file with the final distributions and widths, using weightinguniversal_histograms_finalcombine.py
# If desired, run weightinguniversal_histograms_samples_publicationqualitynotext.py to produce a publication quality plot

import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt
from astropy.io import fits

lens = str(sys.argv[1])
radius = str(sys.argv[2])
inner = str(sys.argv[3])
mag = str(sys.argv[4])
mode = str(sys.argv[5])
photz = str(sys.argv[6])
detect = str(sys.argv[7])
irac = str(sys.argv[8])
zinf = str(sys.argv[9])
zsup = str(sys.argv[10])
bin = int(str(sys.argv[11]))
try: handpicked = '_'+str(sys.argv[12])
except: handpicked = ''
try: specialtest = '_'+str(sys.argv[13])
except: specialtest = ''

plt.clf()

fontlegend = 8
fontsize = 8
fontordonate = 4
fontabsciss = 8
fontlabel = 2
pltrange = 3
samples = 4
#limit = 10**30
limit = np.infty
root = "/Volumes/LaCieSubaru/weightedcounts/%s/" % lens

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
    print '%s/%s' %(nr+1,samples)
    #lstW1_75 = [x for x in os.listdir(root) if ('W1' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.fits' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)] # select from the files in the root directory
    lstW1_75 = [x for x in os.listdir(root) if ('W1' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked) in x) and ('_%s.fits' % nr in x)]
    print lstW1_75
    #lstW2_75 = [x for x in os.listdir(root) if ('W2' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.fits' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW2_75 = [x for x in os.listdir(root) if ('W2' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked) in x) and ('_%s.fits' % nr in x)]
    #lstW3_75 = [x for x in os.listdir(root) if ('W3' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.fits' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW3_75 = [x for x in os.listdir(root) if ('W3' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked) in x) and ('_%s.fits' % nr in x)]
    #lstW4_75 = [x for x in os.listdir(root) if ('W4' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.fits' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW4_75 = [x for x in os.listdir(root) if ('W4' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked) in x) and ('_%s.fits' % nr in x)]
    print lstW1_75

    print "W1..."
#    for i in range(len(lstW1_50)):
#        hdu = fits.open(root+lstW1_50[i]); data = hdu[1].data
#        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
#        if i == 0:
#            q_W1_50read = dataread
#        else: q_W1_50read = np.r_['1',q_W1_50read,dataread]
#        #print np.shape(q_W1_50read)
#        hdu.close()
#        #print np.shape(q_W1_50read)
    for i in range(len(lstW1_75)):
        hdu = fits.open(root+lstW1_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W1_75read = dataread
        else: q_W1_75read = np.r_['1',q_W1_75read,dataread]
        hdu.close()

    print "W2..."
#    for i in range(len(lstW2_50)):
#        hdu = fits.open(root+lstW2_50[i]); data = hdu[1].data
#        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
#        if i == 0:
#            q_W2_50read = dataread
#        else: q_W2_50read = np.r_['1',q_W2_50read,dataread]
#        hdu.close()
#        #print np.shape(q_W2_50read)
    for i in range(len(lstW2_75)):
        hdu = fits.open(root+lstW2_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W2_75read = dataread
        else: q_W2_75read = np.r_['1',q_W2_75read,dataread]
        hdu.close()

    print "W3..."
#    for i in range(len(lstW3_50)):
#        hdu = fits.open(root+lstW3_50[i]); data = hdu[1].data
#        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
#        if i == 0:
#            q_W3_50read = dataread
#        else: q_W3_50read = np.r_['1',q_W3_50read,dataread]
#        hdu.close()
#        #print np.shape(q_W3_50read)
    for i in range(len(lstW3_75)):
        hdu = fits.open(root+lstW3_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W3_75read = dataread
        else: q_W3_75read = np.r_['1',q_W3_75read,dataread]
        hdu.close()

    print "W4..."
#    for i in range(len(lstW4_50)):
#        hdu = fits.open(root+lstW4_50[i]); data = hdu[1].data
#        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
#        if i == 0:
#            q_W4_50read = dataread
#        else: q_W4_50read = np.r_['1',q_W4_50read,dataread]
#        hdu.close()
#        #print np.shape(q_W4_50read)
    for i in range(len(lstW4_75)):
        hdu = fits.open(root+lstW4_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W4_75read = dataread
        else: q_W4_75read = np.r_['1',q_W4_75read,dataread]
        hdu.close()

    for j in range(18):
        #q_W1_50 = q_W1_50read[j][q_W1_50read[j] < limit]
        #if mode == "sum": q_W1_50 = abs(q_W1_50)
        q_W1_75 = q_W1_75read[j][q_W1_75read[j] < limit]
        if mode == "sum": q_W1_75 = abs(q_W1_75)
        #q_W2_50 = q_W2_50read[j][q_W2_50read[j] < limit]
        #if mode == "sum": q_W2_50 = abs(q_W2_50)
        q_W2_75 = q_W2_75read[j][q_W2_75read[j] < limit]
        if mode == "sum": q_W2_75 = abs(q_W2_75)
        #q_W3_50 = q_W3_50read[j][q_W3_50read[j] < limit]
        #if mode == "sum": q_W3_50 = abs(q_W3_50)
        q_W3_75 = q_W3_75read[j][q_W3_75read[j] < limit]
        if mode == "sum": q_W3_75 = abs(q_W3_75)
        #q_W4_50 = q_W4_50read[j][q_W4_50read[j] < limit]
        #if mode == "sum": q_W4_50 = abs(q_W4_50)
        q_W4_75 = q_W4_75read[j][q_W4_75read[j] < limit]
        if mode == "sum": q_W4_75 = abs(q_W4_75)

        if mode == "sum":
            #medsum50W1[j][nr] = np.average(q_W1_50)
            medsum75W1[j][nr] = np.average(q_W1_75)
            #medsum50W2[j][nr] = np.average(q_W2_50)
            medsum75W2[j][nr] = np.average(q_W2_75)
            #medsum50W3[j][nr] = np.average(q_W3_50)
            medsum75W3[j][nr] = np.average(q_W3_75)
            #medsum50W4[j][nr] = np.average(q_W4_50)
            medsum75W4[j][nr] = np.average(q_W4_75)
        if mode == "meds":
            #medsum50W1[j][nr] = np.median(q_W1_50)
            medsum75W1[j][nr] = np.median(q_W1_75)
            #medsum50W2[j][nr] = np.median(q_W2_50)
            medsum75W2[j][nr] = np.median(q_W2_75)
            #medsum50W3[j][nr] = np.median(q_W3_50)
            medsum75W3[j][nr] = np.median(q_W3_75)
            #medsum50W4[j][nr] = np.median(q_W4_50)
            medsum75W4[j][nr] = np.median(q_W4_75)

#std50W1_inf = np.zeros(18)
#std50W1_sup = np.zeros(18)
#std50W2_inf = np.zeros(18)
#std50W2_sup = np.zeros(18)
#std50W3_inf = np.zeros(18)
#std50W3_sup = np.zeros(18)
#std50W4_inf = np.zeros(18)
#std50W4_sup = np.zeros(18)
std75W1_inf = np.zeros(18)
std75W1_sup = np.zeros(18)
std75W2_inf = np.zeros(18)
std75W2_sup = np.zeros(18)
std75W3_inf = np.zeros(18)
std75W3_sup = np.zeros(18)
std75W4_inf = np.zeros(18)
std75W4_sup = np.zeros(18)
std_inf = np.zeros(18)
std_sup = np.zeros(18)

for i in range(18):
    #std50W1_inf[i] = np.percentile(medsum50W1[i], 16)
    #std50W1_sup[i] = np.percentile(medsum50W1[i], 84)
    #std50W2_inf[i] = np.percentile(medsum50W2[i], 16)
    #std50W2_sup[i] = np.percentile(medsum50W2[i], 84)
    #std50W3_inf[i] = np.percentile(medsum50W3[i], 16)
    #std50W3_sup[i] = np.percentile(medsum50W3[i], 84)
    #std50W4_inf[i] = np.percentile(medsum50W4[i], 16)
    #std50W4_sup[i] = np.percentile(medsum50W4[i], 84)
    std75W1_inf[i] = np.percentile(medsum75W1[i], 16)
    std75W1_sup[i] = np.percentile(medsum75W1[i], 84)
    std75W2_inf[i] = np.percentile(medsum75W2[i], 16)
    std75W2_sup[i] = np.percentile(medsum75W2[i], 84)
    std75W3_inf[i] = np.percentile(medsum75W3[i], 16)
    std75W3_sup[i] = np.percentile(medsum75W3[i], 84)
    std75W4_inf[i] = np.percentile(medsum75W4[i], 16)
    std75W4_sup[i] = np.percentile(medsum75W4[i], 84)
    #std_inf[i] = np.percentile([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i],medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]], 16)
    #std_sup[i] = np.percentile([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i],medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]], 84)
    std_inf[i] = np.percentile([medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]], 16)
    std_sup[i] = np.percentile([medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]], 84)

print "Plotting..."

plt.suptitle(r'%s weighted counts histogram W1-W4 %s arcsec %s inner %s %s %s %s %s %s zgap %s %s' % (lens, radius, inner, mag, mode, photz, irac, detect, handpicked, zinf, zsup), fontsize=fontsize, y=0.998)

for i in range(18):
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

    #q_W1_50 = q_W1_50read[i][q_W1_50read[i] < limit]
    #if mode == "sum": q_W1_50 = abs(q_W1_50) # fix the negative halo convergence
    q_W1_75 = q_W1_75read[i][q_W1_75read[i] < limit]
    if mode == "sum": q_W1_75 = abs(q_W1_75)
    #q_W2_50 = q_W2_50read[i][q_W2_50read[i] < limit]
    #if mode == "sum": q_W2_50 = abs(q_W2_50)
    q_W2_75 = q_W2_75read[i][q_W2_75read[i] < limit]
    if mode == "sum": q_W2_75 = abs(q_W2_75)
    #q_W3_50 = q_W3_50read[i][q_W3_50read[i] < limit]
    #if mode == "sum": q_W3_50 = abs(q_W3_50)
    q_W3_75 = q_W3_75read[i][q_W3_75read[i] < limit]
    if mode == "sum": q_W3_75 = abs(q_W3_75)
    #q_W4_50 = q_W4_50read[i][q_W4_50read[i] < limit]
    #if mode == "sum": q_W4_50 = abs(q_W4_50)
    q_W4_75 = q_W4_75read[i][q_W4_75read[i] < limit]
    if mode == "sum": q_W4_75 = abs(q_W4_75)

    #plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    #plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    #plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    #plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange], linestyle='dotted')
    plt.hist(q_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])

#    if mode == "sum":
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W1[i][0],np.average(medsum50W1[i]),std50W1_inf[i],std50W1_sup[i],medsum75W1[i][0],np.average(medsum75W1[i]),std75W1_inf[i],std75W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
#        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W2[i][0],np.average(medsum50W2[i]),std50W2_inf[i],std50W2_sup[i],medsum75W2[i][0],np.average(medsum75W2[i]),std75W2_inf[i],std75W2_sup[i])
#        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W3[i][0],np.average(medsum50W3[i]),std50W3_inf[i],std50W3_sup[i],medsum75W3[i][0],np.average(medsum75W3[i]),std75W3_inf[i],std75W3_sup[i])
#        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W4[i][0],np.average(medsum50W4[i]),std50W4_inf[i],std50W4_sup[i],medsum75W4[i][0],np.average(medsum75W4[i]),std75W4_inf[i],std75W4_sup[i])
#        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
#        s = "W1-4 init %.3f all %.3f (%.3f %.3f)" % (np.average([medsum50W1[i][0],medsum50W2[i][0],medsum50W3[i][0],medsum50W4[i][0],medsum75W1[i][0],medsum75W2[i][0],medsum75W3[i][0],medsum75W4[i][0]]),np.average([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i],medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]]),std_inf[i],std_sup[i])
#        ax.text(0.02, 0.1, s, fontsize=fontlabel+1, color='k',transform=ax.transAxes)
#    if mode == "meds":
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W1[i][0],np.median(medsum50W1[i]),std50W1_inf[i],std50W1_sup[i],medsum75W1[i][0],np.median(medsum75W1[i]),std75W1_inf[i],std75W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
#        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W2[i][0],np.median(medsum50W2[i]),std50W2_inf[i],std50W2_sup[i],medsum75W2[i][0],np.median(medsum75W2[i]),std75W2_inf[i],std75W2_sup[i])
#        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W3[i][0],np.median(medsum50W3[i]),std50W3_inf[i],std50W3_sup[i],medsum75W3[i][0],np.median(medsum75W3[i]),std75W3_inf[i],std75W3_sup[i])
#        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
#        s = "50: init %.3f all %.3f (%.3f %.3f); 75: %.3f %.3f (%.3f %.3f)" % (medsum50W4[i][0],np.median(medsum50W4[i]),std50W4_inf[i],std50W4_sup[i],medsum75W4[i][0],np.median(medsum75W4[i]),std75W4_inf[i],std75W4_sup[i])
#        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
#        s = "W1-4 init %.3f all %.3f (%.3f %.3f)" % (np.median([medsum50W1[i][0],medsum50W2[i][0],medsum50W3[i][0],medsum50W4[i][0],medsum75W1[i][0],medsum75W2[i][0],medsum75W3[i][0],medsum75W4[i][0]]),np.median([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i],medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]]),std_inf[i],std_sup[i])
#        ax.text(0.02, 0.1, s, fontsize=fontlabel+1, color='k',transform=ax.transAxes)

    if mode == "sum":
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W1[i][0],np.average(medsum75W1[i]),std75W1_inf[i],std75W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W2[i][0],np.average(medsum75W2[i]),std75W2_inf[i],std75W2_sup[i])
        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W3[i][0],np.average(medsum75W3[i]),std75W3_inf[i],std75W3_sup[i])
        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W4[i][0],np.average(medsum75W4[i]),std75W4_inf[i],std75W4_sup[i])
        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 init %.3f all %.3f (%.3f %.3f)" % (np.average([medsum75W1[i][0],medsum75W2[i][0],medsum75W3[i][0],medsum75W4[i][0]]),np.average([medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]]),std_inf[i],std_sup[i])
        ax.text(0.02, 0.1, s, fontsize=fontlabel+1, color='k',transform=ax.transAxes)
    if mode == "meds":
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W1[i][0],np.median(medsum75W1[i]),std75W1_inf[i],std75W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W2[i][0],np.median(medsum75W2[i]),std75W2_inf[i],std75W2_sup[i])
        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W3[i][0],np.median(medsum75W3[i]),std75W3_inf[i],std75W3_sup[i])
        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "75: %.3f %.3f (%.3f %.3f)" % (medsum75W4[i][0],np.median(medsum75W4[i]),std75W4_inf[i],std75W4_sup[i])
        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 init %.3f all %.3f (%.3f %.3f)" % (np.median([medsum75W1[i][0],medsum75W2[i][0],medsum75W3[i][0],medsum75W4[i][0]]),np.median([medsum75W1[i],medsum75W2[i],medsum75W3[i],medsum75W4[i]]),std_inf[i],std_sup[i])
        ax.text(0.02, 0.1, s, fontsize=fontlabel+1, color='k',transform=ax.transAxes)
    print i,s
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
    #print "finished subplot %d/18; fraction of points inside the < %s cut: \n W1_50 %.3f\n W2_50 %.3f\n W3_50 %.3f\n W4_50 %.3f" % (subplot, limit, float(q_W1_50.size)/q_W1_50read[0].size, float(q_W2_50.size)/q_W2_50read[0].size, float(q_W3_50.size)/q_W3_50read[0].size, float(q_W4_50.size)/q_W4_50read[0].size)

#np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW1_50%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum50W1.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW1_75%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum75W1.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
#np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW2_50%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum50W2.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW2_75%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum75W2.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
#np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW3_50%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum50W3.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW3_75%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum75W3.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
#np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW4_50%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum50W4.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
np.savetxt('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamplesW4_75%snolim.lst' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), medsum75W4.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')

plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.6)
plt.subplot(5,4,5)
plt.legend(bbox_to_anchor=(5, -5), loc='lower right', borderaxespad=0., fontsize=10)
plt.savefig('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_%ssamples%s.png' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup, samples, specialtest), dpi=500)

# compute the number of fields used
total50 = 0
total75 = 0
good50 = 0
good75 = 0
lst = [x for x in os.listdir(root) if ('_wghtratios_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_zgap%s_%s%s_count%s.cat' %(radius,inner,lens,detect,irac,mode,zinf,zsup,handpicked,specialtest) in x)]
for i in range(len(lst)):
    str = open('%s/%s' %(root,lst[i]),'r').read()
    str = [x.strip() for x in str.split(" ")]
    str_total50 = int(str[9])
    str_total75 = int(str[7][:-1])
    str_good50 = int(str[8])
    str_good75 = int(str[6])
    if i == 0:
        total50 = str_total50
        total75 = str_total75
        good50 = str_good50
        good75 = str_good75
    else:
        total50 += str_total50
        total75 += str_total75
        good50 += str_good50
        good75 += str_good75

print '50%: ', good50, '/', total50, ';', '75%: ', good75, '/', total75
print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
