# The code uses the weighted count ratios derived by weightinguniversal_overlap_sampling_nobeta_WFI2033rethought.py to produce histograms and compute the 16th, 50th and 84th percentiles. It does this by creating an output file where it samples 1000 times from the distributions of the averages/medians. For small widths of the averages/medians distributions, where they are approximately gaussian, use weightinguniversal_histograms_samples_WFI2033.py instead
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

def sample(median,stdlow,stdhigh): # samples from different standard deviation gaussians on each side of the mean
    rand = np.random.uniform(0,1,250)
    result = np.zeros(250)
    result[rand < 0.5] = median - np.abs(np.random.normal(0, np.abs(stdlow), len(rand[rand < 0.5])))
    result[rand > 0.5] = median + np.abs(np.random.normal(0, stdhigh, len(rand[rand > 0.5])))
    return result

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

samplingW1 = np.zeros((18,250))
samplingW2 = np.zeros((18,250))
samplingW3 = np.zeros((18,250))
samplingW4 = np.zeros((18,250))
sampling = np.zeros((18,1000))

for nr in range(samples):
    print '%s...' %nr
    lstW1_50 = [x for x in os.listdir(root) if ('W1' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)] # select from the files in the root directory
    lstW2_50 = [x for x in os.listdir(root) if ('W2' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW3_50 = [x for x in os.listdir(root) if ('W3' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]
    lstW4_50 = [x for x in os.listdir(root) if ('W4' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_%s.lst' %(radius,lens,detect,irac,mode,inner,str(nr)) in x)]

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
    print "W2..."
    for i in range(len(lstW2_50)):
        if i == 0:
            q_W2_50read = np.loadtxt(root+lstW2_50[i], usecols=cols, unpack=True)
        else:
            q_W2_50read = np.r_['1',q_W2_50read,np.loadtxt(root+ lstW2_50[i], usecols=cols, unpack=True)]
    print "W3..."
    for i in range(len(lstW3_50)):
        if i == 0:
            q_W3_50read = np.loadtxt(root+lstW3_50[i], usecols=cols, unpack=True)
        else:
            q_W3_50read = np.r_['1',q_W3_50read,np.loadtxt(root+ lstW3_50[i], usecols=cols, unpack=True)]
    print "W4..."
    for i in range(len(lstW4_50)):
        if i == 0:
            q_W4_50read = np.loadtxt(root+lstW4_50[i], usecols=cols, unpack=True)
        else:
            q_W4_50read = np.r_['1',q_W4_50read,np.loadtxt(root+ lstW4_50[i], usecols=cols, unpack=True)]
    
    for j in range(18):
        q_W1_50 = q_W1_50read[j][q_W1_50read[j] < limit]
        if mode == "sum": q_W1_50 = abs(q_W1_50) # fix the negative halo convergence
        q_W2_50 = q_W2_50read[j][q_W2_50read[j] < limit]
        if mode == "sum": q_W2_50 = abs(q_W2_50)
        q_W3_50 = q_W3_50read[j][q_W3_50read[j] < limit]
        if mode == "sum": q_W3_50 = abs(q_W3_50)
        q_W4_50 = q_W4_50read[j][q_W4_50read[j] < limit]
        if mode == "sum": q_W4_50 = abs(q_W4_50)

        if mode == "sum":
            medsum50W1[j][nr] = np.average(q_W1_50)
            medsum50W2[j][nr] = np.average(q_W2_50)
            medsum50W3[j][nr] = np.average(q_W3_50)
            medsum50W4[j][nr] = np.average(q_W4_50)
        if mode == "meds":
            medsum50W1[j][nr] = np.median(q_W1_50)
            medsum50W2[j][nr] = np.median(q_W2_50)
            medsum50W3[j][nr] = np.median(q_W3_50)
            medsum50W4[j][nr] = np.median(q_W4_50)
                
std50W1_inf = np.zeros(18)
std50W1_sup = np.zeros(18)
std50W2_inf = np.zeros(18)
std50W2_sup = np.zeros(18)
std50W3_inf = np.zeros(18)
std50W3_sup = np.zeros(18)
std50W4_inf = np.zeros(18)
std50W4_sup = np.zeros(18)

for i in range(18):
    std50W1_inf[i] = np.percentile(medsum50W1[i], 16)
    std50W1_sup[i] = np.percentile(medsum50W1[i], 84)
    std50W2_inf[i] = np.percentile(medsum50W2[i], 16)
    std50W2_sup[i] = np.percentile(medsum50W2[i], 84)
    std50W3_inf[i] = np.percentile(medsum50W3[i], 16)
    std50W3_sup[i] = np.percentile(medsum50W3[i], 84)
    std50W4_inf[i] = np.percentile(medsum50W4[i], 16)
    std50W4_sup[i] = np.percentile(medsum50W4[i], 84)

print "Plotting..."

plt.suptitle(r'%s weighted counts histogram W1-W4 %s %s %s %s %s arcsec %s %s' % (lens, radius, mag, mode, photz, irac, inner, detect), fontsize=10, y=0.998)

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
    
    q_W1_50 = q_W1_50read[i][q_W1_50read[i] < limit] # this includes all sampling
    if mode == "sum": q_W1_50 = abs(q_W1_50)
    q_W2_50 = q_W2_50read[i][q_W2_50read[i] < limit]
    if mode == "sum": q_W2_50 = abs(q_W2_50)
    q_W3_50 = q_W3_50read[i][q_W3_50read[i] < limit]
    if mode == "sum": q_W3_50 = abs(q_W3_50)
    q_W4_50 = q_W4_50read[i][q_W4_50read[i] < limit]
    if mode == "sum": q_W4_50 = abs(q_W4_50)
    
    plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])

    if mode == "sum":
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W1[i][0],np.average(medsum50W1[i]),std50W1_inf[i],std50W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W2[i][0],np.average(medsum50W2[i]),std50W2_inf[i],std50W2_sup[i])
        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W3[i][0],np.average(medsum50W3[i]),std50W3_inf[i],std50W3_sup[i])
        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W4[i][0],np.average(medsum50W4[i]),std50W4_inf[i],std50W4_sup[i])
        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 init %.3f all %.3f 16th=%.3f 84th=%.3f" % (np.average([medsum50W1[i][0],medsum50W2[i][0],medsum50W3[i][0],medsum50W4[i][0]]),np.average([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i]]),np.median([std50W1_inf[i],std50W2_inf[i],std50W3_inf[i],std50W4_inf[i]]),np.median([std50W1_sup[i],std50W2_sup[i],std50W3_sup[i],std50W4_sup[i]]))
        ax.text(0.02, 0.1, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        if np.average(medsum50W1[i]) == std50W1_inf[i]: std50W1_inf[i] = np.average(medsum50W1[i]) - 0.01
        if np.average(medsum50W2[i]) == std50W2_inf[i]: std50W2_inf[i] = np.average(medsum50W2[i]) - 0.01
        if np.average(medsum50W3[i]) == std50W3_inf[i]: std50W3_inf[i] = np.average(medsum50W3[i]) - 0.01
        if np.average(medsum50W4[i]) == std50W4_inf[i]: std50W4_inf[i] = np.average(medsum50W4[i]) - 0.01
        if np.average(medsum50W1[i]) == std50W1_sup[i]: std50W1_sup[i] = np.average(medsum50W1[i]) + 0.01
        if np.average(medsum50W2[i]) == std50W2_sup[i]: std50W2_sup[i] = np.average(medsum50W2[i]) + 0.01
        if np.average(medsum50W3[i]) == std50W3_sup[i]: std50W3_sup[i] = np.average(medsum50W3[i]) + 0.01
        if np.average(medsum50W4[i]) == std50W4_sup[i]: std50W4_sup[i] = np.average(medsum50W4[i]) + 0.01
        samplingW1[i] = sample(np.average(medsum50W1[i]), np.average(medsum50W1[i]) - std50W1_inf[i], std50W1_sup[i] - np.average(medsum50W1[i]))
        samplingW2[i] = sample(np.average(medsum50W2[i]), np.average(medsum50W2[i]) - std50W2_inf[i], std50W2_sup[i] - np.average(medsum50W2[i]))
        samplingW3[i] = sample(np.average(medsum50W3[i]), np.average(medsum50W3[i]) - std50W3_inf[i], std50W3_sup[i] - np.average(medsum50W3[i]))
        samplingW4[i] = sample(np.average(medsum50W4[i]), np.average(medsum50W4[i]) - std50W4_inf[i], std50W4_sup[i] - np.average(medsum50W4[i]))
        sampling[i] = np.r_[samplingW1[i],samplingW2[i],samplingW3[i],samplingW4[i]] # takes 250 samples of the sum/median from W1-4 and combines them
    if mode == "meds":
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W1[i][0],np.median(medsum50W1[i]),std50W1_inf[i],std50W1_sup[i]) # init refers to the zeroth sample, all to all samples combined
        ax.text(0.02, 0.9, s, fontsize=fontlabel, color='b',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W2[i][0],np.median(medsum50W2[i]),std50W2_inf[i],std50W2_sup[i])
        ax.text(0.02, 0.7, s, fontsize=fontlabel, color='g',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W3[i][0],np.median(medsum50W3[i]),std50W3_inf[i],std50W3_sup[i])
        ax.text(0.02, 0.5, s, fontsize=fontlabel, color='r',transform=ax.transAxes)
        s = "init %.3f all %.3f 16th=%.3f 84th=%.3f" % (medsum50W4[i][0],np.median(medsum50W4[i]),std50W4_inf[i],std50W4_sup[i])
        ax.text(0.02, 0.3, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        s = "W1-4 init %.3f all %.3f 16th=%.3f 84th=%.3f" % (np.average([medsum50W1[i][0],medsum50W2[i][0],medsum50W3[i][0],medsum50W4[i][0]]),np.median([medsum50W1[i],medsum50W2[i],medsum50W3[i],medsum50W4[i]]),np.median([std50W1_inf[i],std50W2_inf[i],std50W3_inf[i],std50W4_inf[i]]),np.median([std50W1_sup[i],std50W2_sup[i],std50W3_sup[i],std50W4_sup[i]]))
        ax.text(0.02, 0.1, s, fontsize=fontlabel, color='k',transform=ax.transAxes)
        if np.median(medsum50W1[i]) == std50W1_inf[i]: std50W1_inf[i] = np.median(medsum50W1[i]) - 0.01
        if np.median(medsum50W2[i]) == std50W2_inf[i]: std50W2_inf[i] = np.median(medsum50W2[i]) - 0.01
        if np.median(medsum50W3[i]) == std50W3_inf[i]: std50W3_inf[i] = np.median(medsum50W3[i]) - 0.01
        if np.median(medsum50W4[i]) == std50W4_inf[i]: std50W4_inf[i] = np.median(medsum50W4[i]) - 0.01
        if np.median(medsum50W1[i]) == std50W1_sup[i]: std50W1_sup[i] = np.median(medsum50W1[i]) + 0.01
        if np.median(medsum50W2[i]) == std50W2_sup[i]: std50W2_sup[i] = np.median(medsum50W2[i]) + 0.01
        if np.median(medsum50W3[i]) == std50W3_sup[i]: std50W3_sup[i] = np.median(medsum50W3[i]) + 0.01
        if np.median(medsum50W4[i]) == std50W4_sup[i]: std50W4_sup[i] = np.median(medsum50W4[i]) + 0.01
        samplingW1[i] = sample(np.median(medsum50W1[i]), np.median(medsum50W1[i]) - std50W1_inf[i], std50W1_sup[i] - np.median(medsum50W1[i]))
        samplingW2[i] = sample(np.median(medsum50W2[i]), np.median(medsum50W2[i]) - std50W2_inf[i], std50W2_sup[i] - np.median(medsum50W2[i]))
        samplingW3[i] = sample(np.median(medsum50W3[i]), np.median(medsum50W3[i]) - std50W3_inf[i], std50W3_sup[i] - np.median(medsum50W3[i]))
        samplingW4[i] = sample(np.median(medsum50W4[i]), np.median(medsum50W4[i]) - std50W4_inf[i], std50W4_sup[i] - np.median(medsum50W4[i]))
        sampling[i] = np.r_[samplingW1[i],samplingW2[i],samplingW3[i],samplingW4[i]]
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

np.savetxt('%s%s_weightedcountshist_%s_%s_%s_%s_%s_%s_%sarcsec_percentiles_sampled.lst' % (root, lens, radius, mag, mode, photz, detect, irac, inner), sampling.T, fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')
plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.6)
plt.subplot(5,4,5)
plt.legend(bbox_to_anchor=(5, -5), loc='lower right', borderaxespad=0., fontsize=10)
plt.savefig('%s%s_weightedcountshist_%s_%s_%s_%s_%s_%s_%sarcsec_percentiles.png' % (root, lens, radius, mag, mode, photz, detect, irac, inner), dpi=500)
                
print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'

'''
Once all 8 runs have finished for a given radius, do the following:
! cat WFI2033_weightedcountshist_45_23_meds_bpz_deti_IRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_bpz_deti_noIRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_bpz_detir_IRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_bpz_detir_noIRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_eazy_deti_IRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_eazy_deti_noIRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_eazy_detir_IRAC_5arcsec_percentiles_sampled.lst WFI2033_weightedcountshist_45_23_meds_eazy_detir_noIRAC_5arcsec_percentiles_sampled.lst > WFI2033_weightedcountshist_45_23_meds_5arcsec_percentiles_sampled.lst
x = np.loadtxt("WFI2033_weightedcountshist_45_23_meds_5arcsec_percentiles_sampled.lst", unpack=True)
print np.median(x[0]), np.percentile(x[0], 16), np.percentile(x[0], 84)
print np.median(x[1]), np.percentile(x[1], 16), np.percentile(x[1], 84)
print np.median(x[2]), np.percentile(x[2], 16), np.percentile(x[2], 84)
print np.median(x[3]), np.percentile(x[3], 16), np.percentile(x[3], 84)
print np.median(x[4]), np.percentile(x[4], 16), np.percentile(x[4], 84)
print np.median(x[5]), np.percentile(x[5], 16), np.percentile(x[5], 84)
print np.median(x[6]), np.percentile(x[6], 16), np.percentile(x[6], 84)
print np.median(x[7]), np.percentile(x[7], 16), np.percentile(x[7], 84)
print np.median(x[8]), np.percentile(x[8], 16), np.percentile(x[8], 84)
print np.median(x[9]), np.percentile(x[9], 16), np.percentile(x[9], 84)
print np.median(x[10]), np.percentile(x[10], 16), np.percentile(x[10], 84)
print np.median(x[11]), np.percentile(x[11], 16), np.percentile(x[11], 84)
print np.median(x[12]), np.percentile(x[12], 16), np.percentile(x[12], 84)
print np.median(x[13]), np.percentile(x[13], 16), np.percentile(x[13], 84)
print np.median(x[14]), np.percentile(x[14], 16), np.percentile(x[14], 84)
print np.median(x[15]), np.percentile(x[15], 16), np.percentile(x[15], 84)
print np.median(x[16]), np.percentile(x[16], 16), np.percentile(x[16], 84)
print np.median(x[17]), np.percentile(x[17], 16), np.percentile(x[17], 84)
'''
