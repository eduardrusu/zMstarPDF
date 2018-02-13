# The code uses the weighted count ratios derived by weightinguniversal_overlap_sampling_nobeta_WFI2033rethought.py to produce paper-quality histograms without overlapped text
# run as # run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_histograms_publicationqualitynotext.py WFI2033 45 5 23 meds bpz deti IRAC 0.61 0.71 100 handpicked

import numpy as np
import sys
import os
import time
import matplotlib.pyplot as plt

plt.clf()

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

fontabsciss = 14
fontlabel = 4
pltrange = 2.5
if radius == "45":
    vertlimit = 2.5
else:
    vertlimit = 2.5
if mode == "sum":
    vertlimit = 2.5
limit = 10**30
root = "/Volumes/LaCieSubaru/weightedcounts/%s/" % lens

start_time = time.time()

print "Reading..."
''' if I want to plot the benchmark'''
#lstW1_50 = [x for x in os.listdir(root) if ('W1' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_0.lst' %(radius,lens,detect,irac,mode,inner) in x)] # select from the files in the root directory
#lstW2_50 = [x for x in os.listdir(root) if ('W2' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_0.lst' %(radius,lens,detect,irac,mode,inner) in x)]
#lstW3_50 = [x for x in os.listdir(root) if ('W3' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_0.lst' %(radius,lens,detect,irac,mode,inner) in x)]
#lstW4_50 = [x for x in os.listdir(root) if ('W4' in x) and ('_24galphotmstar_50_%s_%s_%s_%s_%s_%sarcsec_0.lst' %(radius,lens,detect,irac,mode,inner) in x)]
''' if I want to plot everything together'''
lstW1_50 = [x for x in os.listdir(root) if ('W1' in x) and ('_24galphotmstar_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,detect,irac,mode,zinf,zsup,handpicked) in x)] # select from the files in the root directory
lstW2_50 = [x for x in os.listdir(root) if ('W2' in x) and ('_24galphotmstar_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,detect,irac,mode,zinf,zsup,handpicked) in x)]
lstW3_50 = [x for x in os.listdir(root) if ('W3' in x) and ('_24galphotmstar_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,detect,irac,mode,zinf,zsup,handpicked) in x)]
lstW4_50 = [x for x in os.listdir(root) if ('W4' in x) and ('_24galphotmstar_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_zgap%s_%s%s' %(radius,inner,lens,detect,irac,mode,zinf,zsup,handpicked) in x)]

if mag == "24" and photz == "bpz": cols=[4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38]
if mag == "24" and photz == "eazy": cols=[40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74]
if mag == "23" and photz == "bpz": cols=[5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]
if mag == "23" and photz == "eazy": cols=[41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75]

for i in range(len(lstW1_50)):
        if i == 0:
            q_W1_50read = np.loadtxt(root+lstW1_50[i], usecols=cols, unpack=True)
        else:
            q_W1_50read = np.r_['1',q_W1_50read,np.loadtxt(root+ lstW1_50[i], usecols=cols, unpack=True)]

for i in range(len(lstW2_50)):
        if i == 0:
            q_W2_50read = np.loadtxt(root+lstW2_50[i], usecols=cols, unpack=True)
        else:
            q_W2_50read = np.r_['1',q_W2_50read,np.loadtxt(root+ lstW2_50[i], usecols=cols, unpack=True)]

for i in range(len(lstW3_50)):
        if i == 0:
            q_W3_50read = np.loadtxt(root+lstW3_50[i], usecols=cols, unpack=True)
        else:
            q_W3_50read = np.r_['1',q_W3_50read,np.loadtxt(root+ lstW3_50[i], usecols=cols, unpack=True)]

for i in range(len(lstW4_50)):
        if i == 0:
            q_W4_50read = np.loadtxt(root+lstW4_50[i], usecols=cols, unpack=True)
        else:
            q_W4_50read = np.r_['1',q_W4_50read,np.loadtxt(root+ lstW4_50[i], usecols=cols, unpack=True)]

print "Plotting..."

fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(6,3,1)
ax1.set_aspect(1, adjustable='datalim')

for i in range(18):
    
    q_W1_50 = q_W1_50read[i][q_W1_50read[i] < limit]
    if mode == "sum": q_W1_50 = abs(q_W1_50) # fix the negative halo convergence
    q_W2_50 = q_W2_50read[i][q_W2_50read[i] < limit]
    if mode == "sum": q_W2_50 = abs(q_W2_50)
    q_W3_50 = q_W3_50read[i][q_W3_50read[i] < limit]
    if mode == "sum": q_W3_50 = abs(q_W3_50)
    q_W4_50 = q_W4_50read[i][q_W4_50read[i] < limit]
    if mode == "sum": q_W4_50 = abs(q_W4_50)
    
    if i == 0: ax=plt.subplot(6,3,1, sharex=ax1, sharey=ax1)
    if i == 1: ax=plt.subplot(6,3,2, sharex=ax1, sharey=ax1)
    if i == 2: ax=plt.subplot(6,3,3, sharex=ax1, sharey=ax1)
    if i == 3: ax=plt.subplot(6,3,4, sharex=ax1, sharey=ax1)
    if i == 4: ax=plt.subplot(6,3,5, sharex=ax1, sharey=ax1)
    if i == 5: ax=plt.subplot(6,3,6, sharex=ax1, sharey=ax1)
    if i == 6: ax=plt.subplot(6,3,7, sharex=ax1, sharey=ax1)
    if i == 7: ax=plt.subplot(6,3,8, sharex=ax1, sharey=ax1)
    if i == 8: ax=plt.subplot(6,3,9, sharex=ax1, sharey=ax1)
    if i == 9: ax=plt.subplot(6,3,10, sharex=ax1, sharey=ax1)
    if i == 10: ax=plt.subplot(6,3,11, sharex=ax1, sharey=ax1)
    if i == 11: ax=plt.subplot(6,3,12, sharex=ax1, sharey=ax1)
    if i == 12: ax=plt.subplot(6,3,13, sharex=ax1, sharey=ax1)
    if i == 13: ax=plt.subplot(6,3,14, sharex=ax1, sharey=ax1)
    if i == 14: ax=plt.subplot(6,3,15, sharex=ax1, sharey=ax1)
    if i == 15: ax=plt.subplot(6,3,16, sharex=ax1, sharey=ax1)
    if i == 16: ax=plt.subplot(6,3,17, sharex=ax1, sharey=ax1)
    if i == 17: ax=plt.subplot(6,3,18, sharex=ax1, sharey=ax1)
    
    plt.hist(q_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])
    plt.hist(q_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=bin, range=[0, pltrange])

    z=np.linspace(0,vertlimit,1000)
    w=np.ones(1000)
    plt.plot(w,z,'k:')
    w=np.ones(1000) * np.median(q_W1_50)
    plt.plot(w,z,'b--')
    w=np.ones(1000) * np.median(q_W2_50)
    plt.plot(w,z,'g--')
    w=np.ones(1000) * np.median(q_W3_50)
    plt.plot(w,z,'r--')
    w=np.ones(1000) * np.median(q_W4_50)
    plt.plot(w,z,'k--')
    plt.xlim(0, 2.5)
    plt.ylim(0, vertlimit)
    if vertlimit == 4:
        ax.set_yticklabels(np.arange(0.0, vertlimit, 0.5))
    else:
        ax.set_yticklabels(np.arange(0.0, vertlimit, 0.5))
    ax.set_xticklabels(np.arange(0.0, 2.5, 0.5))
    if i == 0: ax.text(0.7, 0.7, '$1$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 1: ax.text(0.7, 0.7, '$z$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 2: ax.text(0.7, 0.7, '$M_\star$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 3: ax.text(0.7, 0.7, '$M^2_\star$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 4: ax.text(0.7, 0.7, '$M^3_\star$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 5: ax.text(0.7, 0.7, '$1/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 6: ax.text(0.7, 0.7, '$z/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 7: ax.text(0.7, 0.7, '$M_\star/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 8: ax.text(0.7, 0.7, '$M^2_\star/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 9: ax.text(0.7, 0.7, '$M^3_\star/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 10: ax.text(0.7, 0.7, '$M^2_{\star\mathrm{,rms}}$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 11: ax.text(0.7, 0.7, '$M^3_{\star\mathrm{,rms}}$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 12: ax.text(0.7, 0.7, '$M^2_\star/r_\mathrm{,rms}$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 13: ax.text(0.7, 0.7, '$M^3_\star/r_\mathrm{,rms}$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 14: ax.text(0.7, 0.7, '$M_\star/r^3$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 15: ax.text(0.7, 0.7, '$M_\star/r^2$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 16: ax.text(0.7, 0.7, '$\sqrt{M_\star}/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    if i == 17: ax.text(0.7, 0.7, '$\sqrt{M_h}/r$', fontsize=fontabsciss, color='black',transform=ax.transAxes)
    #ax.set_xticklabels(np.arange(0.0, 3.0, 0.5))
    #if i in [0,4,8,12,15]:
        #plt.ylabel("Normalized counts", fontsize=5)

    plt.tick_params(axis='x', labelsize=10)
    plt.tick_params(axis='y', labelsize=10)
    plt.setp(plt.xticks()[1], rotation=90)
    subplot = i+1
    print "finished subplot %d/18; fraction of points inside the < %s cut: W1_75 %.3f \n W2_75 %.3f \n W3_75 %.3f \n W4_75 %.3f " % (subplot, limit, float(q_W1_50.size)/q_W1_50read[0].size, float(q_W2_50.size)/q_W2_50read[0].size, float(q_W3_50.size)/q_W3_50read[0].size, float(q_W4_50.size)/q_W4_50read[0].size)

ax=plt.subplot(6,3,1, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,2, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,3, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,4, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,5, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,6, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,7, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,8, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,9, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,10, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,11, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,12, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,13, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,14, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,15, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
plt.tick_params(axis='x', labelbottom='off')
ax=plt.subplot(6,3,17, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
ax=plt.subplot(6,3,18, sharex=ax1, sharey=ax1)
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

plt.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.98, wspace=0, hspace=0)
#plt.subplot(6,3,5)
#plt.legend(bbox_to_anchor=(1, 1),loc='lower right', borderaxespad=0., fontsize=10)
#plt.legend(bbox_to_anchor=(5, -5), loc='lower right', borderaxespad=0., fontsize=10)

# for some reason I need to add the lines below, or the last plot is not displayed
#plt.subplots_adjust(top=0.6)
#plt.tight_layout()

#fig.text(0.5, 0.05, '$\zeta^\mathrm{meds,WX}_{q}$', ha='center', va='center', size='20')
if ((radius == "45") & (mag == "23") & (mode == "meds")):
    fig.text(0.5, 0.05, r"$\zeta^\mathrm{meds,WX}_{q,45'',i<23}$", ha='center', va='center', size='20')
if ((radius == "120") & (mag == "23") & (mode == "meds")):
    fig.text(0.5, 0.05, r"$\zeta^\mathrm{meds,WX}_{q,120'',i<23}$", ha='center', va='center', size='20')
if ((radius == "120") & (mag == "24") & (mode == "sum")):
    fig.text(0.5, 0.05, r"$\zeta^\mathrm{sum,WX}_{q,120'',i<24}$", ha='center', va='center', size='20')
#fig.text(0.05, 0.5, 'normalized counts', ha='center', va='center', size='20', rotation='vertical')
plt.savefig('%s%s_weightedcountshist_%sarcsec_%sinner_%s_%s_%s_%s_%s%s_zgap%s_%s_notext.png' % (root, lens, radius, inner, mag, mode, photz, detect, irac, handpicked, zinf, zsup), dpi=500)

print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
