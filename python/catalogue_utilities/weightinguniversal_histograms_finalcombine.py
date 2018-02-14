# CE Rusu, Feb 13 2018
# Combines the results produced by weightinguniversal_histograms_samples.py into a final text file with the final distributions and widths
# run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_histograms_finalcombine.py WFI2033 meds 5 -1.0 -1.0 global handpicked
# At the moment the code only works for mag 23 (WFI2033)
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
try: handpicked = '_'+str(sys.argv[7])
except: handpicked = ''
rootin = "/Volumes/LaCieSubaru/weightedcounts/%s/" % lens
rootout = "/Users/cerusu/Dropbox/Davis_work/code/" % lens

# edit the conditions as desired, to restrict the included criteria:
lst45 = [x for x in os.listdir(rootin) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('45arcsec' in x) and ('%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('%s' %handpicked in x) # CHOOSE WHAT YOU WANT HERE
         and ('W1' in x) and ('W2' in x) and ('W3' in x) and ('W4' in x) and ('50' in x) and ('75' in x) and ('bpz' in x) and ('eazy' in x) and ('IRAC' in x) and ('noIRAC' in x) and ('deti' in x) and ('detir' in x) ]
lst120 = [x for x in os.listdir(root) if ('samples' in x) and ('.lst' in x) and ('%s_weightedcountshist_' %lens in x) and ('120arcsec' in x) and ('%sinner' %inner in x) and ('%s' %mode in x) and ('%s' %zinf in x) and ('%s' %zsup in x) and ('%s' %handpicked in x)
         and ('W1' in x) and ('W2' in x) and ('W3' in x) and ('W4' in x) and ('50' in x) and ('75' in x) and ('bpz' in x) and ('eazy' in x) and ('IRAC' in x) and ('noIRAC' in x) and ('deti' in x) and ('detir' in x) ]
for i in range(len(lst45)):
    if i == 0: x45 = np.loadtxt(lst45[0], unpack=True)
    else: x45 = np.r_[x45,np.loadtxt(lst45[i], unpack=True)]
for i in range(len(lst120)):
    if i == 0: x120 = np.loadtxt(lst120[0], unpack=True)
    else: x120 = np.r_[x120,np.loadtxt(lst120[i], unpack=True)]

if local == 'global':
    f = open('%s/weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(rootout,lens,type,inner,handpicked,zinf,zsup),'w')
    str = '# weight      45_23med    45_23inf  45_23sup 120_23med 120_23inf 120_23sup \n'
    str += 'gal           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84),np.median(x120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))
    str += 'z             %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84),np.median(x120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))
    str += 'mass          %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84),np.median(x120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))
    str += 'mass2         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84),np.median(x120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))
    str += 'mass3         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84),np.median(x120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))
    str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84),np.median(x120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))
    str += 'zoverr        %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84),np.median(x120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))
    str += 'massoverr     %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84),np.median(x120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))
    str += 'mass2overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84),np.median(x120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))
    str += 'mass3overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84),np.median(x120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))
    str += 'mass2rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84),np.median(x120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))
    str += 'mass3rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84),np.median(x120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))
    str += 'mass2overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84),np.median(x120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))
    str += 'mass3overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84),np.median(x120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))
    str += 'flexion       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84),np.median(x120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))
    str += 'tidal         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84),np.median(x120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))
    str += 'SIS           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84),np.median(x120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))
    str += 'SIShalo       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(x45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84),np.median(x120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))
    f.write(str)
    f.close()

else:
    fiducial45 = [x for x in os.listdir(rootin) if ('%s_weightedcountshist_45arcsec_%sinner_23_%s_bpz_detir_IRAC%s_zgap%s_%s_10samples' % (lens,inner,mag,mode,handpicked,zinf,zsup) in x) and ('.lst' in x)] # this will contain W1-W4 and 50/75; CHOOSE WHAT YOU WANT HERE
    fiducial120 = [x for x in os.listdir(rootin) if ('%s_weightedcountshist_120arcsec_%sinner_23_%s_bpz_detir_IRAC%s_zgap%s_%s_10samples' % (lens,inner,mag,mode,handpicked,zinf,zsup) in x) and ('.lst' in x)]
    for i in range(len(fiducial45)):
        if i == 0: f45 = np.loadtxt(fiducial45[0], unpack=True)
        else: f45 = np.r_[f45,np.loadtxt(fiducial45[i], unpack=True)]
    for i in range(len(fiducial120)):
        if i == 0: f120 = np.loadtxt(fiducial120[0], unpack=True)
        else: f120 = np.r_[f120,np.loadtxt(fiducial120[i], unpack=True)]

    f = open('%s/weightedcounts_%s_%s_%sinner%s_zgap%s_%s.cat' %(rootout,lens,type,inner,handpicked,zinf,zsup),'w')
    str = '# weight      45_23med    45_23inf  45_23sup 120_23med 120_23inf 120_23sup \n'
    str += 'gal           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[0]),np.percentile(x45[0], 16),np.percentile(x45[0], 84),np.median(f120[0]),np.percentile(x120[0], 16),np.percentile(x120[0], 84))
    str += 'z             %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[1]),np.percentile(x45[1], 16),np.percentile(x45[1], 84),np.median(f120[1]),np.percentile(x120[1], 16),np.percentile(x120[1], 84))
    str += 'mass          %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[2]),np.percentile(x45[2], 16),np.percentile(x45[2], 84),np.median(f120[2]),np.percentile(x120[2], 16),np.percentile(x120[2], 84))
    str += 'mass2         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[3]),np.percentile(x45[3], 16),np.percentile(x45[3], 84),np.median(f120[3]),np.percentile(x120[3], 16),np.percentile(x120[3], 84))
    str += 'mass3         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[4]),np.percentile(x45[4], 16),np.percentile(x45[4], 84),np.median(f120[4]),np.percentile(x120[4], 16),np.percentile(x120[4], 84))
    str += 'oneoverr      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[5]),np.percentile(x45[5], 16),np.percentile(x45[5], 84),np.median(f120[5]),np.percentile(x120[5], 16),np.percentile(x120[5], 84))
    str += 'zoverr        %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[6]),np.percentile(x45[6], 16),np.percentile(x45[6], 84),np.median(f120[6]),np.percentile(x120[6], 16),np.percentile(x120[6], 84))
    str += 'massoverr     %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[7]),np.percentile(x45[7], 16),np.percentile(x45[7], 84),np.median(f120[7]),np.percentile(x120[7], 16),np.percentile(x120[7], 84))
    str += 'mass2overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[8]),np.percentile(x45[8], 16),np.percentile(x45[8], 84),np.median(f120[8]),np.percentile(x120[8], 16),np.percentile(x120[8], 84))
    str += 'mass3overr    %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[9]),np.percentile(x45[9], 16),np.percentile(x45[9], 84),np.median(f120[9]),np.percentile(x120[9], 16),np.percentile(x120[9], 84))
    str += 'mass2rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[10]),np.percentile(x45[10], 16),np.percentile(x45[10], 84),np.median(f120[10]),np.percentile(x120[10], 16),np.percentile(x120[10], 84))
    str += 'mass3rms      %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[11]),np.percentile(x45[11], 16),np.percentile(x45[11], 84),np.median(f120[11]),np.percentile(x120[11], 16),np.percentile(x120[11], 84))
    str += 'mass2overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[12]),np.percentile(x45[12], 16),np.percentile(x45[12], 84),np.median(f120[12]),np.percentile(x120[12], 16),np.percentile(x120[12], 84))
    str += 'mass3overrrms %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[13]),np.percentile(x45[13], 16),np.percentile(x45[13], 84),np.median(f120[13]),np.percentile(x120[13], 16),np.percentile(x120[13], 84))
    str += 'flexion       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[14]),np.percentile(x45[14], 16),np.percentile(x45[14], 84),np.median(f120[14]),np.percentile(x120[14], 16),np.percentile(x120[14], 84))
    str += 'tidal         %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[15]),np.percentile(x45[15], 16),np.percentile(x45[15], 84),np.median(f120[15]),np.percentile(x120[15], 16),np.percentile(x120[15], 84))
    str += 'SIS           %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[16]),np.percentile(x45[16], 16),np.percentile(x45[16], 84),np.median(f120[16]),np.percentile(x120[16], 16),np.percentile(x120[16], 84))
    str += 'SIShalo       %.2f %.2f %.2f %.2f %.2f %.2f \n' % (np.median(f45[17]),np.percentile(x45[17], 16),np.percentile(x45[17], 84),np.median(f120[17]),np.percentile(x120[17], 16),np.percentile(x120[17], 84))
    f.write(str)
    f.close()

