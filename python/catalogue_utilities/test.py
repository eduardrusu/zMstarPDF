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
samples = 10
limit = 10**30
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
    print '%s/%s' %(nr,samples-1)
    lstW1_50 = [x for x in os.listdir(root) if ('W1' in x) and ('_wghtratios_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)] # select from the files in the root directory
    lstW1_75 = [x for x in os.listdir(root) if ('W1' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW2_50 = [x for x in os.listdir(root) if ('W2' in x) and ('_wghtratios_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW2_75 = [x for x in os.listdir(root) if ('W2' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW3_50 = [x for x in os.listdir(root) if ('W3' in x) and ('_wghtratios_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW3_75 = [x for x in os.listdir(root) if ('W3' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW4_50 = [x for x in os.listdir(root) if ('W4' in x) and ('_wghtratios_50_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]
    lstW4_75 = [x for x in os.listdir(root) if ('W4' in x) and ('_wghtratios_75_msk%sarcsecrad%sarcsecgap_%s_%s_%s_%s_%s_%s_zgap%s_%s%s_%s%s.lst' %(radius,inner,lens,mag,photz,detect,irac,mode,zinf,zsup,handpicked,str(nr),specialtest) in x)]

    print "W1..."
    for i in range(len(lstW1_50)):
        hdu = fits.open(root+lstW1_50[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W1_50read = dataread
        else: q_W1_50read = np.r_['1',q_W1_50read,dataread]
        hdu.close()
        #print np.shape(q_W1_50read)
    for i in range(len(lstW1_75)):
        hdu = fits.open(root+lstW1_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W1_75read = dataread
        else: q_W1_75read = np.r_['1',q_W1_75read,dataread]
        hdu.close()

    print "W2..."
    for i in range(len(lstW2_50)):
        hdu = fits.open(root+lstW2_50[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W2_50read = dataread
        else: q_W2_50read = np.r_['1',q_W2_50read,dataread]
        hdu.close()
        #print np.shape(q_W2_50read)
    for i in range(len(lstW2_75)):
        hdu = fits.open(root+lstW2_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W2_75read = dataread
        else: q_W2_75read = np.r_['1',q_W2_75read,dataread]
        hdu.close()

    print "W3..."
    for i in range(len(lstW3_50)):
        hdu = fits.open(root+lstW3_50[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W3_50read = dataread
        else: q_W3_50read = np.r_['1',q_W3_50read,dataread]
        hdu.close()
        #print np.shape(q_W3_50read)
    for i in range(len(lstW3_75)):
        hdu = fits.open(root+lstW3_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W3_75read = dataread
        else: q_W3_75read = np.r_['1',q_W3_75read,dataread]
        hdu.close()

    print "W4..."
    for i in range(len(lstW4_50)):
        hdu = fits.open(root+lstW4_50[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W4_50read = dataread
        else: q_W4_50read = np.r_['1',q_W4_50read,dataread]
        hdu.close()
        #print np.shape(q_W4_50read)
    for i in range(len(lstW4_75)):
        hdu = fits.open(root+lstW4_75[i]); data = hdu[1].data
        dataread = np.c_[data.field('5_lens_gal'),data.field('6_lens_zweight'),data.field('7_lens_mass'),data.field('8_lens_mass2'),data.field('9_lens_mass3'),data.field('10_lens_oneoverr'),data.field('11_lens_zoverr'),data.field('12_lens_massoverr'),data.field('13_lens_mass2overr'),data.field('14_lens_mass3overr'),data.field('15_lens_mass2rms'),data.field('16_lens_mass3rms'),data.field('17_lens_mass2overrms'),data.field('18_lens_mass3overrms'),data.field('19_lens_flexion'),data.field('20_lens_tidal'),data.field('21_lens_convergence'),data.field('22_lens_convergencehalo')].T
        if i == 0:
            q_W4_75read = dataread
        else: q_W4_75read = np.r_['1',q_W4_75read,dataread]
        hdu.close()
