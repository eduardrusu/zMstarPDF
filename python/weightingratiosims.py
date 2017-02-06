#B1608_sims Millenium
# 25 subsims of 1degx1deg, to emulate CFHTLENS W4
# cells: 4x4arcmin covering each subsim, in a grid
# run as: python weightingratiosims.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_orig_size120_i24.lst quick
#where speed is slow or quick

import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
from scipy import stats
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
from numpy import linspace

start_time_all = time.time()

os.system("wc %s > %s_count.lst" % (str(sys.argv[1]),str(sys.argv[1])[0:len(str(sys.argv[1]))-4])) # need to know how many lines in file
with open('%s_count.lst' % str(sys.argv[1])[0:len(str(sys.argv[1]))-4]) as count:
    for line in count:
        size=int(line.split()[0])
os.system("rm %s_count.lst" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4])
i=0
samples=5
with open(sys.argv[1]) as fileone:
    line=fileone.readlines()
    ID=line[0].split()[1]
    coordsx=line[0].split()[3]
    coordsy=line[0].split()[4]
with open(sys.argv[1]) as fileone:
    for line in fileone:
        if (ID==line.split()[1]) and (coordsx==line.split()[3]) and (coordsy==line.split()[4]):
            i=i+1
        else:
            samples=i
            break
                
if str(sys.argv[2])!="quick":
  os.system("rm %s_ratio.lst" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4])
  print "samples= ", samples
  i=0
  with open(sys.argv[1]) as fileone:
    for line in fileone:
        if (line!="\n"):
          start_time = time.time()
          ID=line.split()[1]
          #print "coordsx", line.split()[3], "coordsy", line.split()[4], "w_gal", line.split()[6], "w_oneoverr", line.split()[8]
          coordsx=float(line.split()[3])
          coordsy=float(line.split()[4])
          w_gal=float(line.split()[6])
          w_oneoverr=float(line.split()[8])
          #w_zweight=float(line.split()[10])
          w_mass=float(line.split()[12])
          #w_mass2=float(line.split()[14])
          #w_mass2rms=float(line.split()[16])
          #w_mass3=float(line.split()[18])
          #w_mass3rms=float(line.split()[20])
          #w_zoverr=float(line.split()[22])
          w_massoverr=float(line.split()[24])
          #w_mass2overr=float(line.split()[26])
          #w_mass3overr=float(line.split()[28])
          #w_mass2overrrms=float(line.split()[30])
          #w_mass3overrrms=float(line.split()[32])
          #w_zmassoverr=float(line.split()[34])
          #w_zmass2overr=float(line.split()[36])
          q_gal=np.array([])
          q_oneoverr=np.array([])
          #q_zweight=np.array([])
          q_mass=np.array([])
          #q_mass2=np.array([])
          #q_mass2rms=np.array([])
          #q_mass3=np.array([])
          #q_mass3rms=np.array([])
          #q_zoverr=np.array([])
          q_massoverr=np.array([])
          #q_mass2overr=np.array([])
          #q_mass3overr=np.array([])
          #q_mass2overrrms=np.array([])
          #q_mass3overrrms=np.array([])
          #q_zmassoverr=np.array([])
          #q_zmass2overr=np.array([])
          pos=i%samples # modulo
          with open(sys.argv[1]) as filetwo:
            linetwo=filetwo.readlines()
            while pos < size:
              if (linetwo!="\n"):
                if ID!=linetwo[pos].split()[1]:
                    q_gal=np.append(q_gal,[w_gal/float(linetwo[pos].split()[6])])
                    q_oneoverr=np.append(q_oneoverr,[w_oneoverr/float(linetwo[pos].split()[8])])
                    #q_zweight=np.append(q_zweight,[w_zweight/float(linetwo[pos].split()[10])])
                    q_mass=np.append(q_mass,[w_mass/float(linetwo[pos].split()[12])])
                    #q_mass2=np.append(q_mass2,[w_mass2/float(linetwo[pos].split()[14])])
                    #q_mass2rms=np.append(q_mass2rms,[w_mass2rms/float(linetwo[pos].split()[16])])
                    #q_mass3=np.append(q_mass3,[w_mass3/float(linetwo[pos].split()[18])])
                    #q_mass3rms=np.append(q_mass3rms,[w_mass3rms/float(linetwo[pos].split()[20])])
                    #q_zoverr=np.append(q_zoverr,[w_zoverr/float(linetwo[pos].split()[22])])
                    q_massoverr=np.append(q_massoverr,[w_massoverr/float(linetwo[pos].split()[24])])
                    #q_mass2overr=np.append(q_mass2overr,[w_mass2overr/float(linetwo[pos].split()[26])])
                    #q_mass3overr=np.append(q_mass3overr,[w_mass3overr/float(linetwo[pos].split()[28])])
                    #q_mass2overrrms=np.append(q_mass2overrrms,[w_mass2overrrms/float(linetwo[pos].split()[30])])
                    #q_mass3overrrms=np.append(q_mass3overrrms,[w_mass3overrrms/float(linetwo[pos].split()[32])])
                    #q_zmassoverr=np.append(q_zmassoverr,[w_zmassoverr/float(linetwo[pos].split()[34])])
                    #q_zmass2overr=np.append(q_zmass2overr,[w_zmass2overr/float(linetwo[pos].split()[36])])
                    #print pos
              pos=pos+samples
          i=i+1
          f=open('%s_ratio.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
          f.write('name= %s coords= %s %s q_gal_med= %.4e q_oneoverr_med= %.4e q_mass_med= %.4e q_massoverr_med= %.4e \n' % (ID, coordsx, coordsy, np.median(q_gal), np.median(q_oneoverr), np.median(q_mass), np.median(q_massoverr)))
        #f.write('name= %s q_gal_med= %.4e q_gal_ave= %.4e q_oneoverr_med= %.4e q_oneoverr_ave= %.4e q_zweight_med= %.4e q_zweight_ave= %.4e q_mass_med= %.4e q_mass_ave= %.4e q_mass2_med= %.4e q_mass2_ave= %.4e q_mass2rms_med= %.4e q_mass2rms_ave= %.4e q_mass3_med= %.4e q_mass3_ave= %.4e q_mass3rms_med= %.4e q_mass3rms_ave= %.4e q_zoverr_med= %.4e q_zoverr_ave= %.4e q_massoverr_med= %.4e q_massoverr_ave= %.4e q_mass2overr_med= %.4e q_mass2overr_ave= %.4e q_mass3overr_med= %.4e q_mass3overr_ave= %.4e q_mass2overrrms_med= %.4e q_mass2overrrms_ave= %.4e q_mass3overrrms_med= %.4e q_mass3overrrms_ave= %.4e q_zmassoverr_med= %.4e  q_zmassoverr_ave= %.4e q_zmass2overr_med= %.4e q_zmass2overr_ave= %.4e \n' % (ID, np.median(q_gal), np.average(q_gal), np.median(q_oneoverr), np.average(q_oneoverr), np.median(q_zweight), np.average(q_zweight), np.median(q_mass), np.average(q_mass), np.median(q_mass2), np.average(q_mass2), np.median(q_mass2rms), np.average(q_mass2rms), np.median(q_mass3), np.average(q_mass3), np.median(q_mass3rms), np.average(q_mass3rms), np.median(q_zoverr), np.average(q_zoverr), np.median(q_massoverr), np.average(q_massoverr), np.median(q_mass2overr), np.average(q_mass2overr), np.median(q_mass3overr), np.average(q_mass3overr), np.median(q_mass2overrrms), np.average(q_mass2overrrms), np.median(q_mass3overrrms), np.average(q_mass3overrrms), np.median(q_zmassoverr), np.average(q_zmassoverr), np.median(q_zmass2overr), np.average(q_zmass2overr)))
          f.close()
          #print i, ("time: --- %s seconds ---" % (time.time() - start_time))
          if i==50:
              print ("time per unit: --- %s seconds ---" % (time.time() - start_time))
          if i%100==0:
              print i,"/",size

  print("Total time: --- %s seconds ---" % (time.time() - start_time_all))




#QUICK WAY OF CALCULATING RATIOS, IN CASE WE DONT USE MASKS, SO WE CAN SIMPLY USE THE AVERAGE; CAREFUL THOUGH, TAKING THE MEDIAN AND THEN RATIOS MAY NOT BE MATHEMATICALLY EQUIVALENT

#ADD AS AN ARGUMENT WHETHER QUICK OR SLOW
# in the quick way I can compute for all weighting types

os.system("rm %s_ratioquick.lst" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4])
field_w_gal = np.loadtxt(sys.argv[1], usecols=[6], unpack=True)
field_w_oneoverr = np.loadtxt(sys.argv[1], usecols=[8], unpack=True)
field_w_zweight = np.loadtxt(sys.argv[1], usecols=[10], unpack=True)
field_w_mass = np.loadtxt(sys.argv[1], usecols=[12], unpack=True)
field_w_mass2 = np.loadtxt(sys.argv[1], usecols=[14], unpack=True)
field_w_mass2rms = np.loadtxt(sys.argv[1], usecols=[16], unpack=True)
field_w_mass3 = np.loadtxt(sys.argv[1], usecols=[18], unpack=True)
field_w_mass3rms = np.loadtxt(sys.argv[1], usecols=[20], unpack=True)
field_w_zoverr = np.loadtxt(sys.argv[1], usecols=[22], unpack=True)
field_w_massoverr = np.loadtxt(sys.argv[1], usecols=[24], unpack=True)
field_w_mass2overr = np.loadtxt(sys.argv[1], usecols=[26], unpack=True)
field_w_mass3overr = np.loadtxt(sys.argv[1], usecols=[28], unpack=True)
field_w_mass2overrrms = np.loadtxt(sys.argv[1], usecols=[30], unpack=True)
field_w_mass3overrrms = np.loadtxt(sys.argv[1], usecols=[32], unpack=True)
field_w_zmassoverr = np.loadtxt(sys.argv[1], usecols=[34], unpack=True)
field_w_zmass2overr = np.loadtxt(sys.argv[1], usecols=[36], unpack=True)
med_w_gal=np.median(field_w_gal)
med_w_oneoverr=np.median(field_w_oneoverr)
med_w_zweight=np.median(field_w_zweight)
med_w_mass=np.median(field_w_mass)
med_w_mass2=np.median(field_w_mass2)
med_w_mass2rms=np.median(field_w_mass2rms)
med_w_mass3=np.median(field_w_mass3)
med_w_mass3rms=np.median(field_w_mass3rms)
med_w_zoverr=np.median(field_w_zoverr)
med_w_massoverr=np.median(field_w_massoverr)
med_w_mass2overr=np.median(field_w_mass2overr)
med_w_mass3overr=np.median(field_w_mass3overr)
med_w_mass2overrrms=np.median(field_w_mass2overrrms)
med_w_mass3overrrms=np.median(field_w_mass3overrrms)
med_w_zmassoverr=np.median(field_w_zmassoverr)
med_w_zmass2overr=np.median(field_w_zmass2overr)
ave_w_gal=np.average(field_w_gal)
ave_w_oneoverr=np.average(field_w_oneoverr)
ave_w_zweight=np.average(field_w_zweight)
ave_w_mass=np.average(field_w_mass)
ave_w_mass2=np.average(field_w_mass2)
ave_w_mass2rms=np.average(field_w_mass2rms)
ave_w_mass3=np.average(field_w_mass3)
ave_w_mass3rms=np.average(field_w_mass3rms)
ave_w_zoverr=np.average(field_w_zoverr)
ave_w_massoverr=np.average(field_w_massoverr)
ave_w_mass2overr=np.average(field_w_mass2overr)
ave_w_mass3overr=np.average(field_w_mass3overr)
ave_w_mass2overrrms=np.average(field_w_mass2overrrms)
ave_w_mass3overrrms=np.average(field_w_mass3overrrms)
ave_w_zmassoverr=np.average(field_w_zmassoverr)
ave_w_zmass2overr=np.average(field_w_zmass2overr)
with open(sys.argv[1]) as fileone:
    for line in fileone:
        if (line!="\n"):
          ID=line.split()[1]
          coordsx=float(line.split()[3])
          coordsy=float(line.split()[4])
          w_gal=float(line.split()[6])
          w_oneoverr=float(line.split()[8])
          w_zweight=float(line.split()[10])
          w_mass=float(line.split()[12])
          w_mass2=float(line.split()[14])
          w_mass2rms=float(line.split()[16])
          w_mass3=float(line.split()[18])
          w_mass3rms=float(line.split()[20])
          w_zoverr=float(line.split()[22])
          w_massoverr=float(line.split()[24])
          w_mass2overr=float(line.split()[26])
          w_mass3overr=float(line.split()[28])
          w_mass2overrrms=float(line.split()[30])
          w_mass3overrrms=float(line.split()[32])
          w_zmassoverr=float(line.split()[34])
          w_zmass2overr=float(line.split()[36])
          q_gal_med=w_gal/med_w_gal
          q_oneoverr_med=w_oneoverr/med_w_oneoverr
          q_zweight_med=w_zweight/med_w_zweight
          q_mass_med=w_mass/med_w_mass
          q_mass2_med=w_mass2/med_w_mass2
          q_mass2rms_med=w_mass2rms/med_w_mass2rms
          q_mass3_med=w_mass3/med_w_mass3
          q_mass3rms_med=w_mass3rms/med_w_mass3rms
          q_zoverr_med=w_zoverr/med_w_zoverr
          q_massoverr_med=w_massoverr/med_w_massoverr
          q_mass2overr_med=w_mass2overr/med_w_mass2overr
          q_mass3overr_med=w_mass3overr/med_w_mass3overr
          q_mass2overrrms_med=w_mass2overrrms/med_w_mass2overrrms
          q_mass3overrrms_med=w_mass3overrrms/med_w_mass3overrrms
          q_zmassoverr_med=w_zmassoverr/med_w_zmassoverr
          q_zmass2overr_med=w_zmass2overr/med_w_zmass2overr
          q_gal_ave=w_gal/ave_w_gal
          q_oneoverr_ave=w_oneoverr/ave_w_oneoverr
          q_zweight_ave=w_zweight/ave_w_zweight
          q_mass_ave=w_mass/ave_w_mass
          q_mass2_ave=w_mass2/ave_w_mass2
          q_mass2rms_ave=w_mass2rms/ave_w_mass2rms
          q_mass3_ave=w_mass3/ave_w_mass3
          q_mass3rms_ave=w_mass3rms/ave_w_mass3rms
          q_zoverr_ave=w_zoverr/ave_w_zoverr
          q_massoverr_ave=w_massoverr/ave_w_massoverr
          q_mass2overr_ave=w_mass2overr/ave_w_mass2overr
          q_mass3overr_ave=w_mass3overr/ave_w_mass3overr
          q_mass2overrrms_ave=w_mass2overrrms/ave_w_mass2overrrms
          q_mass3overrrms_ave=w_mass3overrrms/ave_w_mass3overrrms
          q_zmassoverr_ave=w_zmassoverr/ave_w_zmassoverr
          q_zmass2overr_ave=w_zmass2overr/ave_w_zmass2overr
          f=open('%s_ratioquick.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
          f.write('name= %s coords= %s %s q_gal_med= %.4e q_oneoverr_med= %.4e q_zweight_med= %.4e q_mass_med= %.4e q_mass2_med= %.4e q_mass2rms_med= %.4e q_mass3_med= %.4e q_mass3rms_med= %.4e q_zoverr_med= %.4e q_massoverr_med= %.4e q_mass2overr_med= %.4e q_mass3overr_med= %.4e q_mass2overrrms_med= %.4e q_mass3overrrms_med= %.4e q_zmassoverr_med= %.4e  q_zmass2overr_med= %.4e q_gal_ave= %.4e q_oneoverr_ave= %.4e q_zweight_ave= %.4e q_mass_ave= %.4e q_mass2_ave= %.4e q_mass2rms_ave= %.4e q_mass3_ave= %.4e q_mass3rms_ave= %.4e q_zoverr_ave= %.4e q_massoverr_ave= %.4e q_mass2overr_ave= %.4e q_mass3overr_ave= %.4e q_mass2overrrms_ave= %.4e q_mass3overrrms_ave= %.4e q_zmassoverr_ave= %.4e q_zmass2overr_ave= %.4e \n' % (ID, coordsx, coordsy, q_gal_med, q_oneoverr_med, q_zweight_med, q_mass_med, q_mass2_med, q_mass2rms_med, q_mass3_med, q_mass3rms_med, q_zoverr_med, q_massoverr_med, q_mass2overr_med, q_mass3overr_med, q_mass2overrrms_med, q_mass3overrrms_med, q_zmassoverr_med, q_zmass2overr_med, q_gal_ave, q_oneoverr_ave, q_zweight_ave, q_mass_ave, q_mass2_ave, q_mass2rms_ave, q_mass3_ave, q_mass3rms_ave, q_zoverr_ave, q_massoverr_ave, q_mass2overr_ave, q_mass3overr_ave, q_mass2overrrms_ave, q_mass3overrrms_ave, q_zmassoverr_ave, q_zmass2overr_ave))
          f.close()

# select the first field encountered with overdensity similar to the lens, and plot histogram

lens_q_gal=1.52
lens_q_gal_lim=0.05

if str(sys.argv[2])=="quick":
    filename='%s_ratioquick.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4]
else:
    filename='%s_ratio.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4]
with open(filename) as file:
    for line in file:
        if (float(line.split()[6])>=lens_q_gal-lens_q_gal_lim) and (float(line.split()[6])<=lens_q_gal+lens_q_gal_lim):
            IDmatch=line.split()[1]
            coordsxmatch=float(line.split()[3])
            coordsymatch=float(line.split()[4])
            break

i=0
with open(sys.argv[1]) as fileone:
    for line in fileone:
        if (line!="\n"):
          start_time = time.time()
          ID=line.split()[1]
          coordsx=float(line.split()[3])
          coordsy=float(line.split()[4])
          if ((ID==IDmatch) and (coordsx==coordsxmatch) and (coordsy==coordsymatch)):
            w_gal=float(line.split()[6])
            w_oneoverr=float(line.split()[8])
            w_zweight=float(line.split()[10])
            w_mass=float(line.split()[12])
            w_mass2=float(line.split()[14])
            w_mass2rms=float(line.split()[16])
            w_mass3=float(line.split()[18])
            w_mass3rms=float(line.split()[20])
            w_zoverr=float(line.split()[22])
            w_massoverr=float(line.split()[24])
            w_mass2overr=float(line.split()[26])
            w_mass3overr=float(line.split()[28])
            w_mass2overrrms=float(line.split()[30])
            w_mass3overrrms=float(line.split()[32])
            w_zmassoverr=float(line.split()[34])
            w_zmass2overr=float(line.split()[36])
            q_gal=np.array([])
            q_oneoverr=np.array([])
            q_zweight=np.array([])
            q_mass=np.array([])
            q_mass2=np.array([])
            q_mass2rms=np.array([])
            q_mass3=np.array([])
            q_mass3rms=np.array([])
            q_zoverr=np.array([])
            q_massoverr=np.array([])
            q_mass2overr=np.array([])
            q_mass3overr=np.array([])
            q_mass2overrrms=np.array([])
            q_mass3overrrms=np.array([])
            q_zmassoverr=np.array([])
            q_zmass2overr=np.array([])
            pos=i%samples # modulo
            with open(sys.argv[1]) as filetwo:
                linetwo=filetwo.readlines()
                while pos < size:
                    if (linetwo!="\n"):
                        if ID!=linetwo[pos].split()[1]:
                            q_gal=np.append(q_gal,[w_gal/float(linetwo[pos].split()[6])])
                            q_oneoverr=np.append(q_oneoverr,[w_oneoverr/float(linetwo[pos].split()[8])])
                            q_zweight=np.append(q_zweight,[w_zweight/float(linetwo[pos].split()[10])])
                            q_mass=np.append(q_mass,[w_mass/float(linetwo[pos].split()[12])])
                            q_mass2=np.append(q_mass2,[w_mass2/float(linetwo[pos].split()[14])])
                            q_mass2rms=np.append(q_mass2rms,[w_mass2rms/float(linetwo[pos].split()[16])])
                            q_mass3=np.append(q_mass3,[w_mass3/float(linetwo[pos].split()[18])])
                            q_mass3rms=np.append(q_mass3rms,[w_mass3rms/float(linetwo[pos].split()[20])])
                            q_zoverr=np.append(q_zoverr,[w_zoverr/float(linetwo[pos].split()[22])])
                            q_massoverr=np.append(q_massoverr,[w_massoverr/float(linetwo[pos].split()[24])])
                            q_mass2overr=np.append(q_mass2overr,[w_mass2overr/float(linetwo[pos].split()[26])])
                            q_mass3overr=np.append(q_mass3overr,[w_mass3overr/float(linetwo[pos].split()[28])])
                            q_mass2overrrms=np.append(q_mass2overrrms,[w_mass2overrrms/float(linetwo[pos].split()[30])])
                            q_mass3overrrms=np.append(q_mass3overrrms,[w_mass3overrrms/float(linetwo[pos].split()[32])])
                            q_zmassoverr=np.append(q_zmassoverr,[w_zmassoverr/float(linetwo[pos].split()[34])])
                            q_zmass2overr=np.append(q_zmass2overr,[w_zmass2overr/float(linetwo[pos].split()[36])])
                            #print pos
                    pos=pos+samples
                i=i+1
                for j in range(len(q_gal)):
                    f=open('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
                    #f.write('name= %s q_gal_med= %.4e q_oneoverr_med= %.4e q_mass_med= %.4e q_massoverr_med= %.4e \n' % (ID, np.median(q_gal), np.median(q_oneoverr), np.median(q_mass), np.median(q_massoverr)))
                    f.write('name= %s coords= %s %s q_gal_med= %.4e q_oneoverr_med= %.4e q_zweight_med= %.4e q_mass_med= %.4e q_mass2_med= %.4e q_mass2rms_med= %.4e q_mass3_med= %.4e q_mass3rms_med= %.4e q_zoverr_med= %.4e q_massoverr_med= %.4e q_mass2overr_med= %.4e q_mass3overr_med= %.4e q_mass2overrrms_med= %.4e q_mass3overrrms_med= %.4e q_zmassoverr_med= %.4e q_zmass2overr_med= %.4e \n' % (ID, coordsx, coordsy, q_gal[j], q_oneoverr[j], q_zweight[j], q_mass[j], q_mass2[j], q_mass2rms[j], q_mass3[j], q_mass3rms[j], q_zoverr[j], q_massoverr[j], q_mass2overr[j], q_mass3overr[j], q_mass2overrrms[j], q_mass3overrrms[j], q_zmassoverr[j], q_zmass2overr[j]))
                    f.close()
                #print i, ("time: --- %s seconds ---" % (time.time() - start_time))
                if i==50:
                    print ("time per unit: --- %s seconds ---" % (time.time() - start_time))
                if i%100==0:
                    print i,"/",size
print("Total time: --- %s seconds ---" % (time.time() - start_time_all))
print "Ploting histogram..."
BINS=50

cols=6
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.suptitle(r'Selected simulated weight histogram', fontsize=10, y=0.998)

gauss_q = gaussian_kde(q)


x = linspace(0,2,500)

plt.subplot(451)
rangemax=4
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


#plt.plot(x,gauss_q(x),'b', linewidth=0.5)
ax=plt.subplot(451)
#s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q(x))],np.average(q),np.median(q))
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)

plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 1
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=8
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


x = linspace(0,2,500)

gauss_q = gaussian_kde(q)


plt.subplot(452)
rangemax=4
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


#plt.plot(x,gauss_q(x),'b', linewidth=0.5)
ax=plt.subplot(452)
#s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q(x))],np.average(q),np.median(q))
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)

plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 2
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=10
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(453)
rangemax=4
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(453)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 3
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=12
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(454)
rangemax=4
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(454)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 4
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=14
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(456)
rangemax=7
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(456)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 5
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=16
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(457)
rangemax=6
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(457)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 6
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=18
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(458)
rangemax=6
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(458)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 7
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=20
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(459)
rangemax=7
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(459)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 8
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=22
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,11)
rangemax=4
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,11)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 9
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=24
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,12)
rangemax=7
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,12)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 10
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=26
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,13)
rangemax=6
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,13)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 11
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=28
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,14)
rangemax=6
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,14)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 12
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=30
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


x = linspace(0,3,500)

gauss_q = gaussian_kde(q)


plt.subplot(4,5,16)
rangemax=8
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


#plt.plot(x,gauss_q(x),'b', linewidth=0.5)
ax=plt.subplot(4,5,16)
#s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q(x))],np.average(q),np.median(q))
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)

plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 13
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=32
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


x = linspace(0,3,500)

gauss_q = gaussian_kde(q)


plt.subplot(4,5,17)
rangemax=8
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,17)
#s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q(x))],np.average(q),np.median(q))
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q)) # only the peak depends on the binning
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 14
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=34
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,18)
rangemax=7
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,18)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 15
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

cols=36
qread = np.loadtxt('%s_selecthist.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], usecols=[cols], unpack=True)
q = qread[qread < 10]


plt.subplot(4,5,19)
rangemax=10
n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])


ax=plt.subplot(4,5,19)
s = "%.3f,%.3f,%.3f" % (bins_q[np.argmax(n_q)],np.average(q),np.median(q))
ax.text(0.15, 0.8, s, fontsize=5, color='b',transform=ax.transAxes)


plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Normalized cnts", fontsize=7)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

subplot = 19
print "finished subplot %d/16; fraction of points inside the q < 10 cut: \n W4sim %.3f" % (subplot, float(q.size)/qread.size)

plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()

plt.savefig('%s_selecthist.eps' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4], dpi=500)

#plt.show()

print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'


