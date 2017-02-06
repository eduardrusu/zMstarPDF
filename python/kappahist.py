
# run as: python kappahist.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_orig_size45_i24_ratioquick.lst constraints
#where constraints are eg: gal 1.52 0.05 oneoverr 1.69 0.05 mass 1.9 0.10 z 1.45 0.05 mass2 2.15 0.4 mass2rms 1.5 0.1 mass3 2.1 0.4 mass3rms 1.3 0.1 zoverr 1.6 0.05 massoverr 2.4 0.15 mass2overr 3.3 0.3 mass3overr 3.3 0.3 mass2overrrms 1.9 0.1 mass3overrrms 1.7 0.15 zmassoverr 2.4 0.1 zmass2overr 2.9 0.4

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

start_time=time.time()

q_gal=0
q_oneoverr=0
q_zweight=0
q_mass=0
q_mass2=0
q_mass2rms=0
q_mass3=0
q_mass3rms=0
q_zoverr=0
q_massoverr=0
q_mass2overr=0
q_mass3overr=0
q_mass2overrrms=0
q_mass3overrrms=0
q_zmassoverr=0
q_zmass2overr=0

degree=np.pi/180
L_field=4.0*degree
N_pix_per_dim = 4096
L_pix = L_field / N_pix_per_dim
os.system("ls /Volumes/G-RAIDStudio/simulations/lensing_simulations/ > /Volumes/G-RAIDStudio/simulations/lensing_simulations/ls.lst")
listname=[]
with open("/Volumes/G-RAIDStudio/simulations/lensing_simulations/ls.lst") as filelist:
    for l in filelist:
        listname=listname+[l[0:len(l)-1]]

for i in range(int(len(sys.argv)-1)/3):
    if str(sys.argv[int(3*i+2)])=="gal":
        q_gal=float(sys.argv[int(3*i+3)])
        q_gal_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="oneoverr":
        q_oneoverr=float(sys.argv[int(3*i+3)])
        q_oneoverr_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="z":
        q_zweight=float(sys.argv[int(3*i+3)])
        q_zweight_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass":
        q_mass=float(sys.argv[int(3*i+3)])
        q_mass_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass2":
        q_mass2=float(sys.argv[int(3*i+3)])
        q_mass2_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass2rms":
        q_mass2rms=float(sys.argv[int(3*i+3)])
        q_mass2rms_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass3":
        q_mass3=float(sys.argv[int(3*i+3)])
        q_mass3_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass3rms":
        q_mass3rms=float(sys.argv[int(3*i+3)])
        q_mass3rms_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="zoverr":
        q_zoverr=float(sys.argv[int(3*i+3)])
        q_zoverr_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="massoverr":
        q_massoverr=float(sys.argv[int(3*i+3)])
        q_massoverr_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass2overr":
        q_mass2overr=float(sys.argv[int(3*i+3)])
        q_mass2overr_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass3overr":
        q_mass3overr=float(sys.argv[int(3*i+3)])
        q_mass3overr_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass2overrrms":
        q_mass2overrrms=float(sys.argv[int(3*i+3)])
        q_mass2overrrms_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="mass3overrrms":
        q_mass3overrrm=float(sys.argv[int(3*i+3)])
        q_mass3overrrm_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="zmassoverr":
        q_zmassoverr=float(sys.argv[int(3*i+3)])
        q_gal_lim=float(sys.argv[int(3*i+4)])
    if str(sys.argv[int(3*i+2)])=="zmass2overr":
        q_zmass2overr=float(sys.argv[int(3*i+3)])
        q_zmass2overr_lim=float(sys.argv[int(3*i+4)])
kappa=np.array([])
#kappa_galonly=np.array([])
#kappa_noconstr=np.array([])
cnts=0
with open(str(sys.argv[1])) as file:
    for line in file:
        cnts=cnts+1
        if cnts%1000==0:
            print cnts
        ok=1
        #ok_galonly=1
        ID=line.split()[1]
        posx=float(line.split()[3])
        posy=float(line.split()[4])
        q_gal_file=float(line.split()[6])
        q_oneoverr_file=float(line.split()[8])
        q_zweight_file=float(line.split()[10])
        q_mass_file=float(line.split()[12])
        q_mass2_file=float(line.split()[14])
        q_mass2rms_file=float(line.split()[16])
        q_mass3_file=float(line.split()[18])
        q_mass3rms_file=float(line.split()[20])
        q_zoverr_file=float(line.split()[22])
        q_massoverr_file=float(line.split()[24])
        q_mass2overr_file=float(line.split()[26])
        q_mass3overr_file=float(line.split()[28])
        q_mass2overrrms_file=float(line.split()[30])
        q_mass3overrrms_file=float(line.split()[32])
        q_zmassoverr_file=float(line.split()[34])
        q_zmass2overr_file=float(line.split()[36])
        if (q_gal!=0) and ((q_gal_file < (q_gal - q_gal_lim)) or (q_gal_file > (q_gal + q_gal_lim))):
            ok=0
            #ok_galonly=0
        if (q_oneoverr!=0) and ((q_oneoverr_file < (q_oneoverr - q_oneoverr_lim)) or (q_oneoverr_file > (q_oneoverr + q_oneoverr_lim))):
            ok=0
        if (q_zweight!=0) and ((q_zweight_file < (q_zweight - q_zweight_lim)) or (q_zweight_file > (q_zweight + q_zweight_lim))):
            ok=0
        if (q_mass!=0) and ((q_mass_file < (q_mass - q_mass_lim)) or (q_mass_file > (q_mass + q_mass_lim))):
            ok=0
        if (q_mass2!=0) and ((q_mass2_file < (q_mass2 - q_mass2_lim)) or (q_mass2_file > (q_mass2 + q_mass2_lim))):
            ok=0
        if (q_mass2rms!=0) and ((q_mass2rms_file < (q_mass2rms - q_mass2rms_lim)) or (q_mass2rms_file > (q_mass2rms + q_mass2rms_lim))):
            ok=0
        if (q_mass3!=0) and ((q_mass3_file < (q_mass3 - q_mass3_lim)) or (q_mass3_file > (q_mass3 + q_mass3_lim))):
            ok=0
        if (q_mass3rms!=0) and ((q_mass3rms_file < (q_mass3rms - q_mass3rms_lim)) or (q_mass3rms_file > (q_mass3rms + q_mass3rms_lim))):
            ok=0
        if (q_zoverr!=0) and ((q_zoverr_file < (q_zoverr - q_zoverr_lim)) or (q_zoverr_file > (q_zoverr + q_zoverr_lim))):
            ok=0
        if (q_massoverr!=0) and ((q_massoverr_file < (q_massoverr - q_massoverr_lim)) or (q_massoverr_file > (q_massoverr + q_massoverr_lim))):
            ok=0
        if (q_mass2overr!=0) and ((q_mass2overr_file < (q_mass2overr - q_mass2overr_lim)) or (q_mass2overr_file > (q_mass2overr + q_mass2overr_lim))):
            ok=0
        if (q_mass3overr!=0) and ((q_mass3overr_file < (q_mass3overr - q_mass3overr_lim)) or (q_mass3overr_file > (q_mass3overr + q_mass3overr_lim))):
            ok=0
        if (q_mass2overrrms!=0) and ((q_mass2overrrms_file < (q_mass2overrrms - q_mass2overrrms_lim)) or (q_mass2overrrms_file > (q_mass2overrrms + q_mass2overrrms_lim))):
            ok=0
        if (q_mass3overrrms!=0) and ((q_mass3overrrms_file < (q_mass3overrrms - q_mass3overrrms_lim)) or (q_mass3overrrms_file > (q_mass3overrrms + q_mass3overrrms_lim))):
            ok=0
        if (q_zmassoverr!=0) and ((q_zmassoverr_file < (q_zmassoverr - q_zmassoverr_lim)) or (q_zmassoverr_file > (q_zmassoverr + q_zmassoverr_lim))):
            ok=0
        if (q_zmass2overr!=0) and ((q_zmass2overr_file < (q_zmass2overr - q_zmass2overr_lim)) or (q_zmass2overr_file > (q_zmass2overr + q_zmass2overr_lim))):
            ok=0
        for item in range(len(listname)):
            if (ID[75:80] in listname[item]) and ("dat" in listname[item]):
                kappafileend=listname[item]
        kappafile="/Volumes/G-RAIDStudio/simulations/lensing_simulations/"+kappafileend
        #print ID[75:80],ID,kappafileend
        x=int(round((posx + 0.5*L_field)/L_pix - 0.5))
        y=int(round((posy + 0.5*L_field)/L_pix - 0.5))
        count=int(round((x-1)*4096+y-1))
        if (ok==1):
            with open(kappafile) as kfile:
                l=kfile.readlines()
                kappa=np.append(kappa,[float(l[count].split()[1])])


#print "kappa=", kappa
if len(kappa)==0:
    print "No field satisfies the conditions."
else:
    print "Number of kappa elements:", len(kappa), " ------ Plotting histogram..."
    output=str(sys.argv[1])[0:len(str(sys.argv[1]))-4]
    for i in range(len(sys.argv)-1):
        if i>0:
            output=output+"_"+str(sys.argv[i+1])
    outputfile=output+"_kappa.dat"
    output=output+".eps"
    os.system("rm %s" %outputfile)
    for i in range(len(kappa)):
        kout=open('%s' %outputfile,'a')
        kout.write('%s \n' %kappa[i])
        kout.close()
    BINS=20
    plt.suptitle(r'%s'%output[141:len(output)-4], fontsize=13, y=0.998)
    #x = linspace(0,2,500)
    #plt.subplot(451)
    #n_q, bins_q, patches = plt.hist(q, histtype='step', color='b', label='W4sim', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
    plt.subplot(1,1,1)
    plt.hist(kappa, histtype='step', color='b', label='kappa w/ constr.', linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
    #plt.hist(kappa_galonly, histtype='step', color='r', label='kappa w/ gal constr.', linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
    #plt.hist(kappa_noconstr, histtype='step', color='g', label='kappa w/o constr.', linewidth=1, normed=1, bins=BINS, range=[-0.05, 0.2])
    #plt.hist(kappa_noconstr, histtype='step', color='b', label='kappa w/o constr.', linewidth=0.5, normed=1, bins=BINS, range=[0, rangemax])
    ax=plt.subplot(111)
    s = "med=%.3f, std=%.3f, cnt=%d" % (np.average(kappa),np.std(kappa),len(kappa))
    ax.text(0.55, 0.8, s, fontsize=10, color='b',transform=ax.transAxes)
    #text(0.5, 0.5,'matplotlib',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    #plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
    plt.ylabel("Normalized cnts", fontsize=20)
    plt.tick_params(axis='x', labelsize=13)
    plt.tick_params(axis='y', labelsize=13)
    plt.setp(plt.xticks()[1], rotation=90)
    plt.legend(bbox_to_anchor=(0.55, 0.7), loc='center left', borderaxespad=0., fontsize=10)
    #plt.subplots_adjust(top=0.6)
    #plt.tight_layout()
    plt.savefig(output, dpi=500)
    #plt.show()
print(" --- %s seconds ---" % (time.time() - start_time))

print 'Done!'


