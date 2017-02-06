#B1608_sims Millenium
# 25 subsims of 1degx1deg, to emulate CFHTLENS W4
# cells: 4x4arcmin covering each subsim, in a grid
# run as: python weightingsims.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_8_7_6_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs.cat samplesize# radius maglimit
# where samplesize number is 0, 100 or 1000 and msk_lenssize is 45, 60, 90 or 120; maglimit is 23 23.5 or 24;
#since the weights depend on z_s, I will not run 4 times but instead in the output file output for all 4 redshifts in different columns

import numpy as np
import scipy
import sys
from scipy import special
from scipy import stats
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time

start_time = time.time()

print("Arguments: Catalogue name: %s \n Number of samples to be drawn from P(z) and P(Mstar): %s \n Radius of each cell: %s \n Limiting magnitude: %s" % (str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4])))

deginrad=0.0174532925
zgrid=np.linspace(0.05,3.5,70)
zgridint = np.arange(70) # because stats.rv_discrete only works with integer points
z_s_B1608 = 1.39
z_s_HE0435 = 1.69
z_s_HE1104 = 2.32
z_s_RX1131 = 0.66

limitx,limity = np.loadtxt(sys.argv[1], usecols=[4,5], unpack=True)
cells_on_a_side = [int((round(np.max(limitx/deginrad)-round(np.min(limitx/deginrad))))*3600/240),int((round(np.max(limity/deginrad)-round(np.min(limity/deginrad))))*3600/240)]
print "Position limits:", round(np.min(limitx/deginrad)),round(np.max(limitx/deginrad)),round(np.min(limity/deginrad)),round(np.max(limity/deginrad))
centersx=np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
centersy=np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
minx=round(np.min(limitx/deginrad))
miny=round(np.min(limity/deginrad))
for i in range(cells_on_a_side[0]):
    for j in range(cells_on_a_side[1]):
        centersx[i][j]=minx * deginrad + i * deginrad/cells_on_a_side[0] + deginrad/cells_on_a_side[0]/2
        centersy[i][j]=miny * deginrad + j * deginrad/cells_on_a_side[1] + deginrad/cells_on_a_side[1]/2
#print centersx,centersy
B1608_specsim_gal = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_oneoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_zweight = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_mass = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_mass2 = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_mass3 = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_zoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_massoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_mass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_mass3overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_zmassoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_specsim_zmass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_gal = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_oneoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_zweight = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_mass = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_mass2 = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_mass3 = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_zoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_massoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_mass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_mass3overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_zmassoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_origsim_zmass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1]))
B1608_sim_gal = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_oneoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_zweight = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_mass = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_mass2 = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_mass3 = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_zoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_massoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_mass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_mass3overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_zmassoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_sim_zmass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_gal = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_oneoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_zweight = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_mass = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_mass2 = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_mass3 = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_zoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_massoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_mass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_mass3overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_zmassoverr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))
B1608_tab_sim_zmass2overr = np.zeros((cells_on_a_side[0],cells_on_a_side[1],int(str(sys.argv[2]))))

with open(sys.argv[1]) as file:  #B1608_sims#.lst
    for gal in file:
        i=i+1
        if i in [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000]:
            print i, "objects..."
        if (gal!="\n"):
            ID=str(gal.split()[0])
            posx = float(gal.split()[4])
            posy = float(gal.split()[5])
            MAG_i = float(gal.split()[12])
            catmass = np.log10(float(gal.split()[1]))
            catz = float(gal.split()[2])
            Z_B = float(gal.split()[3])
            Z_B_MIN = Z_B - 0.13 # TEMPORARY solution until i include this value in the _pdz catalogue
            Z_B_MAX = Z_B + 0.13
            Mstar_best_catz = float(gal.split()[93])
            Mstar_inf_catz = float(gal.split()[94])
            Mstar_med_catz = float(gal.split()[95])
            Mstar_sup_catz = float(gal.split()[96])
            Mstar_best_zb = float(gal.split()[97])
            Mstar_inf_zb = float(gal.split()[98])
            Mstar_med_zb = float(gal.split()[99])
            Mstar_sup_zb = float(gal.split()[100])
            x=int((posx/deginrad-minx)*cells_on_a_side[0])
            y=int((posy/deginrad-miny)*cells_on_a_side[1])
            #print x,y #posx/deginrad/cells_on_a_side[0],posy/deginrad/cells_on_a_side[1],minx,miny
            #print ID,posx,posy,MAG_i,catmass,catz,Z_B,Z_B_MIN,Z_B_MAX,Mstar_best_catz,Mstar_inf_catz,Mstar_med_catz, Mstar_sup_catz, Mstar_best_zb, Mstar_inf_zb,Mstar_med_zb,Mstar_sup_zb
            sep=np.sqrt((centersx[x][y]-posx)**2 + (centersy[x][y]-posy)**2)*1/deginrad*3600 # in arcsec
            if (MAG_i <= int(str(sys.argv[4]))) and (sep >= 4) and (sep <= int(str(sys.argv[3]))):

                #spec
                if catz <= z_s_B1608:
                    B1608_specsim_gal[x][y] = B1608_specsim_gal[x][y] + 1
                    B1608_specsim_zweight[x][y] = B1608_specsim_zweight[x][y] + (z_s_B1608 * catz) - (catz * catz)
                    B1608_specsim_mass[x][y] = B1608_specsim_mass[x][y] + 10**Mstar_med_catz
                    B1608_specsim_mass2[x][y] = B1608_specsim_mass2[x][y] + ((10**Mstar_med_catz) * (10**Mstar_med_catz))
                    B1608_specsim_mass3[x][y] = B1608_specsim_mass3[x][y] + ((10**Mstar_med_catz) * (10**Mstar_med_catz) * (10**Mstar_med_catz))
                    if (sep <= 10):
                        B1608_specsim_oneoverr[x][y] = B1608_specsim_oneoverr[x][y] + 0.1
                        B1608_specsim_zoverr[x][y] = B1608_specsim_zoverr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) / 10)
                        B1608_specsim_massoverr[x][y] = B1608_specsim_massoverr[x][y] + ((10**Mstar_med_catz) / 10)
                        B1608_specsim_mass2overr[x][y] = B1608_specsim_mass2overr[x][y] + (((10**Mstar_med_catz) * (10**Mstar_med_catz)) / 10)
                        B1608_specsim_mass3overr[x][y] = B1608_specsim_mass3overr[x][y] + (((10**Mstar_med_catz) * (10**Mstar_med_catz) * (10**Mstar_med_catz)) / 10)
                        B1608_specsim_zmassoverr[x][y] = B1608_specsim_zmassoverr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) * (10**Mstar_med_catz) / 10)
                        B1608_specsim_zmass2overr[x][y] = B1608_specsim_zmass2overr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) * (10**Mstar_med_catz) * (10**Mstar_med_catz) / 10)
                    else:
                        B1608_specsim_oneoverr[x][y] = B1608_specsim_oneoverr[x][y] + 1. / sep
                        B1608_specsim_zoverr[x][y] = B1608_specsim_zoverr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) / sep)
                        B1608_specsim_massoverr[x][y] = B1608_specsim_massoverr[x][y] + ((10**Mstar_med_catz) / sep)
                        B1608_specsim_mass2overr[x][y] = B1608_specsim_mass2overr[x][y] + (((10**Mstar_med_catz) * (10**Mstar_med_catz)) / sep)
                        B1608_specsim_mass3overr[x][y] = B1608_specsim_mass3overr[x][y] + (((10**Mstar_med_catz) * (10**Mstar_med_catz) * (10**Mstar_med_catz)) / sep)
                        B1608_specsim_zmassoverr[x][y] = B1608_specsim_zmassoverr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) * (10**Mstar_med_catz) / sep)
                        B1608_specsim_zmass2overr[x][y] = B1608_specsim_zmass2overr[x][y] + (((z_s_B1608 * catz) - (catz * catz)) * (10**Mstar_med_catz) * (10**Mstar_med_catz) / sep)

                #orig:
                where = 100 + int(np.max([Z_B/0.05, 1])) * 4 #the column in the catalogue
                if Mstar_med_zb < 0:
                    if float(gal.split()[where-1])>0: #mass_med
                        Mstar_med_zb=float(gal.split()[where-1])
                    elif float(gal.split()[where-3])>0: #mass_best
                        Mstar_med_zb=float(gal.split()[where-3])
                    else:
                        pos=int(np.max([Z_B/0.05, 1]))
                        while (pos<68) and (float(gal.split()[100 + pos * 4 - 3])<=0):
                            pos=pos+1
                        Mstar_med_zb=float(gal.split()[100 + pos * 4 - 3])
                        if Mstar_med_zb<=0:
                            pos=int(np.max([Z_B/0.05, 1]))
                            while (pos>1) and (float(gal.split()[100 + pos * 4 - 3])<=0):
                                pos=pos-1
                            Mstar_med_zb=float(gal.split()[100 + pos * 4 - 3])
                    #print Mstar_med_zb
                    if Mstar_med_zb<=0:
                        Mstar_med_zb=9
                if (Mstar_inf_zb < 0) or (Mstar_sup_zb < 0):
                    Mstar_inf_zb=Mstar_med_zb-0.1
                    Mstar_sup_zb=Mstar_med_zb+0.1
                if Z_B <= z_s_B1608:
                   B1608_origsim_gal[x][y] = B1608_origsim_gal[x][y] + 1
                   B1608_origsim_zweight[x][y] = B1608_origsim_zweight[x][y] + (z_s_B1608 * Z_B) - (Z_B * Z_B)
                   B1608_origsim_mass[x][y] = B1608_origsim_mass[x][y] + 10**Mstar_med_zb
                   B1608_origsim_mass2[x][y] = B1608_origsim_mass2[x][y] + ((10**Mstar_med_zb) * (10**Mstar_med_zb))
                   B1608_origsim_mass3[x][y] = B1608_origsim_mass3[x][y] + ((10**Mstar_med_zb) * (10**Mstar_med_zb) * (10**Mstar_med_zb))
                   if (sep <= 10):
                       B1608_origsim_oneoverr[x][y] = B1608_origsim_oneoverr[x][y] + 0.1
                       B1608_origsim_zoverr[x][y] = B1608_origsim_zoverr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) / 10)
                       B1608_origsim_massoverr[x][y] = B1608_origsim_massoverr[x][y] + ((10**Mstar_med_zb) / 10)
                       B1608_origsim_mass2overr[x][y] = B1608_origsim_mass2overr[x][y] + (((10**Mstar_med_zb) * (10**Mstar_med_zb)) / 10)
                       B1608_origsim_mass3overr[x][y] = B1608_origsim_mass3overr[x][y] + (((10**Mstar_med_zb) * (10**Mstar_med_zb) * (10**Mstar_med_zb)) / 10)
                       B1608_origsim_zmassoverr[x][y] = B1608_origsim_zmassoverr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) * (10**Mstar_med_zb) / 10)
                       B1608_origsim_zmass2overr[x][y] = B1608_origsim_zmass2overr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) * (10**Mstar_med_zb) * (10**Mstar_med_zb) / 10)
                   else:
                       B1608_origsim_oneoverr[x][y] = B1608_origsim_oneoverr[x][y] + 1. / sep
                       B1608_origsim_zoverr[x][y] = B1608_origsim_zoverr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) / sep)
                       B1608_origsim_massoverr[x][y] = B1608_origsim_massoverr[x][y] + ((10**Mstar_med_zb) / sep)
                       B1608_origsim_mass2overr[x][y] = B1608_origsim_mass2overr[x][y] + (((10**Mstar_med_zb) * (10**Mstar_med_zb)) / sep)
                       B1608_origsim_mass3overr[x][y] = B1608_origsim_mass3overr[x][y] + (((10**Mstar_med_zb) * (10**Mstar_med_zb) * (10**Mstar_med_zb)) / sep)
                       B1608_origsim_zmassoverr[x][y] = B1608_origsim_zmassoverr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) * (10**Mstar_med_zb) / sep)
                       B1608_origsim_zmass2overr[x][y] = B1608_origsim_zmass2overr[x][y] + (((z_s_B1608 * Z_B) - (Z_B * Z_B)) * (10**Mstar_med_zb) * (10**Mstar_med_zb) / sep)







                                                                  
                samplesup=np.zeros(int(str(sys.argv[2])))
                sampleinf=np.zeros(int(str(sys.argv[2])))
                fracz=1-1.0*(Z_B_MAX-Z_B)/(Z_B_MAX-Z_B_MIN)
                fracmass=1-(1.0*(10**Mstar_sup_zb-10**Mstar_med_zb)/(10**Mstar_sup_zb-10**Mstar_inf_zb))

                #if str(sys.argv[8]) == "samp":
                if (fracz > 0) and (fracz < 1):
                    samplesup=Z_B+abs(np.random.normal(0, Z_B_MAX-Z_B, int(str(sys.argv[2]))))
                    sampleinf=Z_B-abs(np.random.normal(0, Z_B-Z_B_MIN, int(str(sys.argv[2]))))
                        # no negative redshifts
                    while len(sampleinf[sampleinf<0]) > 0:
                        sampleinf[sampleinf<0]=Z_B-abs(np.random.normal(0, Z_B-Z_B_MIN, len(sampleinf[sampleinf<0])))
                    rand=np.random.random(int(str(sys.argv[2])))
                    samplez=sampleinf
                    samplez[np.where(rand>fracz)]=samplesup[np.where(rand>fracz)]
                    #samplez[samplez<0]=0
                if (fracz <= 0) or (fracz >= 1):
                    samplez=np.random.normal(Z_B_MIN+(Z_B_MAX-Z_B_MIN)/2, (Z_B_MAX-Z_B_MIN)/2, int(str(sys.argv[2])))
                    #samplez[samplez<0]=0
                # for mass, assume gaussian distribution in log, not normal space, so I don't get things like negative mass
                if (fracmass > 0) and (fracmass < 1):
                    samplesup=10**(Mstar_med_zb+abs(np.random.normal(0, Mstar_sup_zb-Mstar_med_zb, int(str(sys.argv[2])))))
                    sampleinf=10**(Mstar_med_zb-abs(np.random.normal(0, Mstar_med_zb-Mstar_inf_zb, int(str(sys.argv[2])))))
                    rand=np.random.random(int(str(sys.argv[2])))
                    samplemass=sampleinf
                    samplemass[np.where(rand>fracmass)]=samplesup[np.where(rand>fracmass)]
                if (fracmass <= 0) or (fracmass >= 1):
                    samplemass=10**np.random.normal(Mstar_inf_zb+(Mstar_sup_zb-Mstar_inf_zb)/2, (Mstar_sup_zb-Mstar_inf_zb)/2, int(str(sys.argv[2])))
                #print "samplemass"
                #print samplemass
                #ignore objects with z>z_source for all weights when drawing from the z PDF
                B1608_sim_gal[x][y][samplez<z_s_B1608] = B1608_sim_gal[x][y][samplez<z_s_B1608] + 1
                B1608_sim_zweight[x][y][samplez<z_s_B1608] = B1608_sim_zweight[x][y][samplez<z_s_B1608] + (z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)
                B1608_sim_mass[x][y][samplez<z_s_B1608] = B1608_sim_mass[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]
                B1608_sim_mass2[x][y][samplez<z_s_B1608] = B1608_sim_mass2[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]**2
                B1608_sim_mass3[x][y][samplez<z_s_B1608] = B1608_sim_mass3[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]**3
                #print (catsim['Z_B_MAX'][i] - catsim['Z_B_MIN'][i])/2
                if (sep <= 10):
                   B1608_sim_oneoverr[x][y][samplez<z_s_B1608] = B1608_sim_oneoverr[x][y][samplez<z_s_B1608] + 0.1
                   B1608_sim_zoverr[x][y][samplez<z_s_B1608] = B1608_sim_zoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) / 10)
                   B1608_sim_massoverr[x][y][samplez<z_s_B1608] = B1608_sim_massoverr[x][y][samplez<z_s_B1608] + (samplemass[samplez<z_s_B1608] / 10)
                   B1608_sim_mass2overr[x][y][samplez<z_s_B1608] = B1608_sim_mass2overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**2) / 10)
                   B1608_sim_mass3overr[x][y][samplez<z_s_B1608] = B1608_sim_mass3overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**3) / 10)
                   B1608_sim_zmassoverr[x][y][samplez<z_s_B1608] = B1608_sim_zmassoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * samplemass[samplez<z_s_B1608] / 10)
                   B1608_sim_zmass2overr[x][y][samplez<z_s_B1608] = B1608_sim_zmass2overr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * (samplemass[samplez<z_s_B1608]**2) / 10)
                else:
                   B1608_sim_oneoverr[x][y][samplez<z_s_B1608] = B1608_sim_oneoverr[x][y][samplez<z_s_B1608] + 1/sep
                   B1608_sim_zoverr[x][y][samplez<z_s_B1608] = B1608_sim_zoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) / sep)
                   B1608_sim_massoverr[x][y][samplez<z_s_B1608] = B1608_sim_massoverr[x][y][samplez<z_s_B1608] + (samplemass[samplez<z_s_B1608] / sep)
                   B1608_sim_mass2overr[x][y][samplez<z_s_B1608] = B1608_sim_mass2overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**2) / sep)
                   B1608_sim_mass3overr[x][y][samplez<z_s_B1608] = B1608_sim_mass3overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**3) / sep)
                   B1608_sim_zmassoverr[x][y][samplez<z_s_B1608] = B1608_sim_zmassoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * samplemass[samplez<z_s_B1608] / sep)
                   B1608_sim_zmass2overr[x][y][samplez<z_s_B1608] = B1608_sim_zmass2overr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * (samplemass[samplez<z_s_B1608]**2) / sep)
   
                #if str(sys.argv[8]) == "tab":
                massbest_tab=np.zeros(70)
                massinf_tab=np.zeros(70)
                massmed_tab=np.zeros(70)
                masssup_tab=np.zeros(70)
                pdz_tab=np.zeros(70)
                for m in range(70):
                    pdz_tab[m]=float(gal.split()[22+m])
                    if pdz_tab[m]<0.001:
                        pdz_tab[m]=0
                if len(pdz_tab[pdz_tab!=0])!=0:
                    for m in range(70):
                        if pdz_tab[m]!=0:
                            massbest_tab[m]=float(gal.split()[100 + (m+1) * 4 - 3])
                            #print massbest_tab[m]
                        if massbest_tab[m]<0:
                            massbest_tab[m]=9 # very small number of exceptions
                    for m in range(70):
                        if pdz_tab[m]!=0:
                            massmed_tab[m]=float(gal.split()[100 + (m+1) * 4 - 1])
                            if massmed_tab[m]<0:
                                massmed_tab[m]=massbest_tab[m]
                            massinf_tab[m]=float(gal.split()[100 + (m+1) * 4 - 2])
                            masssup_tab[m]=float(gal.split()[100 + (m+1) * 4])
                            if massinf_tab[m]<0:
                                massinf_tab[m]=massmed_tab[m]-0.1
                            if masssup_tab[m]<0:
                                masssup_tab[m]=massmed_tab[m]+0.1
                            #print ID,m,massbest_tab[m],massinf_tab[m],masssup_tab[m],masssup_tab[m]
                    custm = stats.rv_discrete(name='custm', values=(zgridint, pdz_tab))
                    sample=custm.rvs(size=int(str(sys.argv[2])))
                    iter=0
                    while len(sample[pdz_tab[sample]==0]) != 0: # happens because the probabilities do not sum exactly to 1; first reshuffle 10 times; if this does not solve the problem replace with the value having maximum probability; happens because the probabilities do not sum exactly to 1
                        iter=iter+1
                        sample=custm.rvs(size=int(str(sys.argv[2])))
                        if iter==10:
                            #print pdz_tab,ID,sample,str(sys.argv[1])
                            sample[pdz_tab[sample]==0]=np.where(pdz_tab==np.max(pdz_tab[sample[pdz_tab[sample]!=0]]))[0][0]
                    samplez=zgrid[sample]
                    sample_massinf_tab=massinf_tab[sample] # since "sample" is constant, this insures that Mstar corresponds to z
                    sample_massmed_tab=massmed_tab[sample]
                    #print sample_massmed_tab
                    sample_masssup_tab=masssup_tab[sample]
                    sample_lenssup=np.zeros(int(str(sys.argv[2])))
                    sample_lensinf=np.zeros(int(str(sys.argv[2])))
                    #print sample_masssup_tab, sample_massmed_tab, samplez, catfield['ID'][0]
                    iter=0
                    for l in range(int(str(sys.argv[2]))):
                        if sample_massmed_tab[l]==0:     #I SHOULD NOT HAVE TO DO THIS, THERE IS A BUG; actually I tested and this happens because of missing stellar mass for z=3.5 (the edge)
                            sample_massmed_tab[l]=9
                            sample_massinf_tab[l]=8.9
                            sample_masssup_tab[l]=9.1
                            if iter==0:
                                print "Exception!", ID
                            iter=iter+1
                        sample_lenssup[l]=10**(sample_massmed_tab[l]+abs(np.random.normal(0, sample_masssup_tab[l]-sample_massmed_tab[l], 1)))
                        sample_lensinf[l]=10**(sample_massmed_tab[l]-abs(np.random.normal(0, sample_massmed_tab[l]-sample_massinf_tab[l], 1)))
                        rand=np.random.random(1)
                        samplemass[l]=sample_lensinf[l]
                        fracmass=1-(1.0*(10**sample_masssup_tab[l]-10**sample_massmed_tab[l])/(10**sample_masssup_tab[l]-10**sample_massinf_tab[l]))
                        if rand>fracmass:
                            samplemass[l]=sample_lenssup[l]
                    #print zgrid[sample]
                    #print massinf_tab[sample]
                    #print massmed_tab[sample]
                    #print masssup_tab[sample]
                    #print samplemass
                    #ignore objects with z>z_source for all weights when drawing from the z PDF
                    B1608_tab_sim_gal[x][y][samplez<z_s_B1608] =  B1608_tab_sim_gal[x][y][samplez<z_s_B1608] + 1
                    B1608_tab_sim_zweight[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zweight[x][y][samplez<z_s_B1608] + (z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)
                    B1608_tab_sim_mass[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]
                    B1608_tab_sim_mass2[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass2[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]**2
                    B1608_tab_sim_mass3[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass3[x][y][samplez<z_s_B1608] + samplemass[samplez<z_s_B1608]**3
                    #print (catfield['Z_B_MAX'][i] - catfield['Z_B_MIN'][i])/2
                    if (sep <= 10):
                        B1608_tab_sim_oneoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_oneoverr[x][y][samplez<z_s_B1608] + 0.1
                        B1608_tab_sim_zoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) / 10)
                        B1608_tab_sim_massoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_massoverr[x][y][samplez<z_s_B1608] + (samplemass[samplez<z_s_B1608] / 10)
                        B1608_tab_sim_mass2overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass2overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**2) / 10)
                        B1608_tab_sim_mass3overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass3overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**3) / 10)
                        B1608_tab_sim_zmassoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zmassoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * samplemass[samplez<z_s_B1608] / 10)
                        B1608_tab_sim_zmass2overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zmass2overr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * (samplemass[samplez<z_s_B1608]**2) / 10)
                    else:
                        B1608_tab_sim_oneoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_oneoverr[x][y][samplez<z_s_B1608] + 1/sep
                        B1608_tab_sim_zoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) / sep)
                        B1608_tab_sim_massoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_massoverr[x][y][samplez<z_s_B1608] + (samplemass[samplez<z_s_B1608] / sep)
                        B1608_tab_sim_mass2overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass2overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**2) / sep)
                        B1608_tab_sim_mass3overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_mass3overr[x][y][samplez<z_s_B1608] + ((samplemass[samplez<z_s_B1608]**3) / sep)
                        B1608_tab_sim_zmassoverr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zmassoverr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * samplemass[samplez<z_s_B1608] / sep)
                        B1608_tab_sim_zmass2overr[x][y][samplez<z_s_B1608] =  B1608_tab_sim_zmass2overr[x][y][samplez<z_s_B1608] + (((z_s_B1608 * samplez[samplez<z_s_B1608]) - (samplez[samplez<z_s_B1608]**2)) * (samplemass[samplez<z_s_B1608]**2) / sep)

B1608_specsim_mass2rms = np.sqrt(B1608_origsim_mass2)
B1608_specsim_mass3rms = scipy.special.cbrt(B1608_origsim_mass3)
B1608_specsim_mass2overrrms = np.sqrt(B1608_origsim_mass2overr)
B1608_specsim_mass3overrrms = scipy.special.cbrt(B1608_origsim_mass3overr)
#if str(sys.argv[8]) == "orig":
B1608_origsim_mass2rms = np.sqrt(B1608_origsim_mass2)
B1608_origsim_mass3rms = scipy.special.cbrt(B1608_origsim_mass3)
B1608_origsim_mass2overrrms = np.sqrt(B1608_origsim_mass2overr)
B1608_origsim_mass3overrrms = scipy.special.cbrt(B1608_origsim_mass3overr)
#if str(sys.argv[8]) == "samp":
B1608_sim_mass2rms = np.sqrt(B1608_sim_mass2)
B1608_sim_mass3rms = scipy.special.cbrt(B1608_sim_mass3)
B1608_sim_mass2overrrms = np.sqrt(B1608_sim_mass2overr)
B1608_sim_mass3overrrms = scipy.special.cbrt(B1608_sim_mass3overr)
#if str(sys.argv[8]) == "tab":
B1608_tab_sim_mass2rms = np.sqrt(B1608_tab_sim_mass2)
B1608_tab_sim_mass3rms = scipy.special.cbrt(B1608_tab_sim_mass3)
B1608_tab_sim_mass2overrrms = np.sqrt(B1608_tab_sim_mass2overr)
B1608_tab_sim_mass3overrrms = scipy.special.cbrt(B1608_tab_sim_mass3overr)

s = open('%s_B1608_spec_size%s_i%s.lst' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4],str(sys.argv[3]),str(sys.argv[4])),'w')
f = open('%s_B1608_orig_size%s_i%s.lst' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4],str(sys.argv[3]),str(sys.argv[4])),'w')
#if str(sys.argv[8]) == "samp":
g = open('%s_B1608_samp_size%s_i%s.lst' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4],str(sys.argv[3]),str(sys.argv[4])),'w')
#if str(sys.argv[8]) == "tab":
h = open('%s_B1608_tab_size%s_i%s.lst' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4],str(sys.argv[3]),str(sys.argv[4])),'w')
for i in range(cells_on_a_side[0]):
    for j in range(cells_on_a_side[1]):
            s.write('name= %s coords= %s %s w_gal= %.4e w_oneoverr= %.4e w_zweight= %.4e w_mass= %.4e w_mass2= %.4e w_mass2rms= %.4e w_mass3= %.4e w_mass3rms= %.4e w_zoverr= %.4e w_massoverr= %.4e w_mass2overr= %.4e w_mass3overr= %.4e w_mass2overrrms= %.4e w_mass3overrrms= %.4e w_zmassoverr= %.4e w_zmass2overr= %.4e \n' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4], centersx[i][j], centersy[i][j], B1608_specsim_gal[i][j], B1608_specsim_oneoverr[i][j], B1608_specsim_zweight[i][j], B1608_specsim_mass[i][j], B1608_specsim_mass2[i][j], B1608_specsim_mass2rms[i][j], B1608_specsim_mass3[i][j], B1608_specsim_mass3rms[i][j], B1608_specsim_zoverr[i][j], B1608_specsim_massoverr[i][j], B1608_specsim_mass2overr[i][j], B1608_specsim_mass3overr[i][j], B1608_specsim_mass2overrrms[i][j], B1608_specsim_mass3overrrms[i][j], B1608_specsim_zmassoverr[i][j], B1608_specsim_zmass2overr[i][j]))
            f.write('name= %s coords= %s %s w_gal= %.4e w_oneoverr= %.4e w_zweight= %.4e w_mass= %.4e w_mass2= %.4e w_mass2rms= %.4e w_mass3= %.4e w_mass3rms= %.4e w_zoverr= %.4e w_massoverr= %.4e w_mass2overr= %.4e w_mass3overr= %.4e w_mass2overrrms= %.4e w_mass3overrrms= %.4e w_zmassoverr= %.4e w_zmass2overr= %.4e \n' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4], centersx[i][j], centersy[i][j], B1608_origsim_gal[i][j], B1608_origsim_oneoverr[i][j], B1608_origsim_zweight[i][j], B1608_origsim_mass[i][j], B1608_origsim_mass2[i][j], B1608_origsim_mass2rms[i][j], B1608_origsim_mass3[i][j], B1608_origsim_mass3rms[i][j], B1608_origsim_zoverr[i][j], B1608_origsim_massoverr[i][j], B1608_origsim_mass2overr[i][j], B1608_origsim_mass3overr[i][j], B1608_origsim_mass2overrrms[i][j], B1608_origsim_mass3overrrms[i][j], B1608_origsim_zmassoverr[i][j], B1608_origsim_zmass2overr[i][j]))
            for k in range(int(str(sys.argv[2]))):
                g.write('name= %s coords= %s %s w_gal= %.4e w_oneoverr= %.4e w_zweight= %.4e w_mass= %.4e w_mass2= %.4e w_mass2rms= %.4e w_mass3= %.4e w_mass3rms= %.4e w_zoverr= %.4e w_massoverr= %.4e w_mass2overr= %.4e w_mass3overr= %.4e w_mass2overrrms= %.4e w_mass3overrrms= %.4e w_zmassoverr= %.4e w_zmass2overr= %.4e \n' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4], centersx[i][j], centersy[i][j], B1608_sim_gal[i][j][k], B1608_sim_oneoverr[i][j][k], B1608_sim_zweight[i][j][k], B1608_sim_mass[i][j][k], B1608_sim_mass2[i][j][k], B1608_sim_mass2rms[i][j][k], B1608_sim_mass3[i][j][k], B1608_sim_mass3rms[i][j][k], B1608_sim_zoverr[i][j][k], B1608_sim_massoverr[i][j][k], B1608_sim_mass2overr[i][j][k], B1608_sim_mass3overr[i][j][k], B1608_sim_mass2overrrms[i][j][k], B1608_sim_mass3overrrms[i][j][k], B1608_sim_zmassoverr[i][j][k], B1608_sim_zmass2overr[i][j][k]))
            for k in range(int(str(sys.argv[2]))):
                h.write('name= %s coords= %s %s w_gal= %.4e w_oneoverr= %.4e w_zweight= %.4e w_mass= %.4e w_mass2= %.4e w_mass2rms= %.4e w_mass3= %.4e w_mass3rms= %.4e w_zoverr= %.4e w_massoverr= %.4e w_mass2overr= %.4e w_mass3overr= %.4e w_mass2overrrms= %.4e w_mass3overrrms= %.4e w_zmassoverr= %.4e w_zmass2overr= %.4e \n' % (str(sys.argv[1])[0:len(str(sys.argv[1]))-4], centersx[i][j], centersy[i][j], B1608_tab_sim_gal[i][j][k], B1608_tab_sim_oneoverr[i][j][k], B1608_tab_sim_zweight[i][j][k], B1608_tab_sim_mass[i][j][k], B1608_tab_sim_mass2[i][j][k], B1608_tab_sim_mass2rms[i][j][k], B1608_tab_sim_mass3[i][j][k], B1608_tab_sim_mass3rms[i][j][k], B1608_tab_sim_zoverr[i][j][k], B1608_tab_sim_massoverr[i][j][k], B1608_tab_sim_mass2overr[i][j][k], B1608_tab_sim_mass3overr[i][j][k], B1608_tab_sim_mass2overrrms[i][j][k], B1608_tab_sim_mass3overrrms[i][j][k], B1608_tab_sim_zmassoverr[i][j][k], B1608_tab_sim_zmass2overr[i][j][k]))
s.close()
f.close()
#if str(sys.argv[8]) == "samp":
g.close()
#if str(sys.argv[8]) == "tab":
h.close()


print("Total time forB1608_sim: --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
      