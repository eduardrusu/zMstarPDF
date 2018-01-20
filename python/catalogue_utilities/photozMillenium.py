# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# run from the bpz/test folder as: python /Users/perseus/Dropbox/Davis_work/code/photozMillenium.py #complete address to the simulation     field number
# where the address is, e.g. /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_8_7_6_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt, and number is an integer, corresponding to the order number of each run of this code. Run the code from inside bpz/test

import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time

start_timefield = time.time()

bpznr=1000 # how many objects should bpz be run with
u_millenium=np.zeros(bpznr)
uerr_millenium=np.zeros(bpznr)
g_millenium=np.zeros(bpznr)
gerr_millenium=np.zeros(bpznr)
r_millenium=np.zeros(bpznr)
rerr_millenium=np.zeros(bpznr)
i_millenium=np.zeros(bpznr)
ierr_millenium=np.zeros(bpznr)
z_millenium=np.zeros(bpznr)
zerr_millenium=np.zeros(bpznr)
J_millenium=np.zeros(bpznr)
Jerr_millenium=np.zeros(bpznr)
H_millenium=np.zeros(bpznr)
Herr_millenium=np.zeros(bpznr)
K_millenium=np.zeros(bpznr)
Kerr_millenium=np.zeros(bpznr)
id=np.zeros(bpznr)
pofz=np.zeros((bpznr,70))
mstar=np.zeros(bpznr)
posx=np.zeros(bpznr)
posy=np.zeros(bpznr)
z_best=np.zeros(bpznr)
z_spec=np.zeros(bpznr)
#z=np.linspace(0.05,3.5,70)
os.system("cp millenium.columns millenium_%s.columns" % str(sys.argv[2]))
name_in="/Users/perseus/bpz-1.99.3/test/millenium_%s.cat" % str(sys.argv[2])
name_outzbest="/Users/perseus/bpz-1.99.3/test/millenium_%s_bpz.cat" % str(sys.argv[2])
name_outpdz="/Users/perseus/bpz-1.99.3/test/millenium_%s.probs" % str(sys.argv[2])
name_out="%s.pdz" % (str(sys.argv[1])[0:len(str(sys.argv[1]))-11])
os.system("rm %s" % name_out) # since the code only appends, if we have an incomplete previous output we should remove it
itrue=0 # index of only the objects passing the selection criteria, between 0 and bpznr-1
with open(str(sys.argv[1])) as fields:
    for gal in fields:
      if gal!="\n": # careful to include this, otherwise the objects at the end of file fail to be included
        if gal.split()[0]!="GalID":
            if itrue==bpznr:
                itrue=0
                bpz_in = open(name_in,'w')
                bpz_in.write("#ID u       u_err g       g_err r       r_err i       i_err z       z_err J       J_err H       H_err K       K_err \n")
                for i in range(bpznr):
                    bpz_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' % (i+1,u_millenium[i],uerr_millenium[i],g_millenium[i],gerr_millenium[i],r_millenium[i],rerr_millenium[i],i_millenium[i],ierr_millenium[i],z_millenium[i],zerr_millenium[i],J_millenium[i],Jerr_millenium[i],H_millenium[i],Herr_millenium[i],K_millenium[i],Kerr_millenium[i]))
                bpz_in.close()
                os.system("python $BPZPATH/bpz.py %s -INTERP 2" % name_in)
                os.system("python $BPZPATH/bpzfinalize.py %s" % name_in[0:len(name_in)-4])
                l=0
                with open('%s' % name_outzbest) as bpz_outzbest:
                    for outzbest in bpz_outzbest:
                        if (outzbest.split()[0]!="#") and (outzbest.split()[0]!="#id"):
                            #print outzbest.split()[0], "\n"
                            z_best[l]=outzbest.split()[1]
                            l=l+1
                l=0
                with open('%s' % name_outpdz) as bpz_outpdz:
                    for outpdz in bpz_outpdz:
                        if outpdz.split()[0]!="#":
                            for i in range(70):
                                pofz[l][i]=float(outpdz.split()[1+i])
                            l=l+1
                outfile=open(name_out,'a')
                output=""
                for i in range(bpznr):
                    strid="%15d" % id[i]
                    output=output+strid+"\t"+str(mstar[i])+"\t"+str(z_spec[i])+"\t"+str(z_best[i])+"\t"+str(posx[i])+"\t"+str(posy[i])+"\t"+str(u_millenium[i])+"\t"+str(uerr_millenium[i])+"\t"+str(g_millenium[i])+"\t"+str(gerr_millenium[i])+"\t"+str(r_millenium[i])+"\t"+str(rerr_millenium[i])+"\t"+str(i_millenium[i])+"\t"+str(ierr_millenium[i])+"\t"+str(z_millenium[i])+"\t"+str(zerr_millenium[i])+"\t"+str(J_millenium[i])+"\t"+str(Jerr_millenium[i])+"\t"+str(H_millenium[i])+"\t"+str(Herr_millenium[i])+"\t"+str(K_millenium[i])+"\t"+str(Kerr_millenium[i])+"\t"
                    for j in range(70):
                        output=output+str(pofz[i][j])+"\t"
                    output=output+"\n"
                outfile.write(output)
                outfile.close()
                #print gal.split()[15], "\n"
            if float(gal.split()[15]) <= 24: # if mag_i < 24
                id[itrue]=float(gal.split()[0])
                mstar[itrue]=float(gal.split()[11])
                z_spec[itrue]=float(gal.split()[5])
                posx[itrue]=float(gal.split()[6])
                posy[itrue]=float(gal.split()[7])
                u_millenium[itrue]=float(gal.split()[12])
                if u_millenium[itrue] < 23:
                    uerr_millenium[itrue]=0.01
                if (u_millenium[itrue] < 24) and (u_millenium[itrue] > 23):
                    uerr_millenium[itrue]=0.02
                if (u_millenium[itrue] < 24.5) and (u_millenium[itrue] > 24):
                    uerr_millenium[itrue]=0.03
                if (u_millenium[itrue] < 24.7) and (u_millenium[itrue] > 24.5):
                    uerr_millenium[itrue]=0.04
                if (u_millenium[itrue] < 24.9) and (u_millenium[itrue] > 24.7):
                    uerr_millenium[itrue]=0.05
                if (u_millenium[itrue] < 25.1) and (u_millenium[itrue] > 24.9):
                    uerr_millenium[itrue]=0.06
                if (u_millenium[itrue] < 25.2) and (u_millenium[itrue] > 25.1):
                    uerr_millenium[itrue]=0.07
                if (u_millenium[itrue] < 25.3) and (u_millenium[itrue] > 25.2):
                    uerr_millenium[itrue]=0.08
                if (u_millenium[itrue] < 25.4) and (u_millenium[itrue] > 25.3):
                    uerr_millenium[itrue]=0.09
                if (u_millenium[itrue] < 25.5) and (u_millenium[itrue] > 25.4):
                    uerr_millenium[itrue]=0.10
                if (u_millenium[itrue] < 26) and (u_millenium[itrue] > 25.5):
                    uerr_millenium[itrue]=0.13
                if (u_millenium[itrue] < 26.5) and (u_millenium[itrue] > 26):
                    uerr_millenium[itrue]=0.18
                if (u_millenium[itrue] < 27) and (u_millenium[itrue] > 26.5):
                    uerr_millenium[itrue]=0.25
                if u_millenium[itrue] > 24:
                    uerr_millenium[itrue]=0.30
                g_millenium[itrue]=float(gal.split()[13])
                if g_millenium[itrue] < 24:
                    gerr_millenium[itrue]=0.01
                if (g_millenium[itrue] < 24.5) and (g_millenium[itrue] > 24):
                    gerr_millenium[itrue]=0.02
                if (g_millenium[itrue] < 24.8) and (g_millenium[itrue] > 24.5):
                    gerr_millenium[itrue]=0.03
                if (g_millenium[itrue] < 25) and (g_millenium[itrue] > 24.8):
                    gerr_millenium[itrue]=0.04
                if (g_millenium[itrue] < 25.3) and (g_millenium[itrue] > 25):
                    gerr_millenium[itrue]=0.05
                if (g_millenium[itrue] < 25.5) and (g_millenium[itrue] > 25.3):
                    gerr_millenium[itrue]=0.06
                if (g_millenium[itrue] < 25.7) and (g_millenium[itrue] > 25.5):
                    gerr_millenium[itrue]=0.07
                if (g_millenium[itrue] < 25.8) and (g_millenium[itrue] > 25.7):
                    gerr_millenium[itrue]=0.08
                if (g_millenium[itrue] < 26) and (g_millenium[itrue] > 25.8):
                    gerr_millenium[itrue]=0.09
                if (g_millenium[itrue] < 26.2) and (g_millenium[itrue] > 26):
                    gerr_millenium[itrue]=0.10
                if (g_millenium[itrue] < 26.5) and (g_millenium[itrue] > 26.2):
                    gerr_millenium[itrue]=0.13
                if (g_millenium[itrue] < 27) and (g_millenium[itrue] > 26.5):
                    gerr_millenium[itrue]=0.18
                if g_millenium[itrue] > 27:
                    gerr_millenium[itrue]=0.18
                r_millenium[itrue]=float(gal.split()[14])
                if r_millenium[itrue] < 23.2:
                    rerr_millenium[itrue]=0.01
                if (r_millenium[itrue] < 24.1) and (r_millenium[itrue] > 23.2):
                    rerr_millenium[itrue]=0.02
                if (r_millenium[itrue] < 24.5) and (r_millenium[itrue] > 24.1):
                    rerr_millenium[itrue]=0.03
                if (r_millenium[itrue] < 24.9) and (r_millenium[itrue] > 24.5):
                    rerr_millenium[itrue]=0.04
                if (r_millenium[itrue] < 25.0) and (r_millenium[itrue] > 24.9):
                    rerr_millenium[itrue]=0.05
                if (r_millenium[itrue] < 25.2) and (r_millenium[itrue] > 25.0):
                    rerr_millenium[itrue]=0.06
                if (r_millenium[itrue] < 25.3) and (r_millenium[itrue] > 25.2):
                    rerr_millenium[itrue]=0.07
                if (r_millenium[itrue] < 25.5) and (r_millenium[itrue] > 25.3):
                    rerr_millenium[itrue]=0.08
                if (r_millenium[itrue] < 25.6) and (r_millenium[itrue] > 25.5):
                    rerr_millenium[itrue]=0.09
                if (r_millenium[itrue] < 25.7) and (r_millenium[itrue] > 25.6):
                    rerr_millenium[itrue]=0.10
                if (r_millenium[itrue] < 26.0) and (r_millenium[itrue] > 25.7):
                    rerr_millenium[itrue]=0.12
                if r_millenium[itrue] > 26:
                    rerr_millenium[itrue]=0.12
                i_millenium[itrue]=float(gal.split()[15])
                if i_millenium[itrue] < 22.8:
                    ierr_millenium[itrue]=0.01
                if (i_millenium[itrue] < 23.6) and (i_millenium[itrue] > 22.8):
                    ierr_millenium[itrue]=0.02
                if (i_millenium[itrue] < 24.0) and (i_millenium[itrue] > 23.6):
                    ierr_millenium[itrue]=0.03
                if (i_millenium[itrue] < 24.2) and (i_millenium[itrue] > 24.0):
                    ierr_millenium[itrue]=0.04
                if (i_millenium[itrue] < 24.4) and (i_millenium[itrue] > 24.2):
                    ierr_millenium[itrue]=0.05
                if (i_millenium[itrue] < 24.6) and (i_millenium[itrue] > 24.4):
                    ierr_millenium[itrue]=0.06
                if (i_millenium[itrue] < 24.7) and (i_millenium[itrue] > 24.6):
                    ierr_millenium[itrue]=0.07
                if (i_millenium[itrue] < 24.8) and (i_millenium[itrue] > 24.7):
                    ierr_millenium[itrue]=0.08
                if (i_millenium[itrue] < 24.9) and (i_millenium[itrue] > 24.8):
                    ierr_millenium[itrue]=0.09
                if (i_millenium[itrue] < 25.0) and (i_millenium[itrue] > 24.9):
                    ierr_millenium[itrue]=0.10
                if (i_millenium[itrue] < 25.2) and (i_millenium[itrue] > 25.0):
                    ierr_millenium[itrue]=0.11
                if i_millenium[itrue] > 25.2:
                    ierr_millenium[itrue]=0.11
                z_millenium[itrue]=float(gal.split()[16])
                if z_millenium[itrue] < 22.0:
                    zerr_millenium[itrue]=0.01
                if (z_millenium[itrue] < 22.8) and (z_millenium[itrue] > 22.0):
                    zerr_millenium[itrue]=0.02
                if (z_millenium[itrue] < 23.2) and (z_millenium[itrue] > 22.8):
                    zerr_millenium[itrue]=0.03
                if (z_millenium[itrue] < 24.2) and (z_millenium[itrue] > 23.2):
                    zerr_millenium[itrue]=0.04
                if (z_millenium[itrue] < 24.4) and (z_millenium[itrue] > 24.2):
                    zerr_millenium[itrue]=0.05
                if (z_millenium[itrue] < 24.6) and (z_millenium[itrue] > 24.4):
                    zerr_millenium[itrue]=0.06
                if (z_millenium[itrue] < 24.7) and (z_millenium[itrue] > 24.6):
                    zerr_millenium[itrue]=0.07
                if (z_millenium[itrue] < 24.8) and (z_millenium[itrue] > 24.7):
                    zerr_millenium[itrue]=0.08
                if (z_millenium[itrue] < 24.9) and (z_millenium[itrue] > 24.8):
                    zerr_millenium[itrue]=0.09
                if (z_millenium[itrue] < 25.0) and (z_millenium[itrue] > 24.9):
                    zerr_millenium[itrue]=0.10
                if (z_millenium[itrue] < 25.2) and (z_millenium[itrue] > 25.0):
                    zerr_millenium[itrue]=0.11
                if z_millenium[itrue] > 25.2:
                    zerr_millenium[itrue]=0.11
                J_millenium[itrue]=float(gal.split()[17])
                if J_millenium[itrue] < 20.5:
                    Jerr_millenium[itrue]=0.01
                if (J_millenium[itrue] < 21.2) and (J_millenium[itrue] > 20.5):
                    Jerr_millenium[itrue]=0.02
                if (J_millenium[itrue] < 22.0) and (J_millenium[itrue] > 21.2):
                    Jerr_millenium[itrue]=0.03
                if (J_millenium[itrue] < 22.4) and (J_millenium[itrue] > 22.0):
                    Jerr_millenium[itrue]=0.04
                if (J_millenium[itrue] < 22.6) and (J_millenium[itrue] > 22.4):
                    Jerr_millenium[itrue]=0.05
                if (J_millenium[itrue] < 22.8) and (J_millenium[itrue] > 22.6):
                    Jerr_millenium[itrue]=0.06
                if (J_millenium[itrue] < 22.9) and (J_millenium[itrue] > 22.8):
                    Jerr_millenium[itrue]=0.07
                if (J_millenium[itrue] < 23.1) and (J_millenium[itrue] > 22.9):
                    Jerr_millenium[itrue]=0.08
                if (J_millenium[itrue] < 23.2) and (J_millenium[itrue] > 23.1):
                    Jerr_millenium[itrue]=0.09
                if (J_millenium[itrue] < 23.4) and (J_millenium[itrue] > 23.2):
                    Jerr_millenium[itrue]=0.10
                if (J_millenium[itrue] < 23.7) and (J_millenium[itrue] > 23.4):
                    Jerr_millenium[itrue]=0.15
                if (J_millenium[itrue] < 24.0) and (J_millenium[itrue] > 23.7):
                    Jerr_millenium[itrue]=0.18
                if J_millenium[itrue] > 24.0:
                    Jerr_millenium[itrue]=0.18
                H_millenium[itrue]=float(gal.split()[18]) # I don't have H data yet so I use an average between J and K
                if H_millenium[itrue] < 21.9:
                    Herr_millenium[itrue]=0.01
                if (H_millenium[itrue] < 22.6) and (H_millenium[itrue] > 21.9):
                    Herr_millenium[itrue]=0.02
                if (H_millenium[itrue] < 22.9) and (H_millenium[itrue] > 22.6):
                    Herr_millenium[itrue]=0.03
                if (H_millenium[itrue] < 23.1) and (H_millenium[itrue] > 22.9):
                    Herr_millenium[itrue]=0.04
                if (H_millenium[itrue] < 23.3) and (H_millenium[itrue] > 23.1):
                    Herr_millenium[itrue]=0.05
                if (H_millenium[itrue] < 23.4) and (H_millenium[itrue] > 23.3):
                    Herr_millenium[itrue]=0.06
                if (H_millenium[itrue] < 23.6) and (H_millenium[itrue] > 23.4):
                    Herr_millenium[itrue]=0.07
                if (H_millenium[itrue] < 23.8) and (H_millenium[itrue] > 23.6):
                    Herr_millenium[itrue]=0.08
                if (H_millenium[itrue] < 24.0) and (H_millenium[itrue] > 23.8):
                    Herr_millenium[itrue]=0.09
                if (H_millenium[itrue] < 24.1) and (H_millenium[itrue] > 24.0):
                    Herr_millenium[itrue]=0.10
                if (H_millenium[itrue] < 24.5) and (H_millenium[itrue] > 24.1):
                    Herr_millenium[itrue]=0.15
                if H_millenium[itrue] > 24.5:
                    Herr_millenium[itrue]=0.15
                K_millenium[itrue]=float(gal.split()[19])
                if K_millenium[itrue] < 22.0:
                    Kerr_millenium[itrue]=0.01
                if (K_millenium[itrue] < 22.8) and (K_millenium[itrue] > 22.0):
                    Kerr_millenium[itrue]=0.02
                if (K_millenium[itrue] < 23.3) and (K_millenium[itrue] > 22.8):
                    Kerr_millenium[itrue]=0.03
                if (K_millenium[itrue] < 23.5) and (K_millenium[itrue] > 23.3):
                    Kerr_millenium[itrue]=0.04
                if (K_millenium[itrue] < 23.7) and (K_millenium[itrue] > 23.5):
                    Kerr_millenium[itrue]=0.05
                if (K_millenium[itrue] < 23.9) and (K_millenium[itrue] > 23.7):
                    Kerr_millenium[itrue]=0.06
                if (K_millenium[itrue] < 24.2) and (K_millenium[itrue] > 23.9):
                    Kerr_millenium[itrue]=0.07
                if (K_millenium[itrue] < 24.4) and (K_millenium[itrue] > 24.2):
                    Kerr_millenium[itrue]=0.08
                if (K_millenium[itrue] < 24.5) and (K_millenium[itrue] > 24.4):
                    Kerr_millenium[itrue]=0.09
                if (K_millenium[itrue] < 24.6) and (K_millenium[itrue] > 24.5):
                    Kerr_millenium[itrue]=0.10
                if (K_millenium[itrue] < 25.0) and (K_millenium[itrue] > 24.6):
                    Kerr_millenium[itrue]=0.15
                if K_millenium[itrue] > 24.0:
                    Kerr_millenium[itrue]=0.15
                itrue=itrue+1

    #the code below is necessary to deal with the objects at the end of the file, if there are less than bpznr objects left
bpz_in = open(name_in,'w')
bpz_in.write("#ID u       u_err g       g_err r       r_err i       i_err z       z_err J       J_err H       H_err K       K_err \n")
for i in range(itrue):
    bpz_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' % (i+1,u_millenium[i],uerr_millenium[i],g_millenium[i],gerr_millenium[i],r_millenium[i],rerr_millenium[i],i_millenium[i],ierr_millenium[i],z_millenium[i],zerr_millenium[i],J_millenium[i],Jerr_millenium[i],H_millenium[i],Herr_millenium[i],K_millenium[i],Kerr_millenium[i]))
bpz_in.close()
os.system("python $BPZPATH/bpz.py %s -INTERP 2" % name_in)
os.system("python $BPZPATH/bpzfinalize.py %s" % name_in[0:len(name_in)-4])
l=0
with open('%s' % name_outzbest) as bpz_outzbest:
    for outzbest in bpz_outzbest:
        if (outzbest.split()[0]!="#") and (outzbest.split()[0]!="#id"):
            z_best[l]=outzbest.split()[1]
            l=l+1
l=0
with open('%s' % name_outpdz) as bpz_outpdz:
    for outpdz in bpz_outpdz:
        if outpdz.split()[0]!="#":
            for i in range(70):
                pofz[l][i]=float(outpdz.split()[1+i])
            l=l+1
outfile=open(name_out,'a')
output=""
for i in range(itrue):
    output=output+str(id[i])+"\t"+str(mstar[i])+"\t"+str(z_spec[i])+"\t"+str(z_best[i])+"\t"+str(posx[i])+"\t"+str(posy[i])+"\t"+str(u_millenium[i])+"\t"+str(uerr_millenium[i])+"\t"+str(g_millenium[i])+"\t"+str(gerr_millenium[i])+"\t"+str(r_millenium[i])+"\t"+str(rerr_millenium[i])+"\t"+str(i_millenium[i])+"\t"+str(ierr_millenium[i])+"\t"+str(z_millenium[i])+"\t"+str(zerr_millenium[i])+"\t"+str(J_millenium[i])+"\t"+str(Jerr_millenium[i])+"\t"+str(H_millenium[i])+"\t"+str(Herr_millenium[i])+"\t"+str(K_millenium[i])+"\t"+str(Kerr_millenium[i])+"\t"
    for j in range(70):
        output=output+str(pofz[i][j])+"\t"
    output=output+"\n"
outfile.write(output)
outfile.close()

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

