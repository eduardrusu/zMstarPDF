# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# run from the lephare_dev/test folder as: python /Users/perseus/Dropbox/Davis_work/code/MstarMillenium.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_8_7_7_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.pdz number
# where number should be different for different processes run in the same time

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

lepharenr=500 # how many objects should lephare be run with
u_millenium=np.zeros(lepharenr)
uerr_millenium=np.zeros(lepharenr)
g_millenium=np.zeros(lepharenr)
gerr_millenium=np.zeros(lepharenr)
r_millenium=np.zeros(lepharenr)
rerr_millenium=np.zeros(lepharenr)
i_millenium=np.zeros(lepharenr)
ierr_millenium=np.zeros(lepharenr)
z_millenium=np.zeros(lepharenr)
zerr_millenium=np.zeros(lepharenr)
J_millenium=np.zeros(lepharenr)
Jerr_millenium=np.zeros(lepharenr)
H_millenium=np.zeros(lepharenr)
Herr_millenium=np.zeros(lepharenr)
Ks_millenium=np.zeros(lepharenr)
Kserr_millenium=np.zeros(lepharenr)
zbest=np.zeros(lepharenr)
zspec=np.zeros(lepharenr)
pofz=np.zeros((lepharenr,70))
mass_best=np.zeros((lepharenr,72)) #because I also compute for z_best,z_spec
mass_inf=np.zeros((lepharenr,72))
mass_med=np.zeros((lepharenr,72))
mass_sup=np.zeros((lepharenr,72))
z=np.linspace(0.05,3.5,70)
os.system("rm %s_mstar_noJHKs.cat" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4]) # since the code only appends, if we have an incomplete previous output we should remove it
os.system("rm %s_mstar_withJHKs.cat" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4])
itrue=0
index=0
name=[]
with open(str(sys.argv[1])) as fields:
    for gal in fields:
        index=index+1
        if gal!="\n": # careful to include this, otherwise the objects at the end of file fail to be included
            if itrue==lepharenr:
                name_in="/Users/perseus/lephare_dev/test/millenium_%s.cat" % str(sys.argv[2])
                name_out="/Users/perseus/lephare_dev/test/millenium_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
                for i in range(1): #with and without JHKs
                    lephare_in = open(name_in,'w')
                    lephare_in.write("# \t ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t Ks \t Ks_err \t context \t z-spec \t string \n")
                    list=[]
                    for k in range(lepharenr): #create list of lists
                        list.append([])
                    for k in range(lepharenr):
                        if i==0: #for zbest
                            lephare_in.write('1 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zbest[k]))
                        if i==1:
                            lephare_in.write('1 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zbest[k]))
                        if i==0: #for zspec
                            lephare_in.write('2 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zspec[k]))
                        if i==1:
                            lephare_in.write('2 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zspec[k]))
                        for j in range(70):
                            if pofz[k][j]>0.001:
                               list[k].append(j)
                               if i==0:
                                   lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (len(list[k])+2,u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],z[j]))
                               if i==1:
                                   lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (len(list[k])+2,u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],z[j]))
                    lephare_in.close()
                    os.system("/Users/perseus/lephare_dev/test/runlephare_phys_para_millenium.sh %s" % name_in)
                    l=0
                    with open('%s' % name_out) as lephare_out:
                        lineout=lephare_out.readlines()
                        for k in range(lepharenr):
                            mass_best[k][0]=float(lineout[62+l].split()[33]) # these first 2 are for zbest and zspec
                            mass_inf[k][0]=float(lineout[62+l].split()[34])
                            mass_med[k][0]=float(lineout[62+l].split()[35])
                            mass_sup[k][0]=float(lineout[62+l].split()[36])
                            mass_best[k][1]=float(lineout[63+l].split()[33])
                            mass_inf[k][1]=float(lineout[63+l].split()[34])
                            mass_med[k][1]=float(lineout[63+l].split()[35])
                            mass_sup[k][1]=float(lineout[63+l].split()[36])
                            l=l+2
                            for j in range(len(list[k])):
                               mass_best[k][list[k][j]+2]=float(lineout[62+l].split()[33])
                               mass_inf[k][list[k][j]+2]=float(lineout[62+l].split()[34])
                               mass_med[k][list[k][j]+2]=float(lineout[62+l].split()[35])
                               mass_sup[k][list[k][j]+2]=float(lineout[62+l].split()[36])
                               l=l+1
                    if i==0:
                        outfile=open('%s_mstar_noJHKs.cat' % str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
                    if i==1:
                        outfile=open('%s_mstar_withJHKs.cat' % str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
                    output=""
                    for k in range(lepharenr):
                        output=output+name[k]+"\t"+str(mass_best[k][0])+"\t"+str(mass_inf[k][0])+"\t"+str(mass_med[k][0])+"\t"+str(mass_sup[k][0])+"\t"+str(mass_best[k][1])+"\t"+str(mass_inf[k][1])+"\t"+str(mass_med[k][1])+"\t"+str(mass_sup[k][1])+"\t"
                        for j in range(70):
                            output=output+str(mass_best[k][j+2])+"\t"+str(mass_inf[k][j+2])+"\t"+str(mass_med[k][j+2])+"\t"+str(mass_sup[k][j+2])+"\t"
                        output=output+"\n"
                    outfile.write(output)
                    outfile.close()
                name=[]
                itrue=0
                mass_best=np.zeros((lepharenr,72))
                mass_inf=np.zeros((lepharenr,72))
                mass_med=np.zeros((lepharenr,72))
                mass_sup=np.zeros((lepharenr,72))
            print "------ running object number: \t", index, "------"
            name.append(gal.split()[0])
            u_millenium[itrue]=float(gal.split()[6])
            uerr_millenium[itrue]=float(gal.split()[7])
            g_millenium[itrue]=float(gal.split()[8])
            gerr_millenium[itrue]=float(gal.split()[9])
            r_millenium[itrue]=float(gal.split()[10])
            rerr_millenium[itrue]=float(gal.split()[11])
            i_millenium[itrue]=float(gal.split()[12])
            ierr_millenium[itrue]=float(gal.split()[13])
            z_millenium[itrue]=float(gal.split()[14])
            zerr_millenium[itrue]=float(gal.split()[15])
            J_millenium[itrue]=float(gal.split()[16])
            Jerr_millenium[itrue]=float(gal.split()[17])
            H_millenium[itrue]=float(gal.split()[18])
            Herr_millenium[itrue]=float(gal.split()[19])
            Ks_millenium[itrue]=float(gal.split()[20])
            Kserr_millenium[itrue]=float(gal.split()[21])
            zbest[itrue]=float(gal.split()[3])
            zspec[itrue]=float(gal.split()[2])
            for j in range(70):     # I MODIFIED FROM 69, BUT NEED TO CHECK IF CORRECT
                #print linepdz[i].split()[j+2], "\n"
                #print [x[0:len(linepdz[i].split()[j+2])-1] for x in linepdz[i].split()[j+2]], "\n"
                pofz[itrue][j]=float(str(gal.split()[j+22]))
            itrue=itrue+1
             #the code below is necessary to deal with the objects at the end of the file, if there are less than lepharenr objects left
name_in="/Users/perseus/lephare_dev/test/millenium_%s.cat" % str(sys.argv[2])
name_out="/Users/perseus/lephare_dev/test/millenium_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
for i in range(1):
    lephare_in = open(name_in,'w')
    lephare_in.write("# \t ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t J \t J_err \t H \t H_err \t Ks \t Ks_err \t context \t z-spec \t string \n")
    list=[]
    for k in range(itrue): #create list of lists
        list.append([])
    for k in range(itrue):
        if i==0: #for zbest
            lephare_in.write('1 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zbest[k]))
        if i==1:
            lephare_in.write('1 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zbest[k]))
        if i==0: #for zspec
            lephare_in.write('2 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zspec[k]))
        if i==1:
            lephare_in.write('2 %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],zspec[k]))
        for j in range(70):
            if pofz[k][j]>0.001:
                list[k].append(j)
                if i==0:
                    lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 31 %s\n' % (len(list[k])+2,u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],z[j]))
                if i==1:
                    lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s 255 %s\n' % (len(list[k])+2,u_millenium[k],uerr_millenium[k],g_millenium[k],gerr_millenium[k],r_millenium[k],rerr_millenium[k],i_millenium[k],ierr_millenium[k],z_millenium[k],zerr_millenium[k],J_millenium[k],Jerr_millenium[k],H_millenium[k],Herr_millenium[k],Ks_millenium[k],Kserr_millenium[k],z[j]))
    lephare_in.close()
    os.system("/Users/perseus/lephare_dev/test/runlephare_phys_para_millenium.sh %s" % name_in)
    l=0
    with open('%s' % name_out) as lephare_out:
        lineout=lephare_out.readlines()
        for k in range(itrue):
            mass_best[k][0]=float(lineout[62+l].split()[33])
            mass_inf[k][0]=float(lineout[62+l].split()[34])
            mass_med[k][0]=float(lineout[62+l].split()[35])
            mass_sup[k][0]=float(lineout[62+l].split()[36])
            mass_best[k][1]=float(lineout[63+l].split()[33])
            mass_inf[k][1]=float(lineout[63+l].split()[34])
            mass_med[k][1]=float(lineout[63+l].split()[35])
            mass_sup[k][1]=float(lineout[63+l].split()[36])
            l=l+2
            for j in range(len(list[k])):
                mass_best[k][list[k][j]+2]=float(lineout[62+l].split()[33])
                mass_inf[k][list[k][j]+2]=float(lineout[62+l].split()[34])
                mass_med[k][list[k][j]+2]=float(lineout[62+l].split()[35])
                mass_sup[k][list[k][j]+2]=float(lineout[62+l].split()[36])
                l=l+1
    if i==0:
        outfile=open('%s_mstar_noJHKs.cat' % str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
    if i==1:
        outfile=open('%s_mstar_withJHKs.cat' % str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
    output=""
    for k in range(itrue):
        output=output+name[k]+"\t"+str(mass_best[k][0])+"\t"+str(mass_inf[k][0])+"\t"+str(mass_med[k][0])+"\t"+str(mass_sup[k][0])+"\t"+str(mass_best[k][1])+"\t"+str(mass_inf[k][1])+"\t"+str(mass_med[k][1])+"\t"+str(mass_sup[k][1])+"\t"
        for j in range(70):
            output=output+str(mass_best[k][j+2])+"\t"+str(mass_inf[k][j+2])+"\t"+str(mass_med[k][j+2])+"\t"+str(mass_sup[k][j+2])+"\t"
        output=output+"\n"
    outfile.write(output)
    outfile.close()

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))
                               
print 'Done!'
