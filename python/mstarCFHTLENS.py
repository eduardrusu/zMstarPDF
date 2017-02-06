# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# run from the lephare_dev/test folder as: python /Users/perseus/Dropbox/Davis_work/code/mstarCFHTLENS.py /Users/perseus/Dropbox/Davis_work/code/fieldstry.lst number
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

with open(str(sys.argv[1])) as fields:  # fieldstry#.lst
    listfields = fields.readlines()
lepharenr=10 # how many objects should lephare be run with
u_CFHTLS=np.zeros(lepharenr)
uerr_CFHTLS=np.zeros(lepharenr)
g_CFHTLS=np.zeros(lepharenr)
gerr_CFHTLS=np.zeros(lepharenr)
r_CFHTLS=np.zeros(lepharenr)
rerr_CFHTLS=np.zeros(lepharenr)
i_CFHTLS=np.zeros(lepharenr)
ierr_CFHTLS=np.zeros(lepharenr)
y_CFHTLS=np.zeros(lepharenr)
yerr_CFHTLS=np.zeros(lepharenr)
z_CFHTLS=np.zeros(lepharenr)
zerr_CFHTLS=np.zeros(lepharenr)
pofz=np.zeros((lepharenr,70))
mass_best=np.zeros((lepharenr,70))
mass_inf=np.zeros((lepharenr,70))
mass_med=np.zeros((lepharenr,70))
mass_sup=np.zeros((lepharenr,70))
chigal=np.zeros((lepharenr,70))
chistar=np.zeros((lepharenr,70))
z=np.linspace(0.05,3.5,70)
for count in range(len(listfields)):
    start_timesubfield = time.time()
    i=0 # object index in original catalogue
    itrue=0 # index of only the objects passing the selection criteria
    name=[]
    with open('%s.cat' % [x[0:len(listfields[0])-1] for x in listfields][count]) as orig:
      os.system("rm %s_mstar.cat" % [x[0:len(listfields[0])-1] for x in listfields][count]) # since the code only appends, if we have an incomplete previous output we should remove it
      with open('%s_pdz.cat' % [x[0:len(listfields[0])-1] for x in listfields][count]) as pdz:
        linepdz=pdz.readlines()
        for lineorig in orig:
          if lineorig.split()[0]=="#id":
              check="false"
          else:
              #print lineorig.split()[0], linepdz[i].split()[0], "\n"
            if lineorig.split()[0] != linepdz[i].split()[0]:
                print "Catalogue mismatch!"
                break
            check="false"
            if itrue==lepharenr:
                if i_CFHTLS[0]>0: # the code does not work properly for the < lepharenr objects with y but for which the .sh script considers the mag prior to be in i. Those need to be rerun separately. Actually this is not an issue since it appears that the CFHTLENS subfields do not contain mixed i and y mag objects. 
                    name_in="/Users/perseus/lephare_dev/test/izrgu_%s.cat" % str(sys.argv[2])
                    name_out="/Users/perseus/lephare_dev/test/izrgu_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
                else:
                    name_in="/Users/perseus/lephare_dev/test/yzrgu_%s.cat" % str(sys.argv[2])
                    name_out="/Users/perseus/lephare_dev/test/yzrgu_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
                lephare_in = open(name_in,'w')
                lephare_in.write("# \t ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t y \t y_err \t z \t z_err \t context \t z-spec \t string \n")
                list=[]
                for k in range(lepharenr): #create list of lists
                    list.append([])
                for k in range(lepharenr):
                    for j in range(70):
                        if pofz[k][j]>0.001:
                            list[k].append(j)
                            if i_CFHTLS[0]>0:
                                lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s 47 %s\n' % (len(list[k]),u_CFHTLS[k],uerr_CFHTLS[k],g_CFHTLS[k],gerr_CFHTLS[k],r_CFHTLS[k],rerr_CFHTLS[k],i_CFHTLS[k],ierr_CFHTLS[k],y_CFHTLS[k],yerr_CFHTLS[k],z_CFHTLS[k],zerr_CFHTLS[k],z[j]))
                            else:
                                lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s 55 %s\n' % (len(list[k]),u_CFHTLS[k],uerr_CFHTLS[k],g_CFHTLS[k],gerr_CFHTLS[k],r_CFHTLS[k],rerr_CFHTLS[k],i_CFHTLS[k],ierr_CFHTLS[k],y_CFHTLS[k],yerr_CFHTLS[k],z_CFHTLS[k],zerr_CFHTLS[k],z[j]))
                lephare_in.close()
                os.system("/Users/perseus/lephare_dev/test/runlephare_phys_para_CFHTLS.sh %s" % name_in)
                l=0
                with open('%s' % name_out) as lephare_out:
                    lineout=lephare_out.readlines()
                    for k in range(lepharenr):
                        for j in range(len(list[k])):
                            mass_best[k][list[k][j]]=float(lineout[62+l].split()[29])
                            mass_inf[k][list[k][j]]=float(lineout[62+l].split()[30])
                            mass_med[k][list[k][j]]=float(lineout[62+l].split()[31])
                            mass_sup[k][list[k][j]]=float(lineout[62+l].split()[32])
                            chigal[k][list[k][j]]=float(lineout[62+l].split()[5])
                            chistar[k][list[k][j]]=float(lineout[62+l].split()[6])
                            l=l+1
                outfile=open('%s_mstar.cat' % [x[0:len(listfields[0])-1] for x in listfields][count],'a')
                output=""
                for k in range(lepharenr):
                    output=output+name[k]+"\t"
                    for j in range(70):
                        output=output+str(mass_best[k][j])+"\t"+str(mass_inf[k][j])+"\t"+str(mass_med[k][j])+"\t"+str(mass_sup[k][j])+"\t"+str(chigal[k][j])+"\t"+str(chistar[k][j])+"\t"
                    output=output+"\n"
                outfile.write(output)
                outfile.close()
                name=[]
                itrue=0
                mass_best=np.zeros((lepharenr,70))
                mass_inf=np.zeros((lepharenr,70))
                mass_med=np.zeros((lepharenr,70))
                mass_sup=np.zeros((lepharenr,70))
                chigal=np.zeros((lepharenr,70))
                chistar=np.zeros((lepharenr,70))
            #print lineorig.split()[60], "\n"
            if (int(lineorig.split()[60]) == 0) and ((float(lineorig.split()[79]) <= 24 and float(lineorig.split()[79]) > 0) or (float(lineorig.split()[84]) <= 24 and float(lineorig.split()[84]) > 0)): #star_flag and mag_i, mag_y
                print "------ running object number: \t", i, "------"
                name.append(lineorig.split()[0])
                u_CFHTLS[itrue]=float(lineorig.split()[64])
                uerr_CFHTLS[itrue]=float(lineorig.split()[65])
                g_CFHTLS[itrue]=float(lineorig.split()[69])
                gerr_CFHTLS[itrue]=float(lineorig.split()[70])
                r_CFHTLS[itrue]=float(lineorig.split()[74])
                rerr_CFHTLS[itrue]=float(lineorig.split()[75])
                i_CFHTLS[itrue]=float(lineorig.split()[79])
                ierr_CFHTLS[itrue]=float(lineorig.split()[80])
                y_CFHTLS[itrue]=float(lineorig.split()[84])
                yerr_CFHTLS[itrue]=float(lineorig.split()[85])
                z_CFHTLS[itrue]=float(lineorig.split()[89])
                zerr_CFHTLS[itrue]=float(lineorig.split()[90])
                if u_CFHTLS[itrue]==99.0:
                    u_CFHTLS[itrue]=-99.0
                if uerr_CFHTLS[itrue]==99.0:
                    uerr_CFHTLS[itrue]=-99.0
                if g_CFHTLS[itrue]==99.0:
                    g_CFHTLS[itrue]=-99.0
                if gerr_CFHTLS[itrue]==99.0:
                    gerr_CFHTLS[itrue]=-99.0
                if r_CFHTLS[itrue]==99.0:
                    r_CFHTLS[itrue]=-99.0
                if rerr_CFHTLS[itrue]==99.0:
                    rerr_CFHTLS[itrue]=-99.0
                if i_CFHTLS[itrue]==99.0:
                    i_CFHTLS[itrue]=-99.0
                if ierr_CFHTLS[itrue]==99.0:
                    ierr_CFHTLS[itrue]=-99.0
                if y_CFHTLS[itrue]==99.0:
                    y_CFHTLS[itrue]=-99.0
                if yerr_CFHTLS[itrue]==99.0:
                    yerr_CFHTLS[itrue]=-99.0
                if z_CFHTLS[itrue]==99.0:
                    z_CFHTLS[itrue]=-99.0
                if zerr_CFHTLS[itrue]==99.0:
                    zerr_CFHTLS[itrue]=-99.0
                for j in range(69):
                    #print linepdz[i].split()[j+2], "\n"
                    #print [x[0:len(linepdz[i].split()[j+2])-1] for x in linepdz[i].split()[j+2]], "\n"
                    string=str(linepdz[i].split()[j+2])
                    pofz[itrue][j]=float(string[:-1]) #because the last character is ","
                string=str(linepdz[i].split()[71])
                pofz[itrue][69]=float(string)
                i=i+1
                itrue=itrue+1
                check="true"
          if check=="false":
              i=i+1
    #the code below is necessary to deal with the objects at the end of the file, if there are less than lepharenr objects left
    if i_CFHTLS[0]>0:  # this is a drawback of the current code, in case the lephare input contains mixed i,y all entries are considered as i
        name_in="/Users/perseus/lephare_dev/test/izrgu_%s.cat" % str(sys.argv[2])
        name_out="/Users/perseus/lephare_dev/test/izrgu_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
    else:
        name_in="/Users/perseus/lephare_dev/test/yzrgu_%s.cat" % str(sys.argv[2])
        name_out="/Users/perseus/lephare_dev/test/yzrgu_%s.cat.MAG_BC03_I09.lephareout" % str(sys.argv[2])
    lephare_in = open(name_in,'w')
    lephare_in.write("# \t ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t y \t y_err \t z \t z_err \t context \t z-spec \t string \n")
    list=[]
    for k in range(itrue): #create list of lists
        list.append([])
    for k in range(itrue):
        for j in range(70):
            if pofz[k][j]>0.001:
                list[k].append(j)
                if i_CFHTLS[0]>0:
                    lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s 47 %s\n' % (len(list[k]),u_CFHTLS[k],uerr_CFHTLS[k],g_CFHTLS[k],gerr_CFHTLS[k],r_CFHTLS[k],rerr_CFHTLS[k],i_CFHTLS[k],ierr_CFHTLS[k],y_CFHTLS[k],yerr_CFHTLS[k],z_CFHTLS[k],zerr_CFHTLS[k],z[j]))
                else:
                    lephare_in.write('%s %s %s %s %s %s %s %s %s %s %s %s %s 55 %s\n' % (len(list[k]),u_CFHTLS[k],uerr_CFHTLS[k],g_CFHTLS[k],gerr_CFHTLS[k],r_CFHTLS[k],rerr_CFHTLS[k],i_CFHTLS[k],ierr_CFHTLS[k],y_CFHTLS[k],yerr_CFHTLS[k],z_CFHTLS[k],zerr_CFHTLS[k],z[j]))
    lephare_in.close()
    os.system("/Users/perseus/lephare_dev/test/runlephare_phys_para.sh %s" % name_in)
    l=0
    with open('%s' % name_out) as lephare_out:
        lineout=lephare_out.readlines()
        for k in range(itrue):
            for j in range(len(list[k])):
                mass_best[k][list[k][j]]=float(lineout[62+l].split()[29])
                mass_inf[k][list[k][j]]=float(lineout[62+l].split()[30])
                mass_med[k][list[k][j]]=float(lineout[62+l].split()[31])
                mass_sup[k][list[k][j]]=float(lineout[62+l].split()[32])
                chigal[k][list[k][j]]=float(lineout[62+l].split()[5])
                chistar[k][list[k][j]]=float(lineout[62+l].split()[6])
                l=l+1
    outfile=open('%s_mstar.cat' % [x[0:len(listfields[0])-1] for x in listfields][count],'a')
    output=""
    for k in range(itrue):
        output=output+name[k]+"\t"
        for j in range(70):
            output=output+str(mass_best[k][j])+"\t"+str(mass_inf[k][j])+"\t"+str(mass_med[k][j])+"\t"+str(mass_sup[k][j])+"\t"+str(chigal[k][j])+"\t"+str(chistar[k][j])+"\t"
        output=output+"\n"
    outfile.write(output)
    outfile.close()

    print("Total time for subfield: --- %s seconds ---" % (time.time() - start_timesubfield))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

