# Varies the lens redshift and get resulting parameters. Uses mag2mag to find absolute mag

import sys
import os
import numpy as np
import corner
import pylab as plt

sample=np.linspace(0.1,2.2,22)

listchi = np.array([])
list = np.array([])
listmag = np.array([])
listvel = np.array([])

filein = "inlensz.input"
fileimag = "outlensz_point.dat"
fileout = "outlensz_optresult.dat"
fileres = "samplezveldisp.dat"

for i in range(len(sample)):
    with open(filein, 'r') as f:
        glafic = f.readlines()
    glafic[9-1] = glafic[9-1].replace(glafic[9-1].split()[1],str(sample[i]))
    with open(filein, 'w') as f:
        f.writelines(glafic)
    os.system("glafic %s" % filein)

    with open(fileimag, 'r') as f:
        interest = f.readlines()
    if len(interest) == 5: # accept only if the model creates 4 images
        with open(fileout, 'r') as f:
            interest = f.readlines()
        for j in range(len(interest)):
            if "lens   sie" in interest[j]:
                value = float(interest[j].split()[2])
            if "chi" in interest[j]:
                valuechi = float(interest[j].split()[2])
        listchi = np.append(listchi,valuechi)
        list = np.append(list,value)
#np.savetxt(fileres,np.c_[sample,list,listchi],fmt="%.1f %d %.2f")

import subprocess
for i in range(len(sample)):
    output = subprocess.check_output('$HOME/GITHUB/LensPop/stellarpop/mag2mag.py -vega -T El_cww -m1 17.33 -f1 VISTA_Ks -z1 %s -f2 R_Cousins -z2 0.0; exit 0' % str(sample[i]),stderr=subprocess.STDOUT,shell=True)
    listmag = np.append(listmag,float(output))
    listvel = np.append(listvel,150*((10**((4.65-float(output))/2.5))/(10**10))**0.25)

#Galactic dynamics: log10(sigma/150km/s) = 0.25log10(L_R/(10^10*h_7^-2Lsun)) ->
#sigma = 150km/s*[L_R/(10^10*Lsun)]^0.25
#-2.5log10L_R[Lsun] = M_R - M_Rsun ->
#L_R[Lsun] =  10^[(M_Rsun - M_R)/2.5] ->
#sigma = 150km/s*[10^[(M_Rsun - M_R)/2.5]/10^10]^0.25
#M_R Sun in AB from http://mips.as.arizona.edu/~cnaw/sun.html is 4.65
# scatter in sigma = 0.25*sigma

np.savetxt(fileres,np.c_[sample,list,listchi,listmag,listvel],fmt="%.1f %d %.2f %.2f %d")
plt.plot(sample,list)
plt.plot(sample,listvel)
plt.show()



