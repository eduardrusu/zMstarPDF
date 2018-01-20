# Samples from an observed lens parameter scatter (such as the lens position) and optimizes, then gathers the output of interest. In this particular case I sample from the lens position, and I adopt the limiting chi^2=9.21[ + best original fit chi^2] (Avni 1976 for 2 parametes at 3 sigma). When reading the chi^2 I only consider the point contribution, not the lens. For the lens contribution I compute it using the observed distance prior

import sys
import os
import numpy as np
import corner

sample = 10000
posx = 0.74
posx_sigma = 0.04
posy = 0.75
posy_sigma = 0.10
x = np.random.normal(posx, posx_sigma, sample)
y = np.random.normal(posy, posy_sigma, sample)
limchi = 0.13 + 9.21

listchi = np.array([])
listx = np.array([])
listy = np.array([])
list1 = np.array([])
list2 = np.array([])
list3 = np.array([])
list4 = np.array([])

filein = "insample.input"
fileprior = "priorsample.dat"
fileout = "outsample_optresult.dat"
fileimag = "outsample_point.dat"

for i in range(sample):
    with open(filein, 'r') as f:
        glafic = f.readlines()
    glafic[28-1] = glafic[28-1].replace(glafic[28-1].split()[3],str(x[i]))
    glafic[28-1] = glafic[28-1].replace(glafic[28-1].split()[4],str(y[i]))
    glafic[29-1] = glafic[29-1].replace(glafic[29-1].split()[3],str(x[i]))
    glafic[29-1] = glafic[29-1].replace(glafic[29-1].split()[4],str(y[i]))
    with open(filein, 'w') as f:
        f.writelines(glafic)
    with open(fileprior, 'r') as f:
        prior = f.readlines()
    prior[1-1] = "gauss lens 1 2 %s 0.01 \n" % str(x[i])
    prior[2-1] = "gauss lens 1 3 %s 0.01 \n" % str(y[i])
    with open(fileprior, 'w') as f:
        f.writelines(prior)
    os.system("glafic %s" % filein)
    with open(fileimag, 'r') as f:
        interest = f.readlines()
    if len(interest) == 5: # accept only if the model creates 4 images
        with open(fileout, 'r') as f:
            interest = f.readlines()
        for j in range(len(interest)):
            if "point no 1" in interest[j]:
                value = float(interest[j].split()[4]) + float(interest[j].split()[6]) + ((x[i]-posx)/posx_sigma)**2 + ((y[i]-posy)/posy_sigma)**2
        if value <= limchi:
            for j in range(len(interest)):
                if "sie" in interest[j]:
                    value1 = float(interest[j].split()[5])
                    value2 = float(interest[j].split()[6])
            for j in range(len(interest)):
                if "pert" in interest[j]:
                    value3 = float(interest[j].split()[5])
                    value4 = float(interest[j].split()[6])
            listchi = np.append(listchi,value)
            listx = np.append(listx,x[i])
            listy = np.append(listy,y[i])
            list1 = np.append(list1,value1)
            list2 = np.append(list2,value2)
            list3 = np.append(list3,value3)
            list4 = np.append(list4,value4)
np.savetxt("sample.dat",np.c_[listchi,listx,listy,list1,list2,list3,list4],fmt="%.6e %.6e %.6e %.6e %.6e %.6e %.6e")

sample = np.loadtxt("sample.dat",unpack=True)
#sample1 = sample[0][sample[0] < 2.3+0.13]
#sample2 = sample[1][sample[0] < 2.3+0.13]
#sample3 = sample[2][sample[0] < 2.3+0.13]
#sample4 = sample[3][sample[0] < 2.3+0.13]
#sample5 = sample[4][sample[0] < 2.3+0.13]
#sample6 = sample[5][sample[0] < 2.3+0.13]
#sample7 = sample[6][sample[0] < 2.3+0.13]
#sample=np.vstack((sample1,sample2,sample3,sample4,sample5,sample6,sample7))
figure = corner.corner(sample[1:np.shape(sample)[0]].T, labels=np.linspace(1,np.shape(sample)[0],np.shape(sample)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig("sample.png", dpi=100)



