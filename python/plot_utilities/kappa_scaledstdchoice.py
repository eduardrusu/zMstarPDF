# uses the output of kappa_medsigsim.py to decide the based weight combination

import numpy as np
from functools import reduce

file1 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_0/medstdbias_base120.dat"
file2 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_1/medstdbias_base120.dat"
file3 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_2/medstdbias_base120.dat"
file4 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_0/medstdbias_base45.dat"
file5 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_1/medstdbias_base45.dat"
file6 = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_2/medstdbias_base45.dat"

f1 = np.genfromtxt(file1,usecols=[0,3,4,5],dtype='S200,f,f,f')
f2 = np.genfromtxt(file2,usecols=[0,3,4,5],dtype='S200,f,f,f')
f3 = np.genfromtxt(file3,usecols=[0,3,4,5],dtype='S200,f,f,f')
f4 = np.genfromtxt(file1,usecols=[0,3,4,5],dtype='S200,f,f,f')
f5 = np.genfromtxt(file2,usecols=[0,3,4,5],dtype='S200,f,f,f')
f6 = np.genfromtxt(file3,usecols=[0,3,4,5],dtype='S200,f,f,f')
f1name = np.array([])
f1los = np.array([])
f1lines = np.array([])
f1std = np.array([])
f2name = np.array([])
f2los = np.array([])
f2lines = np.array([])
f2std = np.array([])
f3name = np.array([])
f3los = np.array([])
f3lines = np.array([])
f3std = np.array([])
f4name = np.array([])
f4los = np.array([])
f4lines = np.array([])
f4std = np.array([])
f5name = np.array([])
f5los = np.array([])
f5lines = np.array([])
f5std = np.array([])
f6name = np.array([])
f6los = np.array([])
f6lines = np.array([])
f6std = np.array([])
for i in range(len(f1)):
    f1name = np.append(f1name,f1[i][0])
    f1los = np.append(f1los,f1[i][1])
    f1lines = np.append(f1lines,f1[i][2])
    f1std = np.append(f1std,f1[i][3])
for i in range(len(f2)):
    f2name = np.append(f2name,f2[i][0])
    f2los = np.append(f2los,f2[i][1])
    f2lines = np.append(f2lines,f2[i][2])
    f2std = np.append(f2std,f2[i][3])
for i in range(len(f3)):
    f3name = np.append(f3name,f3[i][0])
    f3los = np.append(f3los,f3[i][1])
    f3lines = np.append(f3lines,f3[i][2])
    f3std = np.append(f3std,f3[i][3])
for i in range(len(f4)):
    f4name = np.append(f4name,f4[i][0])
    f4los = np.append(f4los,f4[i][1])
    f4lines = np.append(f4lines,f4[i][2])
    f4std = np.append(f4std,f4[i][3])
for i in range(len(f5)):
    f5name = np.append(f5name,f5[i][0])
    f5los = np.append(f5los,f5[i][1])
    f5lines = np.append(f5lines,f5[i][2])
    f5std = np.append(f5std,f5[i][3])
for i in range(len(f6)):
    f6name = np.append(f6name,f6[i][0])
    f6los = np.append(f6los,f6[i][1])
    f6lines = np.append(f6lines,f6[i][2])
    f6std = np.append(f6std,f6[i][3])

matches = reduce(np.intersect1d, (f1name,f2name,f3name,f4name,f5name,f6name))
matches = np.sort(matches) # alphabetical sort

out = np.zeros(21)
for i in range(len(matches)):
    med = np.around(np.mean([np.around(float(f1std[np.where(f1name == matches[i])][0]),3),np.around(float(f2std[np.where(f2name == matches[i])][0]),3),np.around(float(f3std[np.where(f3name == matches[i])][0]),3),np.around(float(f4std[np.where(f4name == matches[i])][0]),3),np.around(float(f5std[np.where(f5name == matches[i])][0]),3),np.around(float(f6std[np.where(f6name == matches[i])][0]),3)]),3)
    std = np.around(np.std([np.around(float(f1std[np.where(f1name == matches[i])][0]),3),np.around(float(f2std[np.where(f2name == matches[i])][0]),3),np.around(float(f3std[np.where(f3name == matches[i])][0]),3),np.around(float(f4std[np.where(f4name == matches[i])][0]),3),np.around(float(f5std[np.where(f5name == matches[i])][0]),3),np.around(float(f6std[np.where(f6name == matches[i])][0]),3)]),3)
    out_ = np.array([matches[i],int(f1los[np.where(f1name == matches[i])][0]),int(f2los[np.where(f2name == matches[i])][0]),int(f3los[np.where(f3name == matches[i])][0]),int(f4los[np.where(f4name == matches[i])][0]),int(f5los[np.where(f5name == matches[i])][0]),int(f6los[np.where(f6name == matches[i])][0]),int(f1lines[np.where(f1name == matches[i])][0]),int(f2lines[np.where(f2name == matches[i])][0]),int(f3lines[np.where(f3name == matches[i])][0]),int(f4lines[np.where(f4name == matches[i])][0]),int(f5lines[np.where(f5name == matches[i])][0]),int(f6lines[np.where(f6name == matches[i])][0]),np.around(float(f1std[np.where(f1name == matches[i])][0]),3),np.around(float(f2std[np.where(f2name == matches[i])][0]),3),np.around(float(f3std[np.where(f3name == matches[i])][0]),3),np.around(float(f4std[np.where(f4name == matches[i])][0]),3),np.around(float(f5std[np.where(f5name == matches[i])][0]),3),np.around(float(f6std[np.where(f6name == matches[i])][0]),3),med,std])
    out = np.c_[out,out_]

fout = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/scaledstdchoice.dat"
np.savetxt(fout,out.T,fmt='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s')

# the RMS from inferkappasimbiasscript1.py and inferkappasimbiasscript2.py:
#45: np.std([1.071,1.020,0.995,1.141,1.043,0.901,0.999,0.966,1.108,1.013,0.844,0.972,1.052,1.033,0.935]) = 0.074
#120: np.std([1.168,0.992,1.069,1.008,0.955,0.955,1.071,0.932,1.031,0.977]) = 0.068
rms2 = 2 * 0.045
