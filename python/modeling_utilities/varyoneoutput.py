# Given a out_optresult.dat file output fom glafic varyone, run glafic in order to plot the variation of the corresponding images.

import sys
import os
import numpy as np
import corner

filein = str(sys.argv[1]) # glafic input file
fileout = str(sys.argv[2]) # out_optresult.dat file
with open(fileout, 'r') as f:
    out = f.readlines()

listlens = np.array([])
list = np.array([])
i = 1
index = 12 # first line containing the lens parameters in the file

################# this is in case I'm looking for lens properties:
while index + 2 < len(out):
    listlens = np.append(listlens,float(out[index-1].split()[2]))
    index += 16 # the number of lines until the next lens parameters
print listlens

################# this is only in case I'm looking for image properties:
while index < 0:
#while index + 2 < len(out):
    with open(filein, 'r') as f:
        glafic = f.readlines()
    glafic[10 - 1] = "prefix        out%s" % str(i) + "\n"
    glafic[29 - 1] = out[index-1]
    glafic[30 - 1] = out[index]
    glafic[31 - 1] = out[index+1]
    glafic[32 - 1] = out[index+2]
    glafic[51 - 1] = "findimg" + "\n"
    with open(filein, 'w') as f:
        f.writelines(glafic)
        #f.close()
    os.system("glafic %s" % filein)
    data = np.loadtxt("out%s_point.dat" % str(i))
    list = np.append(list,data[4][3])
    os.system("rm out%s_point.dat" % str(i))
    index += 16 # the number of lines until the next lens parameters
    i += 1

print list


