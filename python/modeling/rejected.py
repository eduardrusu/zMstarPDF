# The code fixes the MCMC chains from Hostlens/Glafic by accounting for the trace of the rejected points. It needs two inputs:
# - the mcmc file exported by hostlens/glafic
# - a copy of the execution screen which displays the mcmc reject/accepted points; edit it to remove the unnecessary lines (header and footer)

import sys
import os
import numpy as np
import corner

prefix = "ilens_out_mcmc"
filescreen = prefix + "screen.dat"
filemcmc = prefix + ".dat"
fileout = prefix + "withrejected.dat"

with open(filescreen, 'r') as f:
    screen = f.readlines()
    f.close()

with open(filemcmc, 'r') as f:
    mcmc = f.readlines()
    f.close()

out = [] # will contain the complete chains
i = 0 # index in mcmc file
init = 0 # whether or not it found the first successful point
for j in range(np.shape(screen)[0]): # index in screen file
    if init == 1 and "[rejected]" in screen[j]:
        out.append(mcmc[i])
    if init == 1 and "[rejected]" not in screen[j]:
        i += 1
        out.append(mcmc[i])
    if init == 0 and "[rejected]" not in screen[j]:
        init = 1
        out.append(mcmc[i])
        i += 1

with open(fileout, 'w') as f:
    f.writelines(out)
    f.close()
