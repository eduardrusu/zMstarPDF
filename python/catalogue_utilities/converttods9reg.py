# read RA and DEC from a file and output a corresponding DS9 regions file
import numpy as np
import sys
input = str(sys.argv[1])
out = input[:-4] + '_reg.cat'

cat = np.loadtxt(input, usecols = [0,1], unpack=True)
np.savetxt(out,cat.T,fmt='fk5;circle(   %.6f , %.6f,    0.0005)')
