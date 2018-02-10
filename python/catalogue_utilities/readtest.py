import numpy as np
import sys
import os
from os import system

def readbinary(replacestr,plane):
    replace = plane + replacestr
    os.system("sed \"11s/.*/  const char kappa_file_name[]   = \\\"\%s\\\";/\" readKappaBinary.c > readKappaBinary_%s.c_" % (replace,plane))
    os.system("sed \"35s/.*/  fpt = fopen (\\\"kappa_values_%s.dat\\\", \\\"w\\\");/\"  readKappaBinary_%s.c_ > readKappaBinary_%s.c" % (plane,plane,plane))
    os.system("rm -f readKappaBinary_%s.c_" % plane)
    os.system("gcc readKappaBinary_%s.c -o compiled_%s.out" % (plane,plane))
    os.system("./compiled_%s.out" % plane)
    os.system("rm -f readKappaBinary_%s.c" % plane)
    os.system("rm -f compiled_%s.out" % plane)

rootkappaplanes = "/lfs08/rusucs/kappaplanes/"
plane = 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_35_f'
pln = 35

if str(pln) in plane:
    os.chdir(rootkappaplanes)
    if 'kappa_values_%s.dat' % plane in os.listdir('.'): readbinary(".kappa",plane)
    pos1D,kappa = np.loadtxt("kappa_values_%s.dat" % plane, unpack=True)
    readbinary(".gamma_1",plane)
    gamma1 = np.loadtxt("kappa_values_%s.dat" % plane, usecols = [1], unpack=True)
    readbinary(".gamma_2",plane)
    gamma2 = np.loadtxt("kappa_values_%s.dat" % plane, usecols = [1], unpack=True)
    kappagamma = np.c_[pos1D,kappa,gamma1,gamma2]
    #os.system("rm -f kappa_values_%s.dat" % plane)
else: sys.exit('Wrong MS plane for this lens!!!')

