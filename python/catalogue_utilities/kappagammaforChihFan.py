# CE Rusu, July 8 2018
# run as e.g.: python /lfs08/rusucs/code/kappamed_insertstarsnobeta.py GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f

import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
import time

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...") # use for debugging

def readbinary(replacestr):
    replace = plane + replacestr
    os.system("sed \"11s/.*/  const char kappa_file_name[]   = \\\"\%s\\\";/\" readKappaBinary.c > readKappaBinary_%s.c_" % (replace,plane))
    os.system("sed \"35s/.*/  fpt = fopen (\\\"kappa_values_%s.dat\\\", \\\"w\\\");/\"  readKappaBinary_%s.c_ > readKappaBinary_%s.c" % (plane,plane,plane))
    os.system("rm -f readKappaBinary_%s.c_" % plane)
    os.system("gcc readKappaBinary_%s.c -o compiled_%s.out" % (plane,plane))
    os.system("./compiled_%s.out" % plane)
    os.system("rm -f readKappaBinary_%s.c" % plane)
    os.system("rm -f compiled_%s.out" % plane)

start_time = time.time()

plane = str(sys.argv[1])

#if lens == "B1608":
#    z_s = 1.39
#    pln = 37 # MS simulation plane
#if lens == "HE0435":
#    z_s = 1.69
#    z_l = 0.455
#if lens == "WFI2033":
#    z_s = 1.66
#    z_l = 0.66
#if lens == "HE1104":
#    z_s = 2.32
#pln = 30 & 31
#if lens == "RX1131":
#    z_s = 0.66
##pln = 45 & 46
#if lens == "J1206":
#    z_s = 1.79
#    pln = 34
#if lens == "PG1115":
#    z_s = 1.72
pln = 34

rootkappaplanes = "/lfs08/rusucs/kappaplanes/"

start_readkappa = time.time()

if str(pln) in plane:
    os.chdir(rootkappaplanes)
    readbinary(".kappa")
    kappa = np.loadtxt("kappa_values_%s.dat" % plane, usecols = [1], unpack=True)
    readbinary(".gamma_1")
    gamma1 = np.loadtxt("kappa_values_%s.dat" % plane, usecols = [1], unpack=True)
    readbinary(".gamma_2")
    gamma2 = np.loadtxt("kappa_values_%s.dat" % plane, usecols = [1], unpack=True)
    kappagamma = np.c_[kappa,np.sqrt(gamma1 ** 2 + gamma2 ** 2)]
    os.system("rm -f kappa_values_%s.dat" % plane)
    np.savetxt("kappagamma_values_%s.dat" % plane,kappagamma,fmt="%e %e")
else: sys.exit('Wrong MS plane for this lens!!!')

print(" Field done in --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
