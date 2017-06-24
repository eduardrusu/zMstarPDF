##########################
# The code is used to combine two LePhare output catalogues: one which considers non-detections, and one which assumes they are non-observations. This is because some non-detections would not execute.
#  Use as: combinelephare.py /Users/eduardrusu/lephare_dev/test/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forlepharewithbpz.cat.MAG_BC03_I09.lephareout
##########################

import numpy as np
import sys

file = str(sys.argv[1])

file_noobs = file[:-28] + "_noobs.cat.MAG_BC03_I09.lephareout"

id = 0
z = 1
chi_best = 5
chi_star = 6
mass_best = 43
mass_inf = 44
mass_med = 45
mass_sup = 46

data = np.loadtxt(file,usecols=(id,z,chi_best,chi_star,mass_best,mass_inf,mass_med,mass_sup),unpack=True)
data_noobs = np.loadtxt(file_noobs,usecols=(id,z,chi_best,chi_star,mass_best,mass_inf,mass_med,mass_sup),unpack=True)

for i in range(len(data[id])):
    if data[:,i][z] == -99.0:
        data[:,i] = data_noobs[:,i]

fileout = file[:-28] + "_combined.cat.MAG_BC03_I09.lephareout"
str = "ID \t z \t chi_best \t chi_star \t mass_best \t mass_inf \t mass_med \t mass_sup"
dataout = np.c_[data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.3f \t %.3f \t %.3f \t %.3f ')

