##########################
# The code is used to convert a photometric catalogue in the format expected by BPZ to the format expected by EaZy. The code handles non-detections (not yet non-exposures).
##########################

import numpy as np
import sys

file = str(sys.argv[1])
#file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.cat"
#file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.cat"

# zeropoint corrections suggested by BPZ:
if 'detectin_iunconv' in file:
    g_corr = -0.077
    r_corr = -0.004
    i_corr =  0.000
    K_corr = -0.022
if 'detectin_iconv' in file:
    g_corr = -0.061
    r_corr = -0.007
    i_corr =  0.000
    K_corr = -0.046

# conversion from Vega to AB; assuming that the input mags are in Vega
K_corr += 1.83

if 'detectin_iunconv' in file:
    id = 0
    g = 22
    g_err = 23
    r = 24
    r_err = 25
    i = 26
    i_err = 27
    K = 28
    K_err = 29
if 'detectin_iconv' in file:
    id = 0
    g = 20
    g_err = 21
    r = 22
    r_err = 23
    i = 24
    i_err = 25
    K = 26
    K_err = 27

data = np.loadtxt(file,unpack=True)

# apply the corrections
data[g][np.abs(data[g]) != 99.0] += g_corr
data[r][np.abs(data[r]) != 99.0] += r_corr
data[i][np.abs(data[i]) != 99.0] += i_corr
data[K][np.abs(data[K]) != 99.0] += K_corr

# make non-detection fluxes -99
data[g][np.abs(data[g]) == 99.0] = -99
data[r][np.abs(data[r]) == 99.0] = -99
data[i][np.abs(data[i]) == 99.0] = -99
data[K][np.abs(data[K]) == 99.0] = -99

# for the non-detections, replace error bars with the flux corresponsing to the corrected limiting mag
data[g_err][np.abs(data[g]) == 99.0] = 3631000000 * 10**(-(data[g_err][np.abs(data[g]) == 99.0] + g_corr)/2.5)
data[r_err][np.abs(data[r]) == 99.0] = 3631000000 * 10**(-(data[r_err][np.abs(data[r]) == 99.0] + r_corr)/2.5)
data[i_err][np.abs(data[i]) == 99.0] = 3631000000 * 10**(-(data[i_err][np.abs(data[i]) == 99.0] + i_corr)/2.5)
data[K_err][np.abs(data[K]) == 99.0] = 3631000000 * 10**(-(data[K_err][np.abs(data[K]) == 99.0] + K_corr)/2.5)

# make minimum delta mag 0.01
data[g_err][(np.abs(data[g]) != 99.0) & (np.abs(data[g_err]) == 0.00)] = 0.01
data[r_err][(np.abs(data[r]) != 99.0) & (np.abs(data[r_err]) == 0.00)] = 0.01
data[i_err][(np.abs(data[i]) != 99.0) & (np.abs(data[i_err]) == 0.00)] = 0.01
data[K_err][(np.abs(data[K]) != 99.0) & (np.abs(data[K_err]) == 0.00)] = 0.01

# a small number of objects in BPZ have good mags, but error on the mag 99; those objects should be good, and Ideally I would fix the errors one by one through closer examination. Here I just replace their errors with 1 mag
data[g_err][np.abs(data[g_err]) == 99.00] = 1.00
data[r_err][np.abs(data[r_err]) == 99.00] = 1.00
data[i_err][np.abs(data[i_err]) == 99.00] = 1.00
data[K_err][np.abs(data[K_err]) == 99.00] = 1.00

# convert AB -> Jky
data[g_err][np.abs(data[g]) != 99.0] = 3631000000 * (10**(-(data[g][np.abs(data[g]) != 99.0] - data[g_err][np.abs(data[g]) != 99.0])/2.5) - 10**(-(data[g][np.abs(data[g]) != 99.0] + data[g_err][np.abs(data[g]) != 99.0])/2.5)) / 2
data[g][np.abs(data[g]) != 99.0] = 3631000000 * 10**(-data[g][np.abs(data[g]) != 99.0]/2.5)
data[r_err][np.abs(data[r]) != 99.0] = 3631000000 * (10**(-(data[r][np.abs(data[r]) != 99.0] - data[r_err][np.abs(data[r]) != 99.0])/2.5) - 10**(-(data[r][np.abs(data[r]) != 99.0] + data[r_err][np.abs(data[r]) != 99.0])/2.5)) / 2
data[r][np.abs(data[r]) != 99.0] = 3631000000 * 10**(-data[r][np.abs(data[r]) != 99.0]/2.5)
data[i_err][np.abs(data[i]) != 99.0] = 3631000000 * (10**(-(data[i][np.abs(data[i]) != 99.0] - data[i_err][np.abs(data[i]) != 99.0])/2.5) - 10**(-(data[i][np.abs(data[i]) != 99.0] + data[i_err][np.abs(data[i]) != 99.0])/2.5)) / 2
data[i][np.abs(data[i]) != 99.0] = 3631000000 * 10**(-data[i][np.abs(data[i]) != 99.0]/2.5)
data[K_err][np.abs(data[K]) != 99.0] = 3631000000 * (10**(-(data[K][np.abs(data[K]) != 99.0] - data[K_err][np.abs(data[K]) != 99.0])/2.5) - 10**(-(data[K][np.abs(data[K]) != 99.0] + data[K_err][np.abs(data[K]) != 99.0])/2.5)) / 2
data[K][np.abs(data[K]) != 99.0] = 3631000000 * 10**(-data[K][np.abs(data[K]) != 99.0]/2.5)
                                                                                                                                                                                                              
fileout = file[:-12] + "_forEazy.cat"
str = "id  F323      E323      F324      E324       F325      E325       F326      E326"
dataout = np.c_[data[id],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[K],data[K_err]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f')
