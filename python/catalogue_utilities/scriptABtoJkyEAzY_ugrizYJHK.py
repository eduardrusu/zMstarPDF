##########################
# The code is used to convert a photometric catalogue in the format expected by BPZ to the format expected by EaZy. The code handles non-detections (not yet non-exposures).
##########################

import numpy as np

#file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.cat"
file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpz.cat"

# zeropoint corrections suggested by BPZ:
u_corr = +1.40
g_corr = -0.03
r_corr = -0.00
i_corr =  0.00
z_corr = -0.01
Y_corr = -0.07
J_corr = -0.10
H_corr = -0.01
K_corr = +0.06

# conversion from Vega to AB; assuming that the input mags are in Vega
J_corr += 0.94 # computed by LePhare
H_corr += 1.35
K_corr += 1.83

id = 0
u = 1
u_err = 2
g = 3
g_err = 4
r = 5
r_err = 6
i = 7
i_err = 8
z = 9
z_err = 10
Y = 11
Y_err = 12
J = 13
J_err = 14
H = 15
H_err = 16
K = 17
K_err = 18

data = np.loadtxt(file,unpack=True)

# apply the corrections
data[u][np.abs(data[u]) != 99.0] += u_corr
data[g][np.abs(data[g]) != 99.0] += g_corr
data[r][np.abs(data[r]) != 99.0] += r_corr
data[i][np.abs(data[i]) != 99.0] += i_corr
data[z][np.abs(data[z]) != 99.0] += z_corr
data[Y][np.abs(data[Y]) != 99.0] += Y_corr
data[J][np.abs(data[J]) != 99.0] += J_corr
data[H][np.abs(data[H]) != 99.0] += H_corr
data[K][np.abs(data[K]) != 99.0] += K_corr

# make non-detection fluxes -99
data[u][np.abs(data[u]) == 99.0] = -99
data[g][np.abs(data[g]) == 99.0] = -99
data[r][np.abs(data[r]) == 99.0] = -99
data[i][np.abs(data[i]) == 99.0] = -99
data[z][np.abs(data[z]) == 99.0] = -99
data[Y][np.abs(data[Y]) == 99.0] = -99
data[J][np.abs(data[J]) == 99.0] = -99
data[H][np.abs(data[H]) == 99.0] = -99
data[K][np.abs(data[K]) == 99.0] = -99

# for the non-detections, replace error bars with the flux corresponsing to the corrected limiting mag
data[u_err][np.abs(data[u]) == 99.0] = 3631000000 * 10**(-(data[u_err][np.abs(data[u]) == 99.0] + u_corr)/2.5)
data[g_err][np.abs(data[g]) == 99.0] = 3631000000 * 10**(-(data[g_err][np.abs(data[g]) == 99.0] + g_corr)/2.5)
data[r_err][np.abs(data[r]) == 99.0] = 3631000000 * 10**(-(data[r_err][np.abs(data[r]) == 99.0] + r_corr)/2.5)
data[i_err][np.abs(data[i]) == 99.0] = 3631000000 * 10**(-(data[i_err][np.abs(data[i]) == 99.0] + i_corr)/2.5)
data[z_err][np.abs(data[z]) == 99.0] = 3631000000 * 10**(-(data[z_err][np.abs(data[z]) == 99.0] + z_corr)/2.5)
data[Y_err][np.abs(data[Y]) == 99.0] = 3631000000 * 10**(-(data[Y_err][np.abs(data[Y]) == 99.0] + Y_corr)/2.5)
data[J_err][np.abs(data[J]) == 99.0] = 3631000000 * 10**(-(data[J_err][np.abs(data[J]) == 99.0] + J_corr)/2.5)
data[H_err][np.abs(data[H]) == 99.0] = 3631000000 * 10**(-(data[H_err][np.abs(data[H]) == 99.0] + H_corr)/2.5)
data[K_err][np.abs(data[K]) == 99.0] = 3631000000 * 10**(-(data[K_err][np.abs(data[K]) == 99.0] + K_corr)/2.5)

# make minimum delta mag 0.01
data[u_err][(np.abs(data[u]) != 99.0) & (np.abs(data[u_err]) == 0.00)] = 0.01
data[g_err][(np.abs(data[g]) != 99.0) & (np.abs(data[g_err]) == 0.00)] = 0.01
data[r_err][(np.abs(data[r]) != 99.0) & (np.abs(data[r_err]) == 0.00)] = 0.01
data[i_err][(np.abs(data[i]) != 99.0) & (np.abs(data[i_err]) == 0.00)] = 0.01
data[z_err][(np.abs(data[z]) != 99.0) & (np.abs(data[z_err]) == 0.00)] = 0.01
data[Y_err][(np.abs(data[Y]) != 99.0) & (np.abs(data[Y_err]) == 0.00)] = 0.01
data[J_err][(np.abs(data[J]) != 99.0) & (np.abs(data[J_err]) == 0.00)] = 0.01
data[H_err][(np.abs(data[H]) != 99.0) & (np.abs(data[H_err]) == 0.00)] = 0.01
data[K_err][(np.abs(data[K]) != 99.0) & (np.abs(data[K_err]) == 0.00)] = 0.01

# a small number of objects in BPZ have good mags, but error on the mag 99; those objects should be good, and Ideally I would fix the errors one by one through closer examination. Here I just replace their errors with 1 mag
data[u_err][np.abs(data[u_err]) == 99.00] = 1.00
data[g_err][np.abs(data[g_err]) == 99.00] = 1.00
data[r_err][np.abs(data[r_err]) == 99.00] = 1.00
data[i_err][np.abs(data[i_err]) == 99.00] = 1.00
data[z_err][np.abs(data[z_err]) == 99.00] = 1.00
data[Y_err][np.abs(data[Y_err]) == 99.00] = 1.00
data[J_err][np.abs(data[J_err]) == 99.00] = 1.00
data[H_err][np.abs(data[H_err]) == 99.00] = 1.00
data[K_err][np.abs(data[K_err]) == 99.00] = 1.00

# convert AB -> Jky
x=data[u][np.abs(data[u]) != 99.0]
y=data[u_err][np.abs(data[u]) != 99.0]
data[u_err][np.abs(data[u]) != 99.0] = 3631000000 * (10**(-(data[u][np.abs(data[u]) != 99.0] - data[u_err][np.abs(data[u]) != 99.0])/2.5) - 10**(-(data[u][np.abs(data[u]) != 99.0] + data[u_err][np.abs(data[u]) != 99.0])/2.5)) / 2
data[u][np.abs(data[u]) != 99.0] = 3631000000 * 10**(-data[u][np.abs(data[u]) != 99.0]/2.5)
data[g_err][np.abs(data[g]) != 99.0] = 3631000000 * (10**(-(data[g][np.abs(data[g]) != 99.0] - data[g_err][np.abs(data[g]) != 99.0])/2.5) - 10**(-(data[g][np.abs(data[g]) != 99.0] + data[g_err][np.abs(data[g]) != 99.0])/2.5)) / 2
data[g][np.abs(data[g]) != 99.0] = 3631000000 * 10**(-data[g][np.abs(data[g]) != 99.0]/2.5)
data[r_err][np.abs(data[r]) != 99.0] = 3631000000 * (10**(-(data[r][np.abs(data[r]) != 99.0] - data[r_err][np.abs(data[r]) != 99.0])/2.5) - 10**(-(data[r][np.abs(data[r]) != 99.0] + data[r_err][np.abs(data[r]) != 99.0])/2.5)) / 2
data[r][np.abs(data[r]) != 99.0] = 3631000000 * 10**(-data[r][np.abs(data[r]) != 99.0]/2.5)
data[i_err][np.abs(data[i]) != 99.0] = 3631000000 * (10**(-(data[i][np.abs(data[i]) != 99.0] - data[i_err][np.abs(data[i]) != 99.0])/2.5) - 10**(-(data[i][np.abs(data[i]) != 99.0] + data[i_err][np.abs(data[i]) != 99.0])/2.5)) / 2
data[i][np.abs(data[i]) != 99.0] = 3631000000 * 10**(-data[i][np.abs(data[i]) != 99.0]/2.5)
data[z_err][np.abs(data[z]) != 99.0] = 3631000000 * (10**(-(data[z][np.abs(data[z]) != 99.0] - data[z_err][np.abs(data[z]) != 99.0])/2.5) - 10**(-(data[z][np.abs(data[z]) != 99.0] + data[z_err][np.abs(data[z]) != 99.0])/2.5)) / 2
data[z][np.abs(data[z]) != 99.0] = 3631000000 * 10**(-data[z][np.abs(data[z]) != 99.0]/2.5)
data[Y_err][np.abs(data[Y]) != 99.0] = 3631000000 * (10**(-(data[Y][np.abs(data[Y]) != 99.0] - data[Y_err][np.abs(data[Y]) != 99.0])/2.5) - 10**(-(data[Y][np.abs(data[Y]) != 99.0] + data[Y_err][np.abs(data[Y]) != 99.0])/2.5)) / 2
data[Y][np.abs(data[Y]) != 99.0] = 3631000000 * 10**(-data[Y][np.abs(data[Y]) != 99.0]/2.5)
data[J_err][np.abs(data[J]) != 99.0] = 3631000000 * (10**(-(data[J][np.abs(data[J]) != 99.0] - data[J_err][np.abs(data[J]) != 99.0])/2.5) - 10**(-(data[J][np.abs(data[J]) != 99.0] + data[J_err][np.abs(data[J]) != 99.0])/2.5)) / 2
data[J][np.abs(data[J]) != 99.0] = 3631000000 * 10**(-data[J][np.abs(data[J]) != 99.0]/2.5)
data[H_err][np.abs(data[H]) != 99.0] = 3631000000 * (10**(-(data[H][np.abs(data[H]) != 99.0] - data[H_err][np.abs(data[H]) != 99.0])/2.5) - 10**(-(data[H][np.abs(data[H]) != 99.0] + data[H_err][np.abs(data[H]) != 99.0])/2.5)) / 2
data[H][np.abs(data[H]) != 99.0] = 3631000000 * 10**(-data[H][np.abs(data[H]) != 99.0]/2.5)
data[K_err][np.abs(data[K]) != 99.0] = 3631000000 * (10**(-(data[K][np.abs(data[K]) != 99.0] - data[K_err][np.abs(data[K]) != 99.0])/2.5) - 10**(-(data[K][np.abs(data[K]) != 99.0] + data[K_err][np.abs(data[K]) != 99.0])/2.5)) / 2
data[K][np.abs(data[K]) != 99.0] = 3631000000 * 10**(-data[K][np.abs(data[K]) != 99.0]/2.5)

                          
# 146  u_DEC: 200
# 147  g_DEC: 300
# 148  r_DEC: 300
# 149  i_DEC: 600
# 150  z_DEC: 400
# 151  Y_DEC: 200
# 152  J_HAWKI: 150 AB-Vega= 1.06
# 153  H_HAWKI: 220 AB-Vega= 1.34
# 154  Ks_HAWKI: 300 AB-Vega= 1.78
                                                                                                                                                                                                              
fileout = file[:-11] + "_forEazy.cat"
str = "id  F146      E146      F147      E147       F148      E148       F149      E149       F150     E150      F151     E151      F152     E152      F153     E153      F154     E154"
dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f')


