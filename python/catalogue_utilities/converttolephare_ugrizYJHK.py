##########################
# The code is used to convert raw matched photometric catalogues (with matched PSFs) from BPZ-expected input into LePhare-expected input. The first part of the code accounts for non-detections but not for non-observations. The second part considers all non-detections as non-observations. This is because some non-detections would not execute.The code requires "_withbpzeazy" files, which include available spectroscopic information.
##########################

import numpy as np 

#file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazy.cat"
#file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmags.cat"
#file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazy.cat"
file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmags.cat"

# zeropoint corrections suggested by lephare:
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

id = 8
u = 9
u_err = 10
g = 11
g_err = 12
r = 13
r_err = 14
i = 15
i_err = 16
z = 17
z_err = 18
Y = 19
Y_err = 20
J = 21
J_err = 22
H = 23
H_err = 24
K = 25
K_err = 26
photoz_bpz = 28
photoz_eazy = 48
specz = 40

#irac1 = 0
#irac1_err = 1
#irac2 = 2
#irac2_err = 3
#irac3 = 4
#irac3_err = 5
#irac4 = 6
#irac4_err = 7
irac1 = 74
irac1_err = 75
irac2 = 76
irac2_err = 77
irac3 = 78
irac3_err = 79
irac4 = 80
irac4_err = 81

data = np.loadtxt(file,unpack=True)

# If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
#data[irac1] = -99
#data[irac1_err] = -99
#data[irac2] = -99
#data[irac2_err] = -99
#data[irac3] = -99
#data[irac3_err] = -99
#data[irac4] = -99
#data[irac4_err] = -99

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

# use photoz where specz is unavailable, or for stars
data[photoz_bpz][data[specz] > 0] = data[specz][data[specz] > 0]
data[photoz_eazy][data[specz] > 0] = data[specz][data[specz] > 0]

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

# correct the limiting mags
data[u_err][np.abs(data[u]) == 99.0] += u_corr
data[g_err][np.abs(data[g]) == 99.0] += g_corr
data[r_err][np.abs(data[r]) == 99.0] += r_corr
data[i_err][np.abs(data[i]) == 99.0] += i_corr
data[z_err][np.abs(data[z]) == 99.0] += z_corr
data[Y_err][np.abs(data[Y]) == 99.0] += Y_corr
data[J_err][np.abs(data[J]) == 99.0] += J_corr
data[H_err][np.abs(data[H]) == 99.0] += H_corr
data[K_err][np.abs(data[K]) == 99.0] += K_corr

# the format for nondetections is error=-1.0 and magnitude at 1-sigma
data[u][np.abs(data[u]) == 99.0] = data[u_err][np.abs(data[u]) == 99.0]
data[g][np.abs(data[g]) == 99.0] = data[g_err][np.abs(data[g]) == 99.0]
data[r][np.abs(data[r]) == 99.0] = data[r_err][np.abs(data[r]) == 99.0]
data[i][np.abs(data[i]) == 99.0] = data[i_err][np.abs(data[i]) == 99.0]
data[z][np.abs(data[z]) == 99.0] = data[z_err][np.abs(data[z]) == 99.0]
data[Y][np.abs(data[Y]) == 99.0] = data[Y_err][np.abs(data[Y]) == 99.0]
data[J][np.abs(data[J]) == 99.0] = data[J_err][np.abs(data[J]) == 99.0]
data[H][np.abs(data[H]) == 99.0] = data[H_err][np.abs(data[H]) == 99.0]
data[K][np.abs(data[K]) == 99.0] = data[K_err][np.abs(data[K]) == 99.0]
data[u_err][np.abs(data[u_err]) > 20] = -1.0
data[g_err][np.abs(data[g_err]) > 20] = -1.0
data[r_err][np.abs(data[r_err]) > 20] = -1.0
data[i_err][np.abs(data[i_err]) > 20] = -1.0
data[z_err][np.abs(data[z_err]) > 20] = -1.0
data[Y_err][np.abs(data[Y_err]) > 20] = -1.0
data[J_err][np.abs(data[J_err]) > 20] = -1.0
data[H_err][np.abs(data[H_err]) > 20] = -1.0
data[K_err][np.abs(data[K_err]) > 20] = -1.0

# LePhare thinks error bars = 0 means non-detection, so fix this
data[u_err][data[u_err] == 0.00] = 0.01
data[g_err][data[g_err] == 0.00] = 0.01
data[r_err][data[r_err] == 0.00] = 0.01
data[i_err][data[i_err] == 0.00] = 0.01
data[z_err][data[z_err] == 0.00] = 0.01
data[Y_err][data[Y_err] == 0.00] = 0.01
data[J_err][data[J_err] == 0.00] = 0.01
data[H_err][data[H_err] == 0.00] = 0.01
data[K_err][data[K_err] == 0.00] = 0.01

#fileout = file[:-16] + "_forlepharewithbpz.cat"
fileout = file[:-4] + "_forlepharewithbpzIRAC.cat"
str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t K_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[photoz_bpz]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')

#fileout = file[:-16] + "_forlepharewitheazy.cat"
fileout = file[:-4] + "_forlepharewitheazyIRAC.cat"
str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t K_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[photoz_eazy]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')

##########################
# The code is used to convert raw matched photometric catalogues (with matched PSFs) from BPZ-expected input into LePhare-expected input. It considers all non-detections as non-observations. This is because some non-detections would not execute.
##########################

data = np.loadtxt(file,unpack=True)

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

# If not observed in a specific band, negative values (-99,-99) can be used for (mag,error)
#data[irac1] = -99
#data[irac1_err] = -99
#data[irac2] = -99
#data[irac2_err] = -99
#data[irac3] = -99
#data[irac3_err] = -99
#data[irac4] = -99
#data[irac4_err] = -99
data[u_err][np.abs(data[u]) == 99.0] = -99.0
data[g_err][np.abs(data[g]) == 99.0] = -99.0
data[r_err][np.abs(data[r]) == 99.0] = -99.0
data[i_err][np.abs(data[i]) == 99.0] = -99.0
data[z_err][np.abs(data[z]) == 99.0] = -99.0
data[Y_err][np.abs(data[Y]) == 99.0] = -99.0
data[J_err][np.abs(data[J]) == 99.0] = -99.0
data[H_err][np.abs(data[H]) == 99.0] = -99.0
data[K_err][np.abs(data[K]) == 99.0] = -99.0
data[u][np.abs(data[u]) == 99.0] = -99.0
data[g][np.abs(data[g]) == 99.0] = -99.0
data[r][np.abs(data[r]) == 99.0] = -99.0
data[i][np.abs(data[i]) == 99.0] = -99.0
data[z][np.abs(data[z]) == 99.0] = -99.0
data[Y][np.abs(data[Y]) == 99.0] = -99.0
data[J][np.abs(data[J]) == 99.0] = -99.0
data[H][np.abs(data[H]) == 99.0] = -99.0
data[K][np.abs(data[K]) == 99.0] = -99.0

# use photoz where specz is unavailable, or for stars
data[photoz_bpz][data[specz] > 0] = data[specz][data[specz] > 0]
data[photoz_eazy][data[specz] > 0] = data[specz][data[specz] > 0]

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

# LePhare thinks error bars = 0 means non-detection, so fix this
data[u_err][data[u_err] == 0.00] = 0.01
data[g_err][data[g_err] == 0.00] = 0.01
data[r_err][data[r_err] == 0.00] = 0.01
data[i_err][data[i_err] == 0.00] = 0.01
data[z_err][data[z_err] == 0.00] = 0.01
data[Y_err][data[Y_err] == 0.00] = 0.01
data[J_err][data[J_err] == 0.00] = 0.01
data[H_err][data[H_err] == 0.00] = 0.01
data[K_err][data[K_err] == 0.00] = 0.01

#fileout = file[:-16] + "_forlepharewithbpz_noobs.cat"
fileout = file[:-4] + "_forlepharewithbpzIRAC_noobs.cat"
str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t K_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[photoz_bpz]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')

#fileout = file[:-16] + "_forlepharewitheazy_noobs.cat"
fileout = file[:-4] + "_forlepharewitheazyIRAC_noobs.cat"
str = "ID \t u \t u_err \t g \t g_err \t r \t r_err \t i \t i_err \t z \t z_err \t Y \t Y_err \t J \t J_err \t H \t H_err \t K \t K_err \t ch1 \t ch1_err \t ch2 \t ch2_err \t ch3 \t ch3_err \t ch4 \t ch4_err \t context z-spec \t string"
dataout = np.c_[data[id],data[u],data[u_err],data[g],data[g_err],data[r],data[r_err],data[i],data[i_err],data[z],data[z_err],data[Y],data[Y_err],data[J],data[J_err],data[H],data[H_err],data[K],data[K_err],data[irac1],data[irac1_err],data[irac2],data[irac2_err],data[irac3],data[irac3_err],data[irac4],data[irac4_err],data[photoz_eazy]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t 8191 \t %.4f ')
