##########################
# Classifies galaxies and stars following the CFHTLenS method and available spectroscopic info
##########################

import numpy as np 

#file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephare.cat"
#file = "i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephare.cat"
#file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephare.cat"
file = "rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephare.cat"
limrad = 2.45

id = 8
#chi_gal_bpz = 58
#chi_star_bpz = 59
#chi_gal_eazy = 65
#chi_star_eazy = 66
chi_gal_bpz = 83 # for IRAC
chi_star_bpz = 84 # for IRAC
chi_gal_eazy = 90 # for IRAC
chi_star_eazy = 91 # for IRAC
radius = 6
mag = 4
spec = 40
classify1_bpz = 8 # first pass, CFHTLenS, fake index
classify2_bpz = 8 # second pass, including spectroscopic info, fake index
classify1_eazy = 8 # first pass, CFHTLenS, fake index
classify2_eazy = 8 # second pass, including spectroscopic info, fake index

data = np.loadtxt(file,usecols=(id,chi_gal_bpz,chi_star_bpz,chi_gal_eazy,chi_star_eazy,radius,mag,spec,classify1_bpz,classify2_bpz,classify1_eazy,classify2_eazy),unpack=True)

# indexing
id = 0
chi_gal_bpz = 1
chi_star_bpz = 2
chi_gal_eazy = 3
chi_star_eazy = 4
radius = 5
mag = 6
spec = 7
classify1_bpz = 8
classify2_bpz = 9
classify1_eazy = 10
classify2_eazy = 11

# CFHTLenS method: star=0 gal=1
for i in range(len(data[id])):
    if data[mag][i] < 21:
        if data[radius][i] <= limrad:
            data[classify1_bpz][i] = 0
            data[classify1_eazy][i] = 0
        else:
            data[classify1_bpz][i] = 1
            data[classify1_eazy][i] = 1
    if data[mag][i] >= 21 and data[mag][i] <= 23:
        if data[radius][i] <= limrad and data[chi_star_bpz][i] < 2.0 * data[chi_gal_bpz][i]:
            data[classify1_bpz][i] = 0
        else: data[classify1_bpz][i] = 1
        if data[radius][i] <= limrad and data[chi_star_eazy][i] < 2.0 * data[chi_gal_eazy][i]:
            data[classify1_eazy][i] = 0
        else: data[classify1_eazy][i] = 1
    if data[mag][i] > 23:
        data[classify1_bpz][i] = 1
        data[classify1_eazy][i] = 1

# include the spectroscopic knowledge:
# nospec star -3
# nospec gal 2
# spec star ok -1
# spec star wrong -2
# spec gal ok 0
# spec gal wrong 1
data[classify2_bpz][(data[classify1_bpz] == 0) & (data[spec] == -1)] = -3
data[classify2_bpz][(data[classify1_bpz] == 1) & (data[spec] == -1)] = 2
data[classify2_bpz][(data[classify1_bpz] == 0) & (data[spec] == -2)] = -1
data[classify2_bpz][(data[classify1_bpz] == 1) & (data[spec] == -2)] = -2
data[classify2_bpz][(data[classify1_bpz] == 1) & (data[spec] >= 0)] = 0
data[classify2_bpz][(data[classify1_bpz] == 0) & (data[spec] >= 0)] = 1

data[classify2_eazy][(data[classify1_eazy] == 0) & (data[spec] == -1)] = -3
data[classify2_eazy][(data[classify1_eazy] == 1) & (data[spec] == -1)] = 2
data[classify2_eazy][(data[classify1_eazy] == 0) & (data[spec] == -2)] = -1
data[classify2_eazy][(data[classify1_eazy] == 1) & (data[spec] == -2)] = -2
data[classify2_eazy][(data[classify1_eazy] == 1) & (data[spec] >= 0)] = 0
data[classify2_eazy][(data[classify1_eazy] == 0) & (data[spec] >= 0)] = 1

print "spectroscopic gal/stars:",len(data[spec][data[spec] >= 0]),"/",len(data[spec][data[spec] == -2])
print "CFHTLenS based gal/stars (BPZ):",len(data[classify1_bpz][data[classify1_bpz]==1]),"/",len(data[classify1_bpz][data[classify1_bpz]==0])," (",len(data[classify1_bpz][data[classify1_bpz]==1]) * 100.0 / (len(data[classify1_bpz][data[classify1_bpz]==1])+len(data[classify1_bpz][data[classify1_bpz]==0])),"%,",len(data[classify1_bpz][data[classify1_bpz]==0])*100.0/(len(data[classify1_bpz][data[classify1_bpz]==1])+len(data[classify1_bpz][data[classify1_bpz]==0])),"%)"
print "incorrectly identified gals (BPZ):",len(data[classify2_bpz][data[classify2_bpz] == 1]),"/",len(data[spec][data[spec]>= 0]),",", 100.0 * len(data[classify2_bpz][data[classify2_bpz] == 1]) / len(data[spec][data[spec]>= 0]), "+/-", np.sqrt(len(data[classify2_bpz][data[classify2_bpz] == 1])), "%"
print "incorrectly identified stars (BPZ):",len(data[classify2_bpz][data[classify2_bpz] == -2]),"/",len(data[spec][data[spec] == -2]),",", 100.0 * len(data[classify2_bpz][data[classify2_bpz] == -2]) / len(data[spec][data[spec] == -2]), "+/-", np.sqrt(len(data[classify2_bpz][data[classify2_bpz] == -2])), "%"

print "CFHTLenS based gal/stars (eazy):",len(data[classify1_eazy][data[classify1_eazy]==1]),"/",len(data[classify1_eazy][data[classify1_eazy]==0])," (",len(data[classify1_eazy][data[classify1_eazy]==1]) * 100.0 / (len(data[classify1_eazy][data[classify1_eazy]==1])+len(data[classify1_eazy][data[classify1_eazy]==0])),"%,",len(data[classify1_eazy][data[classify1_eazy]==0])*100.0/(len(data[classify1_eazy][data[classify1_eazy]==1])+len(data[classify1_eazy][data[classify1_eazy]==0])),"%)"
print "incorrectly identified gals (eazy):",len(data[classify2_eazy][data[classify2_eazy] == 1]),"/",len(data[spec][data[spec]>= 0]),",", 100.0 * len(data[classify2_eazy][data[classify2_eazy] == 1]) / len(data[spec][data[spec]>= 0]), "+/-", np.sqrt(len(data[classify2_eazy][data[classify2_eazy] == 1])), "%"
print "incorrectly identified stars (eazy):",len(data[classify2_eazy][data[classify2_eazy] == -2]),"/",len(data[spec][data[spec] == -2]),",", 100.0 * len(data[classify2_eazy][data[classify2_eazy] == -2]) / len(data[spec][data[spec] == -2]), "+/-", np.sqrt(len(data[classify2_eazy][data[classify2_eazy] == -2])), "%"

fileout = file[:-4] + "_classification.cat"
str = "ID class_bpz class_eazy"
dataout = np.c_[data[id],data[classify2_bpz],data[classify2_eazy]]
np.savetxt(fileout,dataout,header=str,fmt='%d \t %d \t %d \t')
