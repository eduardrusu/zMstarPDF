# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid

import numpy as np
#import scipy
#from scipy import special
#from astropy.io import fits
#from astropy.wcs import WCS
#from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.io import ascii
#from astropy.table import Table, Column
#import time
import matplotlib.pyplot as plt
#from numpy.random import normal
from scipy.stats.kde import gaussian_kde
from numpy import linspace

with open('fieldsforhist50try.lst') as f:
    listfields = f.readlines()

with open('fieldshistW1_50.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_50.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_50.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_50.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldsforhist75try.lst') as f:
    listfields = f.readlines()

with open('fieldshistW1_75.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W1" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW2_75.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W2" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW3_75.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W3" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

with open('fieldshistW4_75.lst', 'w') as outfile:
    for i in range(len(listfields)):
        if "W4" in [x[0:len(listfields[0])-1] for x in listfields][i]:
            with open([x[0:len(listfields[0])-1] for x in listfields][i]) as infile:
                outfile.write(infile.read())

q_gal_W1_50, q_oneoverr_W1_50, q_zweight_W1_50, sigma_zweight_W1_50, q_mass_W1_50, sigma_mass_W1_50, q_mass2_W1_50, sigma_mass2_W1_50, q_mass2rms_W1_50, sigma_mass2rms_W1_50, q_mass3_W1_50, sigma_mass3_W1_50, q_mass3rms_W1_50, sigma_mass3rms_W1_50, q_zoverr_W1_50, sigma_zoverr_W1_50, q_massoverr_W1_50, sigma_massoverr_W1_50, q_mass2overr_W1_50, sigma_mass2overr_W1_50, q_mass3overr_W1_50, sigma_mass3overr_W1_50, q_mass2overrrms_W1_50, sigma_mass2overrrms_W1_50, q_mass3overrrms_W1_50, sigma_mass3overrrms_W1_50, q_zmassoverr_W1_50, sigma_zmassoverr_W1_50, q_zmass2overr_W1_50, sigma_zmass2overr_W1_50 = np.loadtxt('fieldshistW1_50.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W2_50, q_oneoverr_W2_50, q_zweight_W2_50, sigma_zweight_W2_50, q_mass_W2_50, sigma_mass_W2_50, q_mass2_W2_50, sigma_mass2_W2_50, q_mass2rms_W2_50, sigma_mass2rms_W2_50, q_mass3_W2_50, sigma_mass3_W2_50, q_mass3rms_W2_50, sigma_mass3rms_W2_50, q_zoverr_W2_50, sigma_zoverr_W2_50, q_massoverr_W2_50, sigma_massoverr_W2_50, q_mass2overr_W2_50, sigma_mass2overr_W2_50, q_mass3overr_W2_50, sigma_mass3overr_W2_50, q_mass2overrrms_W2_50, sigma_mass2overrrms_W2_50, q_mass3overrrms_W2_50, sigma_mass3overrrms_W2_50, q_zmassoverr_W2_50, sigma_zmassoverr_W2_50, q_zmass2overr_W2_50, sigma_zmass2overr_W2_50 = np.loadtxt('fieldshistW1_50.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W3_50, q_oneoverr_W3_50, q_zweight_W3_50, sigma_zweight_W3_50, q_mass_W3_50, sigma_mass_W3_50, q_mass2_W3_50, sigma_mass2_W3_50, q_mass2rms_W3_50, sigma_mass2rms_W3_50, q_mass3_W3_50, sigma_mass3_W3_50, q_mass3rms_W3_50, sigma_mass3rms_W3_50, q_zoverr_W3_50, sigma_zoverr_W3_50, q_massoverr_W3_50, sigma_massoverr_W3_50, q_mass2overr_W3_50, sigma_mass2overr_W3_50, q_mass3overr_W3_50, sigma_mass3overr_W3_50, q_mass2overrrms_W3_50, sigma_mass2overrrms_W3_50, q_mass3overrrms_W3_50, sigma_mass3overrrms_W3_50, q_zmassoverr_W3_50, sigma_zmassoverr_W3_50, q_zmass2overr_W3_50, sigma_zmass2overr_W3_50 = np.loadtxt('fieldshistW3_50.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W4_50, q_oneoverr_W4_50, q_zweight_W4_50, sigma_zweight_W4_50, q_mass_W4_50, sigma_mass_W4_50, q_mass2_W4_50, sigma_mass2_W4_50, q_mass2rms_W4_50, sigma_mass2rms_W4_50, q_mass3_W4_50, sigma_mass3_W4_50, q_mass3rms_W4_50, sigma_mass3rms_W4_50, q_zoverr_W4_50, sigma_zoverr_W4_50, q_massoverr_W4_50, sigma_massoverr_W4_50, q_mass2overr_W4_50, sigma_mass2overr_W4_50, q_mass3overr_W4_50, sigma_mass3overr_W4_50, q_mass2overrrms_W4_50, sigma_mass2overrrms_W4_50, q_mass3overrrms_W4_50, sigma_mass3overrrms_W4_50, q_zmassoverr_W4_50, sigma_zmassoverr_W4_50, q_zmass2overr_W4_50, sigma_zmass2overr_W4_50 = np.loadtxt('fieldshistW4_50.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W1_75, q_oneoverr_W1_75, q_zweight_W1_75, sigma_zweight_W1_75, q_mass_W1_75, sigma_mass_W1_75, q_mass2_W1_75, sigma_mass2_W1_75, q_mass2rms_W1_75, sigma_mass2rms_W1_75, q_mass3_W1_75, sigma_mass3_W1_75, q_mass3rms_W1_75, sigma_mass3rms_W1_75, q_zoverr_W1_75, sigma_zoverr_W1_75, q_massoverr_W1_75, sigma_massoverr_W1_75, q_mass2overr_W1_75, sigma_mass2overr_W1_75, q_mass3overr_W1_75, sigma_mass3overr_W1_75, q_mass2overrrms_W1_75, sigma_mass2overrrms_W1_75, q_mass3overrrms_W1_75, sigma_mass3overrrms_W1_75, q_zmassoverr_W1_75, sigma_zmassoverr_W1_75, q_zmass2overr_W1_75, sigma_zmass2overr_W1_75 = np.loadtxt('fieldshistW1_75.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W2_75, q_oneoverr_W2_75, q_zweight_W2_75, sigma_zweight_W2_75, q_mass_W2_75, sigma_mass_W2_75, q_mass2_W2_75, sigma_mass2_W2_75, q_mass2rms_W2_75, sigma_mass2rms_W2_75, q_mass3_W2_75, sigma_mass3_W2_75, q_mass3rms_W2_75, sigma_mass3rms_W2_75, q_zoverr_W2_75, sigma_zoverr_W2_75, q_massoverr_W2_75, sigma_massoverr_W2_75, q_mass2overr_W2_75, sigma_mass2overr_W2_75, q_mass3overr_W2_75, sigma_mass3overr_W2_75, q_mass2overrrms_W2_75, sigma_mass2overrrms_W2_75, q_mass3overrrms_W2_75, sigma_mass3overrrms_W2_75, q_zmassoverr_W2_75, sigma_zmassoverr_W2_75, q_zmass2overr_W2_75, sigma_zmass2overr_W2_75 = np.loadtxt('fieldshistW1_75.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W3_75, q_oneoverr_W3_75, q_zweight_W3_75, sigma_zweight_W3_75, q_mass_W3_75, sigma_mass_W3_75, q_mass2_W3_75, sigma_mass2_W3_75, q_mass2rms_W3_75, sigma_mass2rms_W3_75, q_mass3_W3_75, sigma_mass3_W3_75, q_mass3rms_W3_75, sigma_mass3rms_W3_75, q_zoverr_W3_75, sigma_zoverr_W3_75, q_massoverr_W3_75, sigma_massoverr_W3_75, q_mass2overr_W3_75, sigma_mass2overr_W3_75, q_mass3overr_W3_75, sigma_mass3overr_W3_75, q_mass2overrrms_W3_75, sigma_mass2overrrms_W3_75, q_mass3overrrms_W3_75, sigma_mass3overrrms_W3_75, q_zmassoverr_W3_75, sigma_zmassoverr_W3_75, q_zmass2overr_W3_75, sigma_zmass2overr_W3_75 = np.loadtxt('fieldshistW3_75.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)

q_gal_W4_75, q_oneoverr_W4_75, q_zweight_W4_75, sigma_zweight_W4_75, q_mass_W4_75, sigma_mass_W4_75, q_mass2_W4_75, sigma_mass2_W4_75, q_mass2rms_W4_75, sigma_mass2rms_W4_75, q_mass3_W4_75, sigma_mass3_W4_75, q_mass3rms_W4_75, sigma_mass3rms_W4_75, q_zoverr_W4_75, sigma_zoverr_W4_75, q_massoverr_W4_75, sigma_massoverr_W4_75, q_mass2overr_W4_75, sigma_mass2overr_W4_75, q_mass3overr_W4_75, sigma_mass3overr_W4_75, q_mass2overrrms_W4_75, sigma_mass2overrrms_W4_75, q_mass3overrrms_W4_75, sigma_mass3overrrms_W4_75, q_zmassoverr_W4_75, sigma_zmassoverr_W4_75, q_zmass2overr_W4_75, sigma_zmass2overr_W4_75 = np.loadtxt('fieldshistW4_75.lst', usecols=(1, 3, 5, 7, 15, 17, 25, 27, 35, 37, 45, 47, 55, 57, 65, 67, 75, 77, 85, 87, 95, 97, 105, 107, 115, 117, 125, 127, 135, 137), unpack=True)


for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_zweight_W1_50[i], sigma_zweight_W1_50[i], 1000)
    q_zweight_W1_50=np.append(q_zweight_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass_W1_50[i], sigma_mass_W1_50[i], 1000)
    q_mass_W1_50=np.append(q_mass_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass2_W1_50[i], sigma_mass2_W1_50[i], 1000)
    q_mass2_W1_50=np.append(q_mass2_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass2rms_W1_50[i], sigma_mass2rms_W1_50[i], 1000)
    q_mass2rms_W1_50=np.append(q_mass2rms_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass3_W1_50[i], sigma_mass3_W1_50[i], 1000)
    q_mass3_W1_50=np.append(q_mass3_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass3rms_W1_50[i], sigma_mass3rms_W1_50[i], 1000)
    q_mass3rms_W1_50=np.append(q_mass3rms_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_zoverr_W1_50[i], sigma_zoverr_W1_50[i], 1000)
    q_zoverr_W1_50=np.append(q_zoverr_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_massoverr_W1_50[i], sigma_massoverr_W1_50[i], 1000)
    q_massoverr_W1_50=np.append(q_massoverr_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass2overr_W1_50[i], sigma_mass2overr_W1_50[i], 1000)
    q_mass2overr_W1_50=np.append(q_mass2overr_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass3overr_W1_50[i], sigma_mass3overr_W1_50[i], 1000)
    q_mass3overr_W1_50=np.append(q_mass3overr_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass2overrrms_W1_50[i], sigma_mass2overrrms_W1_50[i], 1000)
    q_mass2overrrms_W1_50=np.append(q_mass2overrrms_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_mass3overrrms_W1_50[i], sigma_mass3overrrms_W1_50[i], 1000)
    q_mass3overrrms_W1_50=np.append(q_mass3overrrms_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_zmassoverr_W1_50[i], sigma_zmassoverr_W1_50[i], 1000)
    q_zmassoverr_W1_50=np.append(q_zmassoverr_W1_50,gaussdist)

for i in range(len(q_gal_W1_50)):
    gaussdist = np.random.normal(q_zmass2overr_W1_50[i], sigma_zmass2overr_W1_50[i], 1000)
    q_zmass2overr_W1_50=np.append(q_zmass2overr_W1_50,gaussdist)



for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_zweight_W2_50[i], sigma_zweight_W2_50[i], 1000)
    q_zweight_W2_50=np.append(q_zweight_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass_W2_50[i], sigma_mass_W2_50[i], 1000)
    q_mass_W2_50=np.append(q_mass_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass2_W2_50[i], sigma_mass2_W2_50[i], 1000)
    q_mass2_W2_50=np.append(q_mass2_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass2rms_W2_50[i], sigma_mass2rms_W2_50[i], 1000)
    q_mass2rms_W2_50=np.append(q_mass2rms_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass3_W2_50[i], sigma_mass3_W2_50[i], 1000)
    q_mass3_W2_50=np.append(q_mass3_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass3rms_W2_50[i], sigma_mass3rms_W2_50[i], 1000)
    q_mass3rms_W2_50=np.append(q_mass3rms_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_zoverr_W2_50[i], sigma_zoverr_W2_50[i], 1000)
    q_zoverr_W2_50=np.append(q_zoverr_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_massoverr_W2_50[i], sigma_massoverr_W2_50[i], 1000)
    q_massoverr_W2_50=np.append(q_massoverr_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass2overr_W2_50[i], sigma_mass2overr_W2_50[i], 1000)
    q_mass2overr_W2_50=np.append(q_mass2overr_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass3overr_W2_50[i], sigma_mass3overr_W2_50[i], 1000)
    q_mass3overr_W2_50=np.append(q_mass3overr_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass2overrrms_W2_50[i], sigma_mass2overrrms_W2_50[i], 1000)
    q_mass2overrrms_W2_50=np.append(q_mass2overrrms_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_mass3overrrms_W2_50[i], sigma_mass3overrrms_W2_50[i], 1000)
    q_mass3overrrms_W2_50=np.append(q_mass3overrrms_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_zmassoverr_W2_50[i], sigma_zmassoverr_W2_50[i], 1000)
    q_zmassoverr_W2_50=np.append(q_zmassoverr_W2_50,gaussdist)

for i in range(len(q_gal_W2_50)):
    gaussdist = np.random.normal(q_zmass2overr_W2_50[i], sigma_zmass2overr_W2_50[i], 1000)
    q_zmass2overr_W2_50=np.append(q_zmass2overr_W2_50,gaussdist)



for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_zweight_W3_50[i], sigma_zweight_W3_50[i], 1000)
    q_zweight_W3_50=np.append(q_zweight_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass_W3_50[i], sigma_mass_W3_50[i], 1000)
    q_mass_W3_50=np.append(q_mass_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass2_W3_50[i], sigma_mass2_W3_50[i], 1000)
    q_mass2_W3_50=np.append(q_mass2_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass2rms_W3_50[i], sigma_mass2rms_W3_50[i], 1000)
    q_mass2rms_W3_50=np.append(q_mass2rms_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass3_W3_50[i], sigma_mass3_W3_50[i], 1000)
    q_mass3_W3_50=np.append(q_mass3_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass3rms_W3_50[i], sigma_mass3rms_W3_50[i], 1000)
    q_mass3rms_W3_50=np.append(q_mass3rms_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_zoverr_W3_50[i], sigma_zoverr_W3_50[i], 1000)
    q_zoverr_W3_50=np.append(q_zoverr_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_massoverr_W3_50[i], sigma_massoverr_W3_50[i], 1000)
    q_massoverr_W3_50=np.append(q_massoverr_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass2overr_W3_50[i], sigma_mass2overr_W3_50[i], 1000)
    q_mass2overr_W3_50=np.append(q_mass2overr_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass3overr_W3_50[i], sigma_mass3overr_W3_50[i], 1000)
    q_mass3overr_W3_50=np.append(q_mass3overr_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass2overrrms_W3_50[i], sigma_mass2overrrms_W3_50[i], 1000)
    q_mass2overrrms_W3_50=np.append(q_mass2overrrms_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_mass3overrrms_W3_50[i], sigma_mass3overrrms_W3_50[i], 1000)
    q_mass3overrrms_W3_50=np.append(q_mass3overrrms_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_zmassoverr_W3_50[i], sigma_zmassoverr_W3_50[i], 1000)
    q_zmassoverr_W3_50=np.append(q_zmassoverr_W3_50,gaussdist)

for i in range(len(q_gal_W3_50)):
    gaussdist = np.random.normal(q_zmass2overr_W3_50[i], sigma_zmass2overr_W3_50[i], 1000)
    q_zmass2overr_W3_50=np.append(q_zmass2overr_W3_50,gaussdist)




for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_zweight_W4_50[i], sigma_zweight_W4_50[i], 1000)
    q_zweight_W4_50=np.append(q_zweight_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass_W4_50[i], sigma_mass_W4_50[i], 1000)
    q_mass_W4_50=np.append(q_mass_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass2_W4_50[i], sigma_mass2_W4_50[i], 1000)
    q_mass2_W4_50=np.append(q_mass2_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass2rms_W4_50[i], sigma_mass2rms_W4_50[i], 1000)
    q_mass2rms_W4_50=np.append(q_mass2rms_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass3_W4_50[i], sigma_mass3_W4_50[i], 1000)
    q_mass3_W4_50=np.append(q_mass3_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass3rms_W4_50[i], sigma_mass3rms_W4_50[i], 1000)
    q_mass3rms_W4_50=np.append(q_mass3rms_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_zoverr_W4_50[i], sigma_zoverr_W4_50[i], 1000)
    q_zoverr_W4_50=np.append(q_zoverr_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_massoverr_W4_50[i], sigma_massoverr_W4_50[i], 1000)
    q_massoverr_W4_50=np.append(q_massoverr_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass2overr_W4_50[i], sigma_mass2overr_W4_50[i], 1000)
    q_mass2overr_W4_50=np.append(q_mass2overr_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass3overr_W4_50[i], sigma_mass3overr_W4_50[i], 1000)
    q_mass3overr_W4_50=np.append(q_mass3overr_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass2overrrms_W4_50[i], sigma_mass2overrrms_W4_50[i], 1000)
    q_mass2overrrms_W4_50=np.append(q_mass2overrrms_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_mass3overrrms_W4_50[i], sigma_mass3overrrms_W4_50[i], 1000)
    q_mass3overrrms_W4_50=np.append(q_mass3overrrms_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_zmassoverr_W4_50[i], sigma_zmassoverr_W4_50[i], 1000)
    q_zmassoverr_W4_50=np.append(q_zmassoverr_W4_50,gaussdist)

for i in range(len(q_gal_W4_50)):
    gaussdist = np.random.normal(q_zmass2overr_W4_50[i], sigma_zmass2overr_W4_50[i], 1000)
    q_zmass2overr_W4_50=np.append(q_zmass2overr_W4_50,gaussdist)




for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_zweight_W1_75[i], sigma_zweight_W1_75[i], 1000)
    q_zweight_W1_75=np.append(q_zweight_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass_W1_75[i], sigma_mass_W1_75[i], 1000)
    q_mass_W1_75=np.append(q_mass_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass2_W1_75[i], sigma_mass2_W1_75[i], 1000)
    q_mass2_W1_75=np.append(q_mass2_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass2rms_W1_75[i], sigma_mass2rms_W1_75[i], 1000)
    q_mass2rms_W1_75=np.append(q_mass2rms_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass3_W1_75[i], sigma_mass3_W1_75[i], 1000)
    q_mass3_W1_75=np.append(q_mass3_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass3rms_W1_75[i], sigma_mass3rms_W1_75[i], 1000)
    q_mass3rms_W1_75=np.append(q_mass3rms_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_zoverr_W1_75[i], sigma_zoverr_W1_75[i], 1000)
    q_zoverr_W1_75=np.append(q_zoverr_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_massoverr_W1_75[i], sigma_massoverr_W1_75[i], 1000)
    q_massoverr_W1_75=np.append(q_massoverr_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass2overr_W1_75[i], sigma_mass2overr_W1_75[i], 1000)
    q_mass2overr_W1_75=np.append(q_mass2overr_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass3overr_W1_75[i], sigma_mass3overr_W1_75[i], 1000)
    q_mass3overr_W1_75=np.append(q_mass3overr_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass2overrrms_W1_75[i], sigma_mass2overrrms_W1_75[i], 1000)
    q_mass2overrrms_W1_75=np.append(q_mass2overrrms_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_mass3overrrms_W1_75[i], sigma_mass3overrrms_W1_75[i], 1000)
    q_mass3overrrms_W1_75=np.append(q_mass3overrrms_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_zmassoverr_W1_75[i], sigma_zmassoverr_W1_75[i], 1000)
    q_zmassoverr_W1_75=np.append(q_zmassoverr_W1_75,gaussdist)

for i in range(len(q_gal_W1_75)):
    gaussdist = np.random.normal(q_zmass2overr_W1_75[i], sigma_zmass2overr_W1_75[i], 1000)
    q_zmass2overr_W1_75=np.append(q_zmass2overr_W1_75,gaussdist)




for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_zweight_W2_75[i], sigma_zweight_W2_75[i], 1000)
    q_zweight_W2_75=np.append(q_zweight_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass_W2_75[i], sigma_mass_W2_75[i], 1000)
    q_mass_W2_75=np.append(q_mass_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass2_W2_75[i], sigma_mass2_W2_75[i], 1000)
    q_mass2_W2_75=np.append(q_mass2_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass2rms_W2_75[i], sigma_mass2rms_W2_75[i], 1000)
    q_mass2rms_W2_75=np.append(q_mass2rms_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass3_W2_75[i], sigma_mass3_W2_75[i], 1000)
    q_mass3_W2_75=np.append(q_mass3_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass3rms_W2_75[i], sigma_mass3rms_W2_75[i], 1000)
    q_mass3rms_W2_75=np.append(q_mass3rms_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_zoverr_W2_75[i], sigma_zoverr_W2_75[i], 1000)
    q_zoverr_W2_75=np.append(q_zoverr_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_massoverr_W2_75[i], sigma_massoverr_W2_75[i], 1000)
    q_massoverr_W2_75=np.append(q_massoverr_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass2overr_W2_75[i], sigma_mass2overr_W2_75[i], 1000)
    q_mass2overr_W2_75=np.append(q_mass2overr_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass3overr_W2_75[i], sigma_mass3overr_W2_75[i], 1000)
    q_mass3overr_W2_75=np.append(q_mass3overr_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass2overrrms_W2_75[i], sigma_mass2overrrms_W2_75[i], 1000)
    q_mass2overrrms_W2_75=np.append(q_mass2overrrms_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_mass3overrrms_W2_75[i], sigma_mass3overrrms_W2_75[i], 1000)
    q_mass3overrrms_W2_75=np.append(q_mass3overrrms_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_zmassoverr_W2_75[i], sigma_zmassoverr_W2_75[i], 1000)
    q_zmassoverr_W2_75=np.append(q_zmassoverr_W2_75,gaussdist)

for i in range(len(q_gal_W2_75)):
    gaussdist = np.random.normal(q_zmass2overr_W2_75[i], sigma_zmass2overr_W2_75[i], 1000)
    q_zmass2overr_W2_75=np.append(q_zmass2overr_W2_75,gaussdist)




for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_zweight_W3_75[i], sigma_zweight_W3_75[i], 1000)
    q_zweight_W3_75=np.append(q_zweight_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass_W3_75[i], sigma_mass_W3_75[i], 1000)
    q_mass_W3_75=np.append(q_mass_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass2_W3_75[i], sigma_mass2_W3_75[i], 1000)
    q_mass2_W3_75=np.append(q_mass2_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass2rms_W3_75[i], sigma_mass2rms_W3_75[i], 1000)
    q_mass2rms_W3_75=np.append(q_mass2rms_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass3_W3_75[i], sigma_mass3_W3_75[i], 1000)
    q_mass3_W3_75=np.append(q_mass3_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass3rms_W3_75[i], sigma_mass3rms_W3_75[i], 1000)
    q_mass3rms_W3_75=np.append(q_mass3rms_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_zoverr_W3_75[i], sigma_zoverr_W3_75[i], 1000)
    q_zoverr_W3_75=np.append(q_zoverr_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_massoverr_W3_75[i], sigma_massoverr_W3_75[i], 1000)
    q_massoverr_W3_75=np.append(q_massoverr_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass2overr_W3_75[i], sigma_mass2overr_W3_75[i], 1000)
    q_mass2overr_W3_75=np.append(q_mass2overr_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass3overr_W3_75[i], sigma_mass3overr_W3_75[i], 1000)
    q_mass3overr_W3_75=np.append(q_mass3overr_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass2overrrms_W3_75[i], sigma_mass2overrrms_W3_75[i], 1000)
    q_mass2overrrms_W3_75=np.append(q_mass2overrrms_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_mass3overrrms_W3_75[i], sigma_mass3overrrms_W3_75[i], 1000)
    q_mass3overrrms_W3_75=np.append(q_mass3overrrms_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_zmassoverr_W3_75[i], sigma_zmassoverr_W3_75[i], 1000)
    q_zmassoverr_W3_75=np.append(q_zmassoverr_W3_75,gaussdist)

for i in range(len(q_gal_W3_75)):
    gaussdist = np.random.normal(q_zmass2overr_W3_75[i], sigma_zmass2overr_W3_75[i], 1000)
    q_zmass2overr_W3_75=np.append(q_zmass2overr_W3_75,gaussdist)




for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_zweight_W4_75[i], sigma_zweight_W4_75[i], 1000)
    q_zweight_W4_75=np.append(q_zweight_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass_W4_75[i], sigma_mass_W4_75[i], 1000)
    q_mass_W4_75=np.append(q_mass_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass2_W4_75[i], sigma_mass2_W4_75[i], 1000)
    q_mass2_W4_75=np.append(q_mass2_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass2rms_W4_75[i], sigma_mass2rms_W4_75[i], 1000)
    q_mass2rms_W4_75=np.append(q_mass2rms_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass3_W4_75[i], sigma_mass3_W4_75[i], 1000)
    q_mass3_W4_75=np.append(q_mass3_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass3rms_W4_75[i], sigma_mass3rms_W4_75[i], 1000)
    q_mass3rms_W4_75=np.append(q_mass3rms_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_zoverr_W4_75[i], sigma_zoverr_W4_75[i], 1000)
    q_zoverr_W4_75=np.append(q_zoverr_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_massoverr_W4_75[i], sigma_massoverr_W4_75[i], 1000)
    q_massoverr_W4_75=np.append(q_massoverr_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass2overr_W4_75[i], sigma_mass2overr_W4_75[i], 1000)
    q_mass2overr_W4_75=np.append(q_mass2overr_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass3overr_W4_75[i], sigma_mass3overr_W4_75[i], 1000)
    q_mass3overr_W4_75=np.append(q_mass3overr_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass2overrrms_W4_75[i], sigma_mass2overrrms_W4_75[i], 1000)
    q_mass2overrrms_W4_75=np.append(q_mass2overrrms_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_mass3overrrms_W4_75[i], sigma_mass3overrrms_W4_75[i], 1000)
    q_mass3overrrms_W4_75=np.append(q_mass3overrrms_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_zmassoverr_W4_75[i], sigma_zmassoverr_W4_75[i], 1000)
    q_zmassoverr_W4_75=np.append(q_zmassoverr_W4_75,gaussdist)

for i in range(len(q_gal_W4_75)):
    gaussdist = np.random.normal(q_zmass2overr_W4_75[i], sigma_zmass2overr_W4_75[i], 1000)
    q_zmass2overr_W4_75=np.append(q_zmass2overr_W4_75,gaussdist)

plt.suptitle(r'HE0435 weight histogram test W1-W4', fontsize=10, y=0.998)

gauss_q_gal_W1_50 = gaussian_kde(q_gal_W1_50)
gauss_q_gal_W2_50 = gaussian_kde(q_gal_W2_50)
gauss_q_gal_W3_50 = gaussian_kde(q_gal_W3_50)
gauss_q_gal_W4_50 = gaussian_kde(q_gal_W4_50)
gauss_q_gal_W1_75 = gaussian_kde(q_gal_W1_75)
gauss_q_gal_W2_75 = gaussian_kde(q_gal_W2_75)
gauss_q_gal_W3_75 = gaussian_kde(q_gal_W3_75)
gauss_q_gal_W4_75 = gaussian_kde(q_gal_W4_75)

x = linspace(0,2,500)

plt.subplot(451)
#plt.hist(q_gal_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#n_q_gal_W1_50, bins_q_gal_W1_50, patches = plt.hist(q_gal_W1_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_gal_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_gal_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_gal_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_gal_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_gal_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_gal_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_gal_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.plot(x,gauss_q_gal_W1_50(x),'b', linewidth=0.5)
#ave1 = (np.abs(gauss_q_gal_W1_50(x) - np.average(gauss_q_gal_W1_50(x)))).argmin()
#second = np.delete(gauss_q_gal_W1_50(x),(np.abs(gauss_q_gal_W1_50(x) - np.average(gauss_q_gal_W1_50(x)))).argmin())
#ave2 = (np.abs(second - np.average(gauss_q_gal_W1_50(x)))).argmin()
#med1 = (np.abs(gauss_q_gal_W1_50(x) - np.median(gauss_q_gal_W1_50(x)))).argmin()
#second = np.delete(gauss_q_gal_W1_50(x),(np.abs(gauss_q_gal_W1_50(x) - np.median(gauss_q_gal_W1_50(x)))).argmin())
#med2 = (np.abs(second - np.median(gauss_q_gal_W1_50(x)))).argmin()
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_gal_W1_50(x))],np.average(q_gal_W1_50),np.median(q_gal_W1_50))
plt.text(0.7, 2.5, s, fontsize=5, color='b')
plt.plot(x,gauss_q_gal_W2_50(x),'g', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_gal_W2_50(x))],np.average(q_gal_W2_50),np.median(q_gal_W2_50))
plt.text(0.7, 2.0, s, fontsize=5, color='g')
plt.plot(x,gauss_q_gal_W3_50(x),'r', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_gal_W3_50(x))],np.average(q_gal_W3_50),np.median(q_gal_W3_50))
plt.text(0.7, 1.5, s, fontsize=5, color='r')
plt.plot(x,gauss_q_gal_W4_50(x),'k', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_gal_W4_50(x))],np.average(q_gal_W4_50),np.median(q_gal_W4_50))
plt.text(0.7, 1.0, s, fontsize=5, color='k')
#s = "%.3f" % bins_q_gal_W2_50[np.argmax(n_q_gal_W2_50)]
#plt.text(1.3, 2.0, s, fontsize=5, color='g')
plt.xlabel(r'$\zeta_{gal}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

x = linspace(0,2,500)

gauss_q_oneoverr_W1_50 = gaussian_kde(q_oneoverr_W1_50)
gauss_q_oneoverr_W2_50 = gaussian_kde(q_oneoverr_W2_50)
gauss_q_oneoverr_W3_50 = gaussian_kde(q_oneoverr_W3_50)
gauss_q_oneoverr_W4_50 = gaussian_kde(q_oneoverr_W4_50)
gauss_q_oneoverr_W1_75 = gaussian_kde(q_oneoverr_W1_75)
gauss_q_oneoverr_W2_75 = gaussian_kde(q_oneoverr_W2_75)
gauss_q_oneoverr_W3_75 = gaussian_kde(q_oneoverr_W3_75)
gauss_q_oneoverr_W4_75 = gaussian_kde(q_oneoverr_W4_75)

plt.subplot(452)
#plt.hist(q_oneoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_oneoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_oneoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_oneoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
#plt.hist(q_oneoverr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_oneoverr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_oneoverr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
#plt.hist(q_oneoverr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.plot(x,gauss_q_oneoverr_W1_50(x),'b', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_oneoverr_W1_50(x))],np.average(q_oneoverr_W1_50),np.median(q_oneoverr_W1_50))
plt.text(0.6, 2.5, s, fontsize=5, color='b')
plt.plot(x,gauss_q_oneoverr_W2_50(x),'g', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_oneoverr_W2_50(x))],np.average(q_oneoverr_W2_50),np.median(q_oneoverr_W2_50))
plt.text(0.6, 2.0, s, fontsize=5, color='g')
plt.plot(x,gauss_q_oneoverr_W3_50(x),'r', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_oneoverr_W3_50(x))],np.average(q_oneoverr_W3_50),np.median(q_oneoverr_W3_50))
plt.text(0.6, 1.5, s, fontsize=5, color='r')
plt.plot(x,gauss_q_oneoverr_W4_50(x),'k', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_oneoverr_W4_50(x))],np.average(q_oneoverr_W4_50),np.median(q_oneoverr_W4_50))
plt.text(0.6, 1.0, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{1}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(453)
n_q_zweight_W1_50, bins_q_zweight_W1_50, patches = plt.hist(q_zweight_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zweight_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zweight_W2_50, bins_q_zweight_W2_50, patches = plt.hist(q_zweight_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zweight_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zweight_W3_50, bins_q_zweight_W3_50, patches = plt.hist(q_zweight_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zweight_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zweight_W4_50, bins_q_zweight_W4_50, patches = plt.hist(q_zweight_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zweight_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zweight_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zweight_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zweight_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zweight_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_zweight_W1_50[np.argmax(n_q_zweight_W1_50)],np.average(q_zweight_W1_50),np.median(q_zweight_W1_50))
plt.text(0.6, 2.0, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_zweight_W2_50[np.argmax(n_q_zweight_W2_50)],np.average(q_zweight_W2_50),np.median(q_zweight_W2_50))
plt.text(0.6, 1.5, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_zweight_W3_50[np.argmax(n_q_zweight_W3_50)],np.average(q_zweight_W3_50),np.median(q_zweight_W3_50))
plt.text(0.6, 1.0, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_zweight_W4_50[np.argmax(n_q_zweight_W4_50)],np.average(q_zweight_W4_50),np.median(q_zweight_W4_50))
plt.text(0.6, 0.5, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{z}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(454)
n_q_mass_W1_50, bins_q_mass_W1_50, patches = plt.hist(q_mass_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass_W2_50, bins_q_mass_W2_50, patches = plt.hist(q_mass_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass_W3_50, bins_q_mass_W3_50, patches = plt.hist(q_mass_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass_W4_50, bins_q_mass_W4_50, patches = plt.hist(q_mass_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass_W1_50[np.argmax(n_q_mass_W1_50)],np.average(q_mass_W1_50),np.median(q_mass_W1_50))
plt.text(0.8, 1.4, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass_W2_50[np.argmax(n_q_mass_W2_50)],np.average(q_mass_W2_50),np.median(q_mass_W2_50))
plt.text(0.8, 1.2, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass_W3_50[np.argmax(n_q_mass_W3_50)],np.average(q_mass_W3_50),np.median(q_mass_W3_50))
plt.text(0.8, 1.0, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass_W4_50[np.argmax(n_q_mass_W4_50)],np.average(q_mass_W4_50),np.median(q_mass_W4_50))
plt.text(0.8, 0.8, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{M}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(456)
n_q_mass2_W1_50, bins_q_mass2_W1_50, patches = plt.hist(q_mass2_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_mass2_W2_50, bins_q_mass2_W2_50, patches = plt.hist(q_mass2_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_mass2_W3_50, bins_q_mass2_W3_50, patches = plt.hist(q_mass2_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_mass2_W4_50, bins_q_mass2_W4_50, patches = plt.hist(q_mass2_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_mass2_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_mass2_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_mass2_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_mass2_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_mass2_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_mass2_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_mass2_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_mass2_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2_W1_50[np.argmax(n_q_mass2_W1_50)],np.average(q_mass2_W1_50),np.median(q_mass2_W1_50))
plt.text(1.0, 0.65, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2_W2_50[np.argmax(n_q_mass2_W2_50)],np.average(q_mass2_W2_50),np.median(q_mass2_W2_50))
plt.text(1.0, 0.50, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2_W3_50[np.argmax(n_q_mass2_W3_50)],np.average(q_mass2_W3_50),np.median(q_mass2_W3_50))
plt.text(1.0, 0.35, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2_W4_50[np.argmax(n_q_mass2_W4_50)],np.average(q_mass2_W4_50),np.median(q_mass2_W4_50))
plt.text(1.0, 0.2, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{M^2}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(457)
n_q_mass2rms_W1_50, bins_q_mass2rms_W1_50, patches = plt.hist(q_mass2rms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass2rms_W2_50, bins_q_mass2rms_W2_50, patches = plt.hist(q_mass2rms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass2rms_W3_50, bins_q_mass2rms_W3_50, patches = plt.hist(q_mass2rms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass2rms_W4_50, bins_q_mass2rms_W4_50, patches = plt.hist(q_mass2rms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2rms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2rms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2rms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2rms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2rms_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass2rms_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass2rms_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass2rms_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2rms_W1_50[np.argmax(n_q_mass2rms_W1_50)],np.average(q_mass2rms_W1_50),np.median(q_mass2rms_W1_50))
plt.text(1.0, 1.0, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2rms_W2_50[np.argmax(n_q_mass2rms_W2_50)],np.average(q_mass2rms_W2_50),np.median(q_mass2rms_W2_50))
plt.text(1.0, 0.8, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2rms_W3_50[np.argmax(n_q_mass2rms_W3_50)],np.average(q_mass2rms_W3_50),np.median(q_mass2rms_W3_50))
plt.text(1.0, 0.6, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2rms_W4_50[np.argmax(n_q_mass2rms_W4_50)],np.average(q_mass2rms_W4_50),np.median(q_mass2rms_W4_50))
plt.text(1.0, 0.4, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{M^2_{rms}}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(458)
n_q_mass3_W1_50, bins_q_mass3_W1_50, patches = plt.hist(q_mass3_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass3_W2_50, bins_q_mass3_W2_50, patches = plt.hist(q_mass3_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass3_W3_50, bins_q_mass3_W3_50, patches = plt.hist(q_mass3_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass3_W4_50, bins_q_mass3_W4_50, patches = plt.hist(q_mass3_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass3_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass3_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass3_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass3_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass3_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass3_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass3_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass3_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3_W1_50[np.argmax(n_q_mass3_W1_50)],np.average(q_mass3_W1_50),np.median(q_mass3_W1_50))
plt.text(1.0, 0.8, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3_W2_50[np.argmax(n_q_mass3_W2_50)],np.average(q_mass3_W2_50),np.median(q_mass3_W2_50))
plt.text(1.0, 0.6, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3_W3_50[np.argmax(n_q_mass3_W3_50)],np.average(q_mass3_W3_50),np.median(q_mass3_W3_50))
plt.text(1.0, 0.4, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3_W4_50[np.argmax(n_q_mass3_W4_50)],np.average(q_mass3_W4_50),np.median(q_mass3_W4_50))
plt.text(1.0, 0.2, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{M^3}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(459)
n_q_mass3rms_W1_50, bins_q_mass3rms_W1_50, patches = plt.hist(q_mass3rms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass3rms_W2_50, bins_q_mass3rms_W2_50, patches = plt.hist(q_mass3rms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass3rms_W3_50, bins_q_mass3rms_W3_50, patches = plt.hist(q_mass3rms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_mass3rms_W4_50, bins_q_mass3rms_W4_50, patches = plt.hist(q_mass3rms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3rms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3rms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3rms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3rms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3rms_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass3rms_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass3rms_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_mass3rms_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3rms_W1_50[np.argmax(n_q_mass3rms_W1_50)],np.average(q_mass3rms_W1_50),np.median(q_mass3rms_W1_50))
plt.text(1.0, 1.2, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3rms_W2_50[np.argmax(n_q_mass3rms_W2_50)],np.average(q_mass3rms_W2_50),np.median(q_mass3rms_W2_50))
plt.text(1.0, 1.0, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3rms_W3_50[np.argmax(n_q_mass3rms_W3_50)],np.average(q_mass3rms_W3_50),np.median(q_mass3rms_W3_50))
plt.text(1.0, 0.8, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3rms_W4_50[np.argmax(n_q_mass3rms_W4_50)],np.average(q_mass3rms_W4_50),np.median(q_mass3rms_W4_50))
plt.text(1.0, 0.6, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_{M^3_{rms}}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,11)
n_q_zoverr_W1_50, bins_q_zoverr_W1_50, patches = plt.hist(q_zoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zoverr_W2_50, bins_q_zoverr_W2_50, patches = plt.hist(q_zoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zoverr_W3_50, bins_q_zoverr_W3_50, patches = plt.hist(q_zoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
n_q_zoverr_W4_50, bins_q_zoverr_W4_50, patches = plt.hist(q_zoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 2])
plt.hist(q_zoverr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zoverr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zoverr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
plt.hist(q_zoverr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 2], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_zoverr_W1_50[np.argmax(n_q_zoverr_W1_50)],np.average(q_zoverr_W1_50),np.median(q_zoverr_W1_50))
plt.text(0.5, 1.6, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_zoverr_W2_50[np.argmax(n_q_zoverr_W2_50)],np.average(q_zoverr_W2_50),np.median(q_zoverr_W2_50))
plt.text(0.5, 1.3, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_zoverr_W3_50[np.argmax(n_q_zoverr_W3_50)],np.average(q_zoverr_W3_50),np.median(q_zoverr_W3_50))
plt.text(0.5, 1.0, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_zoverr_W4_50[np.argmax(n_q_zoverr_W4_50)],np.average(q_zoverr_W4_50),np.median(q_zoverr_W4_50))
plt.text(0.5, 0.7, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{z}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,12)
n_q_massoverr_W1_50, bins_q_massoverr_W1_50, patches = plt.hist(q_massoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_massoverr_W2_50, bins_q_massoverr_W2_50, patches = plt.hist(q_massoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_massoverr_W3_50, bins_q_massoverr_W3_50, patches = plt.hist(q_massoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
n_q_massoverr_W4_50, bins_q_massoverr_W4_50, patches = plt.hist(q_massoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_massoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_massoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_massoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_massoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_massoverr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_massoverr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_massoverr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
plt.hist(q_massoverr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_massoverr_W1_50[np.argmax(n_q_massoverr_W1_50)],np.average(q_massoverr_W1_50),np.median(q_massoverr_W1_50))
plt.text(1.0, 1.0, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_massoverr_W2_50[np.argmax(n_q_massoverr_W2_50)],np.average(q_massoverr_W2_50),np.median(q_massoverr_W2_50))
plt.text(1.0, 0.8, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_massoverr_W3_50[np.argmax(n_q_massoverr_W3_50)],np.average(q_massoverr_W3_50),np.median(q_massoverr_W3_50))
plt.text(1.0, 0.6, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_massoverr_W4_50[np.argmax(n_q_massoverr_W4_50)],np.average(q_massoverr_W4_50),np.median(q_massoverr_W4_50))
plt.text(1.0, 0.4, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{M}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,13)
n_q_mass2overr_W1_50, bins_q_mass2overr_W1_50, patches = plt.hist(q_mass2overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass2overr_W2_50, bins_q_mass2overr_W2_50, patches = plt.hist(q_mass2overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass2overr_W3_50, bins_q_mass2overr_W3_50, patches = plt.hist(q_mass2overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
n_q_mass2overr_W4_50, bins_q_mass2overr_W4_50, patches = plt.hist(q_mass2overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass2overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass2overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass2overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass2overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 5])
plt.hist(q_mass2overr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass2overr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass2overr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
plt.hist(q_mass2overr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 5], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2overr_W1_50[np.argmax(n_q_mass2overr_W1_50)],np.average(q_mass2overr_W1_50),np.median(q_mass2overr_W1_50))
plt.text(1.5, 0.7, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2overr_W2_50[np.argmax(n_q_mass2overr_W2_50)],np.average(q_mass2overr_W2_50),np.median(q_mass2overr_W2_50))
plt.text(1.5, 0.5, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2overr_W3_50[np.argmax(n_q_mass2overr_W3_50)],np.average(q_mass2overr_W3_50),np.median(q_mass2overr_W3_50))
plt.text(1.5, 0.3, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass2overr_W4_50[np.argmax(n_q_mass2overr_W4_50)],np.average(q_mass2overr_W4_50),np.median(q_mass2overr_W4_50))
plt.text(1.5, 0.1, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{M^2}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,14)
n_q_mass3overr_W1_50, bins_q_mass3overr_W1_50, patches = plt.hist(q_mass3overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_mass3overr_W2_50, bins_q_mass3overr_W2_50, patches = plt.hist(q_mass3overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_mass3overr_W3_50, bins_q_mass3overr_W3_50, patches = plt.hist(q_mass3overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_mass3overr_W4_50, bins_q_mass3overr_W4_50, patches = plt.hist(q_mass3overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_mass3overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_mass3overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_mass3overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_mass3overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_mass3overr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_mass3overr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_mass3overr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_mass3overr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3overr_W1_50[np.argmax(n_q_mass3overr_W1_50)],np.average(q_mass3overr_W1_50),np.median(q_mass3overr_W1_50))
plt.text(1.5, 0.7, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3overr_W2_50[np.argmax(n_q_mass3overr_W2_50)],np.average(q_mass3overr_W2_50),np.median(q_mass3overr_W2_50))
plt.text(1.5, 0.5, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3overr_W3_50[np.argmax(n_q_mass3overr_W3_50)],np.average(q_mass3overr_W3_50),np.median(q_mass3overr_W3_50))
plt.text(1.5, 0.3, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_mass3overr_W4_50[np.argmax(n_q_mass3overr_W4_50)],np.average(q_mass3overr_W4_50),np.median(q_mass3overr_W4_50))
plt.text(1.5, 0.1, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{M^3}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

x = linspace(0,3,500)

gauss_q_mass2overrrms_W1_50 = gaussian_kde(q_mass2overrrms_W1_50)
gauss_q_mass2overrrms_W2_50 = gaussian_kde(q_mass2overrrms_W2_50)
gauss_q_mass2overrrms_W3_50 = gaussian_kde(q_mass2overrrms_W3_50)
gauss_q_mass2overrrms_W4_50 = gaussian_kde(q_mass2overrrms_W4_50)
gauss_q_mass2overrrms_W1_75 = gaussian_kde(q_mass2overrrms_W1_75)
gauss_q_mass2overrrms_W2_75 = gaussian_kde(q_mass2overrrms_W2_75)
gauss_q_mass2overrrms_W3_75 = gaussian_kde(q_mass2overrrms_W3_75)
gauss_q_mass2overrrms_W4_75 = gaussian_kde(q_mass2overrrms_W4_75)

plt.subplot(4,5,16)
plt.hist(q_mass2overrrms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2overrrms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2overrrms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass2overrrms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
#plt.hist(q_mass2overrrms_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass2overrrms_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass2overrrms_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass2overrrms_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.plot(x,gauss_q_mass2overrrms_W1_50(x),'b', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass2overrrms_W1_50(x))],np.average(q_mass2overrrms_W1_50),np.median(q_mass2overrrms_W1_50))
plt.text(1.0, 1.4, s, fontsize=5, color='b')
#plt.plot(x,gauss_q_mass2overrrms_W2_50(x),'g', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass2overrrms_W2_50(x))],np.average(q_mass2overrrms_W2_50),np.median(q_mass2overrrms_W2_50))
plt.text(1.0, 1.2, s, fontsize=5, color='g')
#plt.plot(x,gauss_q_mass2overrrms_W3_50(x),'r', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass2overrrms_W3_50(x))],np.average(q_mass2overrrms_W3_50),np.median(q_mass2overrrms_W3_50))
plt.text(1.0, 1.0, s, fontsize=5, color='r')
#plt.plot(x,gauss_q_mass2overrrms_W4_50(x),'k', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass2overrrms_W4_50(x))],np.average(q_mass2overrrms_W4_50),np.median(q_mass2overrrms_W4_50))
plt.text(1.0, 0.8, s, fontsize=5, color='k')
plt.xlabel(r'${\zeta_\frac{M_{rms}^2}{r}}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

x = linspace(0,3,500)

gauss_q_mass3overrrms_W1_50 = gaussian_kde(q_mass3overrrms_W1_50)
gauss_q_mass3overrrms_W2_50 = gaussian_kde(q_mass3overrrms_W2_50)
gauss_q_mass3overrrms_W3_50 = gaussian_kde(q_mass3overrrms_W3_50)
gauss_q_mass3overrrms_W4_50 = gaussian_kde(q_mass3overrrms_W4_50)
gauss_q_mass3overrrms_W1_75 = gaussian_kde(q_mass3overrrms_W1_75)
gauss_q_mass3overrrms_W2_75 = gaussian_kde(q_mass3overrrms_W2_75)
gauss_q_mass3overrrms_W3_75 = gaussian_kde(q_mass3overrrms_W3_75)
gauss_q_mass3overrrms_W4_75 = gaussian_kde(q_mass3overrrms_W4_75)

plt.subplot(4,5,17)
plt.hist(q_mass3overrrms_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3overrrms_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3overrrms_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
plt.hist(q_mass3overrrms_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 3])
#plt.hist(q_mass3overrrms_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass3overrrms_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass3overrrms_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.hist(q_mass3overrrms_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 3], linestyle='dotted')
#plt.plot(x,gauss_q_mass3overrrms_W1_50(x),'b', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass3overrrms_W1_50(x))],np.average(q_mass3overrrms_W1_50),np.median(q_mass3overrrms_W1_50))
plt.text(1.0, 1.4, s, fontsize=5, color='b')
#plt.plot(x,gauss_q_mass3overrrms_W2_50(x),'g', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass3overrrms_W2_50(x))],np.average(q_mass3overrrms_W2_50),np.median(q_mass3overrrms_W2_50))
plt.text(1.0, 1.2, s, fontsize=5, color='g')
#plt.plot(x,gauss_q_mass3overrrms_W3_50(x),'r', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass3overrrms_W3_50(x))],np.average(q_mass3overrrms_W3_50),np.median(q_mass3overrrms_W3_50))
plt.text(1.0, 1.0, s, fontsize=5, color='r')
#plt.plot(x,gauss_q_mass3overrrms_W4_50(x),'k', linewidth=0.5)
s = "%.3f,%.3f,%.3f" % (x[np.argmax(gauss_q_mass3overrrms_W4_50(x))],np.average(q_mass3overrrms_W4_50),np.median(q_mass3overrrms_W4_50))
plt.text(1.0, 0.8, s, fontsize=5, color='k')
plt.xlabel(r'${\zeta_\frac{M_{rms}^3}{r}}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,18)
n_q_zmassoverr_W1_50, bins_q_zmassoverr_W1_50, patches = plt.hist(q_zmassoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_zmassoverr_W2_50, bins_q_zmassoverr_W2_50, patches = plt.hist(q_zmassoverr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_zmassoverr_W3_50, bins_q_zmassoverr_W3_50, patches = plt.hist(q_zmassoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
n_q_zmassoverr_W4_50, bins_q_zmassoverr_W4_50, patches = plt.hist(q_zmassoverr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_zmassoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_zmassoverr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_zmassoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_zmassoverr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 4])
plt.hist(q_zmassoverr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_zmassoverr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_zmassoverr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
plt.hist(q_zmassoverr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 4], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_zmassoverr_W1_50[np.argmax(n_q_zmassoverr_W1_50)],np.average(q_zmassoverr_W1_50),np.median(q_massoverr_W1_50))
plt.text(1.2, 0.8, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_zmassoverr_W2_50[np.argmax(n_q_zmassoverr_W2_50)],np.average(q_zmassoverr_W2_50),np.median(q_massoverr_W2_50))
plt.text(1.2, 0.6, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_zmassoverr_W3_50[np.argmax(n_q_zmassoverr_W3_50)],np.average(q_zmassoverr_W3_50),np.median(q_massoverr_W3_50))
plt.text(1.2, 0.4, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_zmassoverr_W4_50[np.argmax(n_q_zmassoverr_W4_50)],np.average(q_zmassoverr_W4_50),np.median(q_massoverr_W4_50))
plt.text(1.2, 0.2, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{zM}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6)
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)

plt.subplot(4,5,19)
n_q_zmass2overr_W1_50, bins_q_zmass2overr_W1_50, patches = plt.hist(q_zmass2overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_zmass2overr_W2_50, bins_q_zmass2overr_W2_50, patches = plt.hist(q_zmass2overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_zmass2overr_W3_50, bins_q_zmass2overr_W3_50, patches = plt.hist(q_zmass2overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
n_q_zmass2overr_W4_50, bins_q_zmass2overr_W4_50, patches = plt.hist(q_zmass2overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_zmass2overr_W1_50, histtype='step', color='b', label='W1_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_zmass2overr_W2_50, histtype='step', color='g', label='W2_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_zmass2overr_W3_50, histtype='step', color='r', label='W3_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_zmass2overr_W4_50, histtype='step', color='k', label='W4_50', linewidth=0.5, normed=1, bins=500, range=[0, 7])
plt.hist(q_zmass2overr_W1_75, histtype='step', color='b', label='W1_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_zmass2overr_W2_75, histtype='step', color='g', label='W2_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_zmass2overr_W3_75, histtype='step', color='r', label='W3_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
plt.hist(q_zmass2overr_W4_75, histtype='step', color='k', label='W4_75', linewidth=0.5, normed=1, bins=500, range=[0, 7], linestyle='dotted')
s = "%.3f,%.3f,%.3f" % (bins_q_zmass2overr_W1_50[np.argmax(n_q_zmass2overr_W1_50)],np.average(q_zmass2overr_W1_50),np.median(q_mass3overrrms_W1_50))
plt.text(2, 0.4, s, fontsize=5, color='b')
s = "%.3f,%.3f,%.3f" % (bins_q_zmass2overr_W2_50[np.argmax(n_q_zmass2overr_W2_50)],np.average(q_zmass2overr_W2_50),np.median(q_zmass2overr_W2_50))
plt.text(2, 0.3, s, fontsize=5, color='g')
s = "%.3f,%.3f,%.3f" % (bins_q_zmass2overr_W3_50[np.argmax(n_q_zmass2overr_W3_50)],np.average(q_zmass2overr_W3_50),np.median(q_zmass2overr_W3_50))
plt.text(2, 0.2, s, fontsize=5, color='r')
s = "%.3f,%.3f,%.3f" % (bins_q_zmass2overr_W4_50[np.argmax(n_q_zmass2overr_W4_50)],np.average(q_zmass2overr_W4_50),np.median(q_zmass2overr_W4_50))
plt.text(2, 0.1, s, fontsize=5, color='k')
plt.xlabel(r'$\zeta_\frac{zM^2}{r}$', fontsize=15)
plt.ylabel("Frequency", fontsize=10)
plt.tick_params(axis='x', labelsize=6, direction='up')
plt.tick_params(axis='y', labelsize=6)
plt.setp(plt.xticks()[1], rotation=90)


plt.legend(bbox_to_anchor=(1.5, 4), loc='center left', borderaxespad=0., fontsize=10)

#plt.subplots_adjust(top=0.6)

plt.tight_layout()



plt.savefig('HE0435overdensities.png', dpi=1000)

#plt.show()



print 'Done!'
