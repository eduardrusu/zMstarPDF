# the code transforms Mstar to Mhalo using Behroozi et al. 2010. I used this to test if the MS De Lucia, Guo or Bower galaxies can be used to recover Mhalo from Mstar

import numpy as np
import sys

file = "/Users/perseus/Desktop/GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt"
fileout = file[:-4] + "plot.txt"

zspec = 5
mhalo = 9
mstar = 11

data = np.loadtxt(file,usecols=[zspec,mhalo,mstar],unpack=True)

zspec = 0
mhalo = 1
mstar = 2

data[mstar] = np.log10(data[mstar])
data[mhalo] = np.log10(data[mhalo])

# Behroozi et al 2010 parameters for z < 1:
M10_ = 12.35
M1a_ = 0.28
Ms00_ = 10.72
Ms0a_ = 0.55
b0_ = 0.44
ba_ = 0.18
d0_ = 0.57
da_ = 0.17
g0_ = 1.56
ga_ = 2.51
# z >= 1:
M10 = 12.27
M1a = -0.84
Ms00 = 11.09
Ms0a = 0.56
b0 = 0.65
ba = 0.31
d0 = 0.56
da = -0.12
g0 = 1.12
ga = -0.53

datahalo = np.zeros(len(data[zspec]))

a = 1 / (1 + data[zspec][data[zspec] <= 1])
logM1a = M10_ + M1a_ * (a - 1)
logMs0a = Ms00_ + Ms0a_ * (a - 1)
notlogMs0a = 10 ** logMs0a
b = b0_ + ba_ * (a - 1)
d = d0_ + da_ * (a - 1)
g = g0_ + ga_ * (a - 1)
datahalo[data[zspec] <= 1] = logM1a + b * (data[mstar][data[zspec] <= 1] - logMs0a) + ((10 ** data[mstar][data[zspec] <= 1]/notlogMs0a)**d)/(1+(10 ** data[mstar][data[zspec] <= 1]/notlogMs0a)**(-g)) - 1/2
del logM1a
del logMs0a
del notlogMs0a
del b
del d
del g

a = 1 / (1 + data[zspec][data[zspec] > 1])
logM1a = M10 + M1a * (a-1)
logMs0a = Ms00 + Ms0a * (a-1)
notlogMs0a = 10 ** logMs0a
b = b0 + ba * (a-1)
d = d0 + da * (a-1)
g = g0 + ga * (a-1)
datahalo[data[zspec] > 1] = logM1a + b * (data[mstar][data[zspec] > 1] - logMs0a) + ((10 ** data[mstar][data[zspec] > 1]/notlogMs0a)**d)/(1+(10 ** data[mstar][data[zspec] > 1]/notlogMs0a)**(-g)) - 1/2
del logM1a
del logMs0a
del notlogMs0a
del b
del d
del g

dataout = np.c_[data[mhalo],data[mstar],datahalo]
head = "M_Halo \t M_Stellar \t M_Halo_Behroozi"
np.savetxt(fileout,dataout,header=head,fmt='%.3e \t %.3e \t %.3e')

print ' Done!'
