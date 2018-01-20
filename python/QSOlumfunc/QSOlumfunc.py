# The code uses the QSO luminosity function from Croom et al. 2009 (SDSS + 2dF), as referenced in Oguri et al. 2012 (AJ,143:120) to compute a redshift probability for a quasar of given observed magnitude. I made sure that those papers use AB. However, I am not completely satisfied because the code also needs to use mag2mag.py with the AGN template from ............to convert from observed magnitude: M_i(z=0)=m_i-DM(z)-K'(z)
# The luminosity function is given in terms of rest-frame g-band absolute magnitudes at z = 2 (i.e., M_g = M_g(z = 2)).

import numpy as np
import os

# observed, magnification corrected QSO magnitudes:
qso_g = 18.92
qso_r = 18.80
qso_i = 18.66
qso_z = 18.53
qso_y = 18.24

Kcorr = np.loadtxt("datafile4.txt",unpack=True) # K-correction from Richards et al. 2006
zlist = Kcorr[0] # the list of redshifts z = np.linspace(0,5.49,550)

h=0.7 # Hubble
alpha =-0.5
beta_h=3.33
beta_1=1.42
Phi_star=1.45*(10**-6)*((h/0.7)**3) # Mpc^-3 mag-1, h=0.7
M_g_star_zero=-22.18 + 5*np.log10(h/0.7)
k_1=1.44
k_2=-0.315
M_g_star = lambda z : M_g_star_zero-2.5*(k_1*z+k_2*(z**2))
Phi = lambda M_g,z: Phi_star/(10**(0.4*(1-beta_h)*(M_g-M_g_star(z)))+10**(0.4*(1-beta_1)*(M_g-M_g_star(z))))

import distances
d = distances.Distance()
DM = lambda z: d.distance_modulus(z)

Phi_list = np.zeros(len(zlist))
for z in range(len(zlist)):
    if z != 0: # error in distance_modulus
        M_i_z0 = qso_i - DM(zlist[z]) + 2.5*(1+alpha)*np.log10(1+zlist[z]) # absolute magnitude at redshift zero assuming a quasar continuum (formula 1 in Richards et al. 2006, and m=M+DM+K)
        M_i_z = M_i_z0 - 2.5*(1+alpha)*np.log10(1+zlist[z])
        M_g = M_i_z + 0.255 - Kcorr[1][z] # for each redshift, find M_g(z=2) which corresponds to m_i and apply the K-correction computed for M_i(z=2)
        #mag2mag = float(os.popen("mag2mag.py -T agn -m1 %s -f1 i_ps1 -z1 %s -f2 i_SDSS -z2 0.00" %(qso_i,zlist[z])).read())
        #M_g = mag2mag + 2.5*(1+alpha)*np.log10(1+zlist[z]) + 0.255 - Kcorr[1][z]
        Phi_list[z] = Phi(M_g,zlist[z])

#Phi_list = np.zeros(len(zlist))
#for z in range(len(zlist)):
#    Phi_list[z] = Phi(-28,zlist[z])

np.savetxt("QLF(z)_iobserved.cat",np.c_[zlist,Phi_list])

import pylab as plt
plt.clf()
plt.xlim([0.4,3])
Phi_list[40:300] = Phi_list[40:300]/np.sum(Phi_list[40:300])
plt.plot(zlist[40:300],Phi_list[40:300])
plt.show()



#import os
#"mag2mag.py -T agn -m1 18.66 -f1 i_ps1 -z1 2.00 -f2 g_SDSS -z2 0.00"
#mag2mag=float(os.popen("mag2mag.py -T agn -m1 18.66 -f1 i_ps1 -z1 2.00 -f2 g_SDSS -z2 0.00").read()) # redirect standard output to variable
# Phi(M_g,z) = Phi_star/(10**(0.4*(1-beta_h)*(M_g-M_g_star(z)))+10**(0.4*(1-beta_1)*(M_g-M_g_star(z))))
# M_g_star(z) = M_g_star_zero-2.5*(k_1*z+k_2*(z**2))
