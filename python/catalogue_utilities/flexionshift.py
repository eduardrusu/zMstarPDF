# calculates flexion shift given a catalogue of stellar masses, redshifts and positions on the sky. Check Behroozi et al 2010, Auger et al 2009 and Jabran Zahid et al 2017 for the applicability of the formulas.
# To produce the input, run awk '{if (($5<=23) && ($98>=0)) print $1,$2,$3,$4,$5,$6,$29,$30,$31,$41,$86,$87,$88,$89}' rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat > WFI2033gal_KenDominique.cat

import numpy as np
import scipy
from scipy import stats
import sys
import os
from os import system
from astropy.coordinates import SkyCoord
from astropy import units as u
import distances

####################
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
d0 = 0.55
da = -0.12
g0 = 1.12
ga = -0.53

####################
# Auger et al 2009 parameters:
a_aug = 0.23
b_aug = 0.14

####################
# Jabran Zahid et al 2017 parameters:
sb = 10**(2.07)
Mb = 10**(10.26)
a1 = 0.403
a2 = 0.293

####################
lens = SkyCoord('20:33:42.16 -47:23:44.20', frame='fk5', unit=(u.hourangle, u.deg)) # center of the lensing galaxy
z_s = 1.66
z_d = 0.66
theta_lens = 0.96 # arcsec # Einstein radius from the Vuissoz et al. 2008 SIE+gamma model
dist = distances.Distance()
dist.OMEGA_M = 0.31 # Planck
dist.OMEGA_L = 0.69 # Planck
dist.h = 0.68 # Planck
D_S = dist.angular_diameter_distance(0,z_s) # in Mpc
D_DS = dist.angular_diameter_distance(z_d,z_s)



const = 9.626*(10**(-20)) # 4GM_sol/c^2 in Mpc
radinsec = 206265  # radian in arcsec


# read from file
ra_ = 2
dec_ = 3
z_b_ = 6
spec_ = 9
mstarbest_ = 10
mstarmed_ = 12

file = "WFI2033gal_KenDominique.cat"
ra,dec,z_b,spec,mstarbest,mstarmed,sep,flex_halo,flex_sluse,flex_zahid = np.loadtxt(file,usecols=[ra_,dec_,z_b_,spec_,mstarbest_,mstarmed_,ra_,ra_,ra_,ra_],unpack=True)
mstarmed[mstarmed < 0] = mstarbest[mstarmed < 0]
mstar = mstarmed
mstar_sluse = mstar - np.log10(0.55) # convert to Salpeter, but only for Auger et al 2009 (Sluse), becasue I computed the masses assuming Chabrier, and so did Zahid et al.
sigma_zahid = np.zeros(len(mstar))
z_b[spec>0] = spec[spec>0]

for i in range(len(ra)):
    coord=SkyCoord(ra=ra[i]*u.degree, dec=dec[i]*u.degree, frame='fk5')
    sep_lens = coord.separation(lens).arcsec
    sep[i] = sep_lens
    z = z_b[i]
    a = 1 / (1 + z)
    if z <= 1:
        logM1a = M10_ + M1a_ * (a - 1)
        logMs0a = Ms00_ + Ms0a_ * (a - 1)
        notlogMs0a = 10 ** logMs0a
        b = b0_ + ba_ * (a-1)
        d = d0_ + da_ * (a-1)
        g = g0_ + ga_ * (a-1)
        mhalo = logM1a + b * (mstar[i] - logMs0a) + ((10 ** mstar[i]/notlogMs0a)**d)/(1+(10 ** mstar[i]/notlogMs0a)**(-g)) - 1/2
    else:
        logM1a = M10 + M1a * (a - 1)
        logMs0a = Ms00 + Ms0a * (a - 1)
        notlogMs0a = 10 ** logMs0a
        b = b0 + ba * (a-1)
        d = d0 + da * (a-1)
        g = g0 + ga * (a-1)
        mhalo = logM1a + b * (mstar[i] - logMs0a) + ((10 ** mstar[i]/notlogMs0a)**d)/(1+(10 ** mstar[i]/notlogMs0a)**(-g)) - 1/2
    if mstar_sluse[i] > 10.5:
        fDM = a_aug * (mstar_sluse[i] - 11) + b_aug
        #print fDM
        mtotal = (10**mstar_sluse[i]) / (1 - fDM)
    beta = dist.angular_diameter_distance(z_d,z) * D_S / (dist.angular_diameter_distance(0,z) * D_DS)
    if z > z_d:
        fbeta = (1 - beta)**2
        #print fbeta
    else:
        fbeta = 1
    if 10**(mstar[i]) < Mb:
        sigma_zahid[i] = sb * ((10**mstar[i] / Mb)**a1)
    else:
        sigma_zahid[i] = sb * ((10**mstar[i] / Mb)**a2)
    theta_halo = np.sqrt((10 ** mhalo) / (10 ** 11.09)) * np.sqrt(1000 * dist.angular_diameter_distance(z,z_s)/(dist.angular_diameter_distance(0,z) * D_S)) # https://en.wikipedia.org/wiki/Einstein_radius I already checked this formula manually
    theta_zahid = radinsec * 4 * 3.14 * (sigma_zahid[i]/300000)**2 * dist.angular_diameter_distance(z,z_s) / D_S
    #print theta_zahid
    flex_halo[i] = fbeta * (theta_lens * theta_halo)**2 / (sep[i]**3)
    flex_zahid[i] = fbeta * (theta_lens * theta_zahid)**2 / (sep[i]**3)
    if mstar_sluse[i] > 10.5:
        theta_sluse = np.sqrt(mtotal / (10 ** 11.09)) * np.sqrt(1000 * dist.angular_diameter_distance(z,z_s)/(dist.angular_diameter_distance(0,z) * D_S)) # https://en.wikipedia.org/wiki/Einstein_radius
        flex_sluse[i] = fbeta * (theta_lens * theta_sluse)**2 / (sep[i]**3)
    else:
        flex_sluse[i] = flex_zahid[i]
    if z >= z_s:
        flex_halo[i] = 0
        flex_sluse[i] = 0
        flex_zahid[i] = 0
head = "sep     flex_halo   flex_sluse  flex_zahid veldisp"
np.savetxt(file[:-4] + "_flexionshift.cat",np.c_[sep,flex_halo,flex_sluse,flex_zahid,sigma_zahid],fmt='%.2f %.4e %.4e %.4e %.2f',header=head)


