# Run this code in order to sample possible galaxies part of an incomplete spectroscopic group
# no uncertainties assumed on group centroid because they are small enough to ignore

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
#import scipy
#from scipy import stats
#import sys
#import os
#from os import system
#import time
#zgroup = 0.6592
zgroup = 0.4956
limmag = 22.5
faintmagspec = 22.5
center_lens = SkyCoord('20:33:42.080 -47:23:43.00', frame='fk5', unit=(u.hourangle, u.deg))
#center_group = SkyCoord('308.434659 -47.384838', frame='fk5', unit=(u.deg, u.deg)) #  lens redshift
center_group = SkyCoord('308.428143 -47.383722', frame='fk5', unit=(u.deg, u.deg))
sep_group = center_lens.separation(center_group).arcsec
virrad = 240 # arcsec
#observed_members = 19
observed_members = 10
file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat"
filephotozpool = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/grouppool.cat"
data = np.loadtxt(file,usecols=[2,3,4,8,28,29,30,40,97])
ra = 0
dec = 1
i = 2
id = 3
z = 4
zinf = 5
zsup = 6
spec = 7
cls = 8
coord = SkyCoord(ra=data[:,ra]*u.degree, dec=data[:,dec]*u.degree, frame='fk5')
sep = coord.separation(center_lens).arcsec
print 'all',len(data[:,id][(sep <= 120) & (data[:,i] <= limmag)])
print 'specs',len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,spec] > 0)])
#print data[:,i][(sep <= 120) & (data[:,i] >= limmag) & (data[:,spec] > 0)]
#print data[:,spec][(sep <= 120) & (data[:,i] >= limmag) & (data[:,spec] > 0)]
#print data[:,id][(sep <= 120) & (data[:,i] >= limmag) & (data[:,spec] > 0)]
observed120_membersID = data[:,id][(data[:,spec] <= zgroup + 0.01) & (data[:,spec] >= zgroup - 0.01) & (sep <= 120) & (data[:,i] <= limmag)]
print 'members',len(observed120_membersID)
pool = data[:,id][(data[:,spec] == -1) & (sep <= 120) & (data[:,z] - 1.5 * (data[:,z] - data[:,zinf]) <= zgroup) & (data[:,z] + 1.5 * (data[:,zsup] - data[:,z]) >= zgroup) & (data[:,i] <= faintmagspec)]
print pool

# sample uniformly a cube of size 2 x unit radius
x = np.random.uniform(-1,1,10000)
y = np.random.uniform(-1,1,10000)
z = np.random.uniform(-1,1,10000)
cx = 1.0 * sep_group / virrad
cy = 0
rad = 120.0 / virrad
frac = 1.0*len(x[(x**2 + y**2 + z**2 <= 1) & ((x-cx)**2 + (y-cy)**2 <=rad**2)])/len(x[x**2 + y**2 + z**2 <= 1])
# fraction of volume out of the virial sphere which contains the 120"-radius cylinder centered on the lens
def expected_members():
    while True:
        x = np.abs(np.random.normal(50, 17, 1)).astype(int)[0]
        #  fig 16 Berlind et al. 2006 based on velocity dispersion
        if (x > observed_members) and (frac * x > len(observed120_membersID)): break
    return x
for i in range(10):
    missing120_membersID = np.random.choice(a=pool, size=int(frac * expected_members() - len(observed120_membersID)), replace=False)
    missing120_membersra = np.array([])
    missing120_membersdec = np.array([])
    for j in range(len(missing120_membersID)):
        missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==missing120_membersID[j]][0])
        missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==missing120_membersID[j]][0])
    for k in range(len(observed120_membersID)):
        missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==observed120_membersID[k]][0])
        missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==observed120_membersID[k]][0])
    #np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgrouphandpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
    np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgroup049handpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
