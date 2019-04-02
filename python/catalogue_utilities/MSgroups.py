
import numpy as np
import os
from stellarpop import distances

mingroup = 5 # I would normally use minimum number of members of any group spectroscopically identified by Dominique in the real data, but that is problematic due tot different sigma_8
mingrouplimmag = 5
magmingroup = 25.5
h = 0.73 # MS
omegaL = 0.75
omega0 = 0.25
zs = 1.662
zl = 0.661
thetaL = 0.95 # COSMOGRAIL VII
dist = distances.Distance()
dist.OMEGA_M = omega0
dist.OMEGA_L = omegaL
dist.h = h
def fbeta(z,zs,zl,dist):
    beta = (dist.angular_diameter_distance(zl, z) * dist.angular_diameter_distance(zs)) / (dist.angular_diameter_distance(z) * dist.angular_diameter_distance(zl, zs))
    if z > zl: return (1-beta)**2
    else: return 1
# indices
halo = 0
specz = 1
posx = 2
posy = 3
mhalo = 4
imag = 5

os.system('rm -f /Users/cerusu/Dropbox/Davis_work/code/WFI2033/8_0_0groups.cat')
f = open("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/8_0_0groups.cat", "a")
root = "/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/original/"
lst = [x for x in os.listdir(root) if (('GGL_los_8_0_0' in x) and ('_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt' in x))]
#print lst[0]
nr = 0
for j in range(len(lst)):
    print j+1,'/',16
    cat = np.loadtxt(root+lst[j],usecols=[1,5,6,7,9,15],comments='GalID')
    catobs = cat[cat[:,imag] <= magmingroup]
    haloID = np.unique(catobs[:,halo],return_counts=True)
    groupID = haloID[0][np.where(haloID[1] >= mingrouplimmag)]
    if j == 0: f.write("# nr galcount posx posy logMhalo specz veldisp thetaE radius \n")
    for i in range(len(groupID)):
        if len(cat[:,halo][cat[:,halo] == groupID[i]]) > mingroup:
            nr += 1
            mhalogroup = np.max(cat[:,mhalo][cat[:,halo] == groupID[i]])
            speczhalo = np.mean(cat[:,specz][cat[:,halo] == groupID[i]])
            sigma = (np.sqrt(omegaL + omega0 * (1+speczhalo)**3) * (mhalogroup/h) / 1200000) ** (1.0/3) # eq. 10 from Finn et al. 2005
            posxgroup = np.mean(cat[:,posx][cat[:,halo] == groupID[i]])
            posygroup = np.mean(cat[:,posy][cat[:,halo] == groupID[i]])
            Dgs = dist.angular_diameter_distance(speczhalo, zs)
            Ds = dist.angular_diameter_distance(zs)
            thetaE = Dgs/Ds * 30 * (sigma/1000) ** 2 # https://ned.ipac.caltech.edu/level5/Mellier/Mellier2_3.html in arcsec
            radius = (((thetaE * thetaL) ** 2) * fbeta(speczhalo,zs,zl,dist)/0.0001) ** (1.0/3) # arcsec
            str = '%d %d %.8f %.8f %.2f %.4f %.2f %.2f %.2f \n' %(nr, len(cat[:,halo][cat[:,halo] == groupID[i]]), posxgroup, posygroup, np.log10(mhalogroup), speczhalo, sigma, thetaE, radius)
            #print nr, len(cat[:,halo][cat[:,halo] == groupID[i]]), posxgroup, posygroup, np.log10(mhalogroup), speczhalo, sigma, thetaE, radius
            f.write(str)
f.close()
