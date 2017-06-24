# Creates figure 4 from Rusu et al. 2017

import numpy as np
import sys
import os
from os import system
import time
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection

plt.clf()

def plot(ind,px,py,str,mult):
    ax=plt.subplot(4,4,ind, sharex=ax1, sharey=ax1)
    plt.axis([-130, 130, -130, 130])
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='both', which='major', labelsize=6)
    if (ind == 10) | (ind == 11) | (ind == 12) | (ind == 13) | (ind == 14):
        plt.xticks(rotation='vertical')
    ax.text(px, py, str, fontsize=fontlabel, color='k',transform=ax.transAxes)
    xx = data[:,ind+2][data[:,2]>maglim2].astype(float)*mult # multiplier to make the bubbles have an appropriate radius
    xx[xx<2]=2 # minimum visible bubble size
    yy = data[:,ind+2][data[:,2]<=maglim2].astype(float)*mult
    yy[yy<2]=2
    plt.scatter((-pix/2 + data[:,0][data[:,2]>maglim2])*scale, (-pix/2 + data[:,1][data[:,2]>maglim2])*scale, s=xx.astype(int), c='black', lw = 0, alpha=0.5) # typically red
    plt.scatter((-pix/2 + data[:,0][data[:,2]<=maglim2])*scale, (-pix/2 + data[:,1][data[:,2]<=maglim2])*scale, s=yy.astype(int), c='black', lw = 0, alpha=0.5) # typically blue

fontabsciss = 9
fontlabel = 7
nRows = 4
nCols = 4
nPlots = 14

z_s = 1.66
pix = 915 # number of pixels in 4 arcmin
scale = 240.0 / pix
maglim1 = 23 # for plots when I want blue and red points I would use 24
maglim2 = 23

x,y,i,classify,z,mstar,mhalo = np.loadtxt("/Users/eduardrusu/Dropbox/Davis_work/code/WFI2033/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W_WFI2033IRACbpz_nobeta.cat",usecols=[0,1,4,7,8,9,10],unpack=True)
sep = np.sqrt((x - 457.5)**2 + (y - 457.5)**2) * scale

x = x[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
y = y[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
z = z[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
i_ = i[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
sep_ = sep[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
mstar = mstar[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
mhalo = mhalo[(sep <= 120) & (i <= maglim1) & (classify >= 0)]
sep = sep_
i = i_

wht_1 = np.ones(len(x))
wht_z = z_s * z - z * z
wht_M = 10 ** mstar
wht_M2 = (10 ** mstar) ** 2
wht_M3 = (10 ** mstar) ** 3
wht_1r = 1/sep
wht_zr = wht_z * wht_1r
wht_Mr = wht_M * wht_1r
wht_M2r = wht_M2 * wht_1r
wht_M3r = wht_M3 * wht_1r
wht_flexion = wht_M * (wht_1r ** 3)
wht_tidal = wht_M * (wht_1r ** 2)
wht_SIS = np.sqrt(wht_M) * wht_1r
wht_SIShalo = np.sqrt(10 ** mhalo) * wht_1r
data = np.c_[x,y,i,wht_1,wht_z,wht_M,wht_M2,wht_M3,wht_1r,wht_zr,wht_Mr,wht_M2r,wht_M3r,wht_flexion,wht_tidal,wht_SIS,wht_SIShalo]

# 0 x
# 1 y
# 2 i
# 3 wht_1
# 4 wht_z
# 5 wht_M
# 6 wht_M2
# 7 wht_M3
# 8 wht_1r
# 9 wht_zr
# 10 wht_Mr
# 11 wht_M2r
# 12 wht_M3r
# 13 wht_flexion
# 14 wht_tidal
# 15 wht_SIS
# 16 wht_SIShalo

plt.clf()

plt.axis([-130, 130, -130, 130])

fig = plt.figure()
ax1 = fig.add_subplot(4,4,1)
ax = plt.subplot(4,4,1, sharex=ax1, sharey=ax1)
ax.set_aspect(1, adjustable='datalim')

for i in range(14):
    
    if i == 0: plot(i+1,0.9,0.85,"$1$",5)
    if i == 1: plot(i+1,0.9,0.85,"$z$",10)
    if i == 2: plot(i+1,0.9,0.85,"$M_\star$",1.0/1000000000)
    if i == 3: plot(i+1,0.85,0.85,"$M^2_\star$",1.0/100000000000000000000)
    if i == 4: plot(i+1,0.85,0.85,"$M^3_\star$",1.0/20000000000000000000000000000000)
    if i == 5: plot(i+1,0.8,0.85,"$1/r$",200)
    if i == 6: plot(i+1,0.8,0.85,"$z/r$",400)
    if i == 7: plot(i+1,0.8,0.85,"$M_\star/r$",1.0/100000000)
    if i == 8: plot(i+1,0.8,0.85,"$M^2_\star/r$",1.0/10000000000000000000)
    if i == 9: plot(i+1,0.8,0.85,"$M^3_\star/r$",1.0/1000000000000000000000000000000)
    if i == 10: plot(i+1,0.8,0.85,"$M_\star/r^3$",1.0/2000000)
    if i == 11: plot(i+1,0.8,0.85,"$M_\star/r^2$",1.0/10000000)
    if i == 12: plot(i+1,0.75,0.85,"$\sqrt{M_\star}/r$",1.0/500)
    if i == 13: plot(i+1,0.7,0.85,"$\sqrt{M_h}/r$",1.0/10000)

    circle120=plt.Circle((0,0),120,color='k',fill=False)
    circle45=plt.Circle((0,0),45,color='k',fill=False)
    fig = plt.gcf()
    fig.gca().add_artist(circle120)
    fig.gca().add_artist(circle45)

# hide the plots with no data in the grid
ax=plt.subplot(4,4,15, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
ax.set_frame_on(False)
ax=plt.subplot(4,4,16, sharex=ax1, sharey=ax1)
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
ax.set_frame_on(False)

index = 1
for r in range(1, nRows +1):
    for c in range(1, nCols + 1):
        ax = plt.subplot(nRows, nCols, index, sharex=ax1, sharey=ax1)
        index += 1
         # Turn off y tick labels for all but the first column.
        if ((c != 1) and (index <= nPlots)):
             plt.setp(ax.get_yticklabels(), visible=False)
        #if c == 1:
            #plt.ylabel('arcsec',fontsize=fontabsciss)
          # Turn off x tick lables for all but the bottom plot in each 
          # column. 
        if ((nPlots - index) >= nCols):
             plt.setp(ax.get_xticklabels(), visible=False)
             #plt.set_aspect('equal')
        #if (index == 15) or (index == 16) or (index == 17) or (index == 18):
             #plt.xlabel('arcsec',fontsize=fontabsciss)
        if index == 15:
            plt.setp(ax.get_yticklabels(), visible=False)

fig.text(0.5, 0.04, 'radius [arcsec]', ha='center', va='center')
fig.text(0.06, 0.5, 'radius [arcsec]', ha='center', va='center', rotation='vertical')
plt.subplots_adjust(wspace=0, hspace=0)

#plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.6)

plt.savefig('bubbles.png', dpi=500, bbox_inches='tight')


print 'Done!'
