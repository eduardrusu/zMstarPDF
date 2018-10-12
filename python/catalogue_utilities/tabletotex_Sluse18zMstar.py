# Takes a table and returns a custom latex version content
import numpy as np
import os
from os import system
from astropy import units as u
from astropy.coordinates import SkyCoord

################## read catalogue
filein = '/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat'
ra = 2
dec = 3
itot = 4
spec = 40
z_eazy = 48
z_inf = 50
z_sup = 51
mass_best = 92
mass_inf = 93
mass_med = 94
mass_sup = 95
class_bpz = 98
data = np.loadtxt(filein,usecols=[ra,dec,itot,spec,z_eazy,z_inf,z_sup,mass_best,mass_inf,mass_med,mass_sup,class_bpz],unpack=False)
print np.shape(data)
ra = 0
dec = 1
itot = 2
spec = 3
z_eazy = 4
z_inf = 5
z_sup = 6
mass_best = 7
mass_inf = 8
mass_med = 9
mass_sup = 10
class_bpz = 11
lens = SkyCoord(308.4253, -47.39528, unit='deg')
x = SkyCoord(243.079849, 53.012886, unit='deg')
x = SkyCoord(data[:,ra], data[:,dec], unit='deg')
sep = x.separation(lens).arcsec
data = np.c_[data,sep]
sep = 12

################## impose conditions
# for missing mass_med take its error bars to be typical for that mag range
errbar18 = np.median(data[:,mass_med][(data[:,itot]>15) & (data[:,itot]<18) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>15) & (data[:,itot]<18) & (data[:,mass_med]!=-99)])
errbar19 = np.median(data[:,mass_med][(data[:,itot]>18) & (data[:,itot]<19) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>18) & (data[:,itot]<19) & (data[:,mass_med]!=-99)])
errbar20 = np.median(data[:,mass_med][(data[:,itot]>19) & (data[:,itot]<20) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>19) & (data[:,itot]<20) & (data[:,mass_med]!=-99)])
errbar21 = np.median(data[:,mass_med][(data[:,itot]>20) & (data[:,itot]<21) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>20) & (data[:,itot]<21) & (data[:,mass_med]!=-99)])
errbar22 = np.median(data[:,mass_med][(data[:,itot]>21) & (data[:,itot]<22) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>21) & (data[:,itot]<22) & (data[:,mass_med]!=-99)])
errbar23 = np.median(data[:,mass_med][(data[:,itot]>21) & (data[:,mass_med]!=-99)] - data[:,mass_inf][(data[:,itot]>21) & (data[:,mass_med]!=-99)])
data[:,mass_best][data[:,mass_best] < 0] = 9
data[:,mass_best][data[:,mass_med] > 0] = data[:,mass_med][data[:,mass_med] > 0]
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>15) & (data[:,itot]<18)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>15) & (data[:,itot]<18)] - errbar18
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>18) & (data[:,itot]<19)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>18) & (data[:,itot]<19)] - errbar19
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>19) & (data[:,itot]<20)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>19) & (data[:,itot]<20)] - errbar20
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>20) & (data[:,itot]<21)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>20) & (data[:,itot]<21)] - errbar21
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>21) & (data[:,itot]<22)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>21) & (data[:,itot]<22)] - errbar22
data[:,mass_inf][(data[:,mass_inf] < 0) & (data[:,itot]>21)] = data[:,mass_best][(data[:,mass_inf] < 0) & (data[:,itot]>21)] - errbar23
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>15) & (data[:,itot]<18)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>15) & (data[:,itot]<18)] + errbar18
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>18) & (data[:,itot]<19)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>18) & (data[:,itot]<19)] + errbar19
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>19) & (data[:,itot]<20)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>19) & (data[:,itot]<20)] + errbar20
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>20) & (data[:,itot]<21)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>20) & (data[:,itot]<21)] + errbar21
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>21) & (data[:,itot]<22)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>21) & (data[:,itot]<22)] + errbar22
data[:,mass_sup][(data[:,mass_sup] < 0) & (data[:,itot]>21)] = data[:,mass_best][(data[:,mass_sup] < 0) & (data[:,itot]>21)] + errbar23
data = data[data[:,itot] <= 23]
data = data[data[:,class_bpz] >= 0]
data = data[data[:,sep] <= 120]
data[:,z_eazy][data[:,spec] > 0] = data[:,spec][data[:,spec] > 0]
data[:,z_inf][data[:,spec] > 0] = data[:,spec][data[:,spec] > 0]
data[:,z_sup][data[:,spec] > 0] = data[:,spec][data[:,spec] > 0]

print np.shape(data)

################## write tex file
fileout = '/Users/cerusu/GITHUB/H0LiCOW/papers/WFI2033Environment/table_zMstar.tex'
f = open(fileout,'w')
for i in range(np.shape(data)[0]/2):
    f.write('%.5f & $%.5f$ & %.1f & $%.3f_{-%.3f}^{+%.3f}$ & $%.2f_{-%.2f}^{+%.2f}$ & %.5f & $%.5f$ & %.1f & $%.3f_{-%.3f}^{+%.3f}$ & $%.2f_{-%.2f}^{+%.2f}$\\\\\n' % (data[:,ra][2*i],data[:,dec][2*i],data[:,sep][2*i],data[:,z_eazy][2*i],data[:,z_eazy][2*i]-data[:,z_inf][2*i],data[:,z_sup][2*i]-data[:,z_eazy][2*i],data[:,mass_best][2*i],data[:,mass_best][2*i]-data[:,mass_inf][2*i],data[:,mass_sup][2*i]-data[:,mass_best][2*i], data[:,ra][2*i+1],data[:,dec][2*i+1],data[:,sep][2*i+1],data[:,z_eazy][2*i+1],data[:,z_eazy][2*i+1]-data[:,z_inf][2*i+1],data[:,z_sup][2*i+1]-data[:,z_eazy][2*i+1],data[:,mass_best][2*i+1],data[:,mass_best][2*i+1]-data[:,mass_inf][2*i+1],data[:,mass_sup][2*i+1]-data[:,mass_best][2*i+1]))
if np.shape(data)[0] % 2 == 1:
    f.write('%.5f & $%.5f$ & %.1f & $%.3f_{-%.3f}^{+%.3f}$ & $%.2f_{-%.2f}^{+%.2f}$ \\\\\n' % (data[:,ra][2*i+2],data[:,dec][2*i+2],data[:,sep][2*i+2],data[:,z_eazy][2*i+2],data[:,z_eazy][2*i+2]-data[:,z_inf][2*i+2],data[:,z_sup][2*i+2]-data[:,z_eazy][2*i+2],data[:,mass_best][2*i+2],data[:,mass_best][2*i+2]-data[:,mass_inf][2*i+2],data[:,mass_sup][2*i+2]-data[:,mass_best][2*i+2]))
f.close()
