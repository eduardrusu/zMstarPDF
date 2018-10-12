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
itot_err = 5
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
ch1 = 74
ch1_err = 75
ch2 = 76
ch2_err = 77
ch3 = 78
ch3_err = 79
ch4 = 80
ch4_err = 81
class_eazy = 98
data = np.loadtxt(filein,usecols=[ra,dec,itot,itot_err,u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err,ch1,ch1_err,ch2,ch2_err,ch3,ch3_err,ch4,ch4_err,class_eazy],unpack=False)
print np.shape(data)
ra = 0
dec = 1
itot = 2
itot_err = 3
u = 4
u_err = 5
g = 6
g_err = 7
r = 8
r_err = 9
i = 10
i_err = 11
z = 12
z_err = 13
Y = 14
Y_err = 15
J = 16
J_err = 17
H = 18
H_err = 19
K = 20
K_err = 21
ch1 = 22
ch1_err = 23
ch2 = 24
ch2_err = 25
ch3 = 26
ch3_err = 27
ch4 = 28
ch4_err = 29
class_eazy = 30
lens = SkyCoord(308.4253, -47.39528, unit='deg')
x = SkyCoord(243.079849, 53.012886, unit='deg')
x = SkyCoord(data[:,ra], data[:,dec], unit='deg')
sep = x.separation(lens).arcsec
data = np.c_[data,sep]
sep = 31

################## impose conditions
data = data[data[:,itot] <= 23]
data = data[data[:,class_eazy] >= 0]
data = data[data[:,sep] <= 120]
print np.shape(data)
fileout = '/Users/cerusu/GITHUB/H0LiCOW/papers/WFI2033Environment/table_phot.tex'
np.savetxt(fileout,np.c_[data[:,ra],data[:,dec],data[:,itot],data[:,itot_err],data[:,u],data[:,u_err],data[:,g],data[:,g_err],data[:,r],data[:,r_err],data[:,i],data[:,i_err],data[:,z],data[:,z_err],data[:,Y],data[:,Y_err],data[:,J],data[:,J_err],data[:,H],data[:,H_err],data[:,K],data[:,K_err],data[:,ch1],data[:,ch1_err],data[:,ch2],data[:,ch2_err],data[:,ch3],data[:,ch3_err],data[:,ch4],data[:,ch4_err]],fmt='%.5f & $%.5f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ \\\\')
strin = '-99.00 \\\pm -99.00'
strout = '-'
os.system("sed -i -e \'s/%s/%s/g\' %s" % (strin,strout,fileout))

#np.savetxt(masterfile[:-4] + "_WFI2033noIRACeazy_nobeta_testduplicate.cat",data.T,fmt='%s %s %s %s %.2f %.2f %d %d %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f %.2f %.4f %.4f')
