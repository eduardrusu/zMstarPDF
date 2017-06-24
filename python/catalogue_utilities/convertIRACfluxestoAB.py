##########################
# Takes outputs from Tphot/Tfit and converts them to Lephare-expected input
##########################


import numpy as np  

file1 = "/Users/eduardrusu/Desktop/WFI2033/WFI2033_IRAC/modifiedforTPHOT/ch1_4amin_nolens_SingleFit_1/ch1_ir_tphot.cat_pass2_best"
file2 = "/Users/eduardrusu/Desktop/WFI2033/WFI2033_IRAC/modifiedforTPHOT/ch2_4amin_nolens_SingleFit_1/ch2_ir_tphot.cat_pass2_best"
file3 = "/Users/eduardrusu/Desktop/WFI2033/WFI2033_IRAC/modifiedforTPHOT/ch3_4amin_nolens_SingleFit_1/ch3_ir_tphot.cat_pass2_best"
file4 = "/Users/eduardrusu/Desktop/WFI2033/WFI2033_IRAC/modifiedforTPHOT/ch4_4amin_nolens_SingleFit_1/ch4_ir_tphot.cat_pass2_best"

data1 = np.loadtxt(file1,usecols=(1,2,7,8),unpack=True)
data2 = np.loadtxt(file2,usecols=(7,8),unpack=True)
data3 = np.loadtxt(file3,usecols=(7,8),unpack=True)
data4 = np.loadtxt(file4,usecols=(7,8),unpack=True)

temperrfunc1 = 0.23
temperrfunc2 = 0.27
temperrfunc3 = 0.30
temperrfunc4 = 0.35

x = data1[0] * 3.0
y = data1[1] * 3.0

ch1 = data1[2]
ch1err = data1[3]
ch1[(ch1 < 0) | (ch1 - ch1err <= 0)] = -99.0
ch1err[ch1 < 0] = -99.0
ch1err[ch1 > 0] = np.sqrt(((21.58 - 2.5 * np.log10(ch1[ch1 > 0]-ch1err[ch1 > 0]) - (21.58 - 2.5 * np.log10(ch1[ch1 > 0]+ch1err[ch1 > 0]))) / 2.0) ** 2 + temperrfunc1 ** 2)
ch1[ch1 > 0] = 21.58 - 2.5 * np.log10(ch1[ch1 > 0])
ch1[ch1err >= 1] = -99.0 # 1 sigma detection
ch1err[ch1err >= 1] = -99.0
#ch1err[ch1 > 0] = 21.58 - 2.5 * np.log10(ch1[ch1 > 0]-ch1err[ch1 > 0]) - (21.58 - 2.5 * np.log10(ch1[ch1 > 0]+ch1err[ch1 > 0])) / 2.0
#ch1err[ch1 > 0] = 21.58 - 2.5 * np.log10(ch1[ch1 > 0]-ch1err[ch1 > 0])
#print ch1[ch1 > 0]
#print ch1err[ch1 > 0]

ch2 = data2[0]
ch2err = data2[1]
ch2[(ch2 < 0) | (ch2 - ch2err <= 0)] = -99.0
ch2err[ch2 < 0] = -99.0
ch2err[ch2 > 0] = np.sqrt(((21.58 - 2.5 * np.log10(ch2[ch2 > 0]-ch2err[ch2 > 0]) - (21.58 - 2.5 * np.log10(ch2[ch2 > 0]+ch2err[ch2 > 0]))) / 2.0) ** 2 + temperrfunc2 ** 2)
ch2[ch2 > 0] = 21.58 - 2.5 * np.log10(ch2[ch2 > 0])
ch2[ch2err >= 1] = -99.0
ch2err[ch2err >= 1] = -99.0

ch3 = data3[0]
ch3err = data3[1]
ch3[(ch3 < 0) | (ch3 - ch3err <= 0)] = -99.0
ch3err[ch3 < 0] = -99.0
ch3err[ch3 > 0] = np.sqrt(((21.58 - 2.5 * np.log10(ch3[ch3 > 0]-ch3err[ch3 > 0]) - (21.58 - 2.5 * np.log10(ch3[ch3 > 0]+ch3err[ch3 > 0]))) / 2.0) ** 2 + temperrfunc3 ** 2)
ch3[ch3 > 0] = 21.58 - 2.5 * np.log10(ch3[ch3 > 0])
ch3[ch3err >= 1] = -99.0
ch3err[ch3err >= 1] = -99.0

ch4 = data4[0]
ch4err = data4[1]
ch4[(ch4 < 0) | (ch4 - ch4err <= 0)] = -99.0
ch4err[ch4 < 0] = -99.0
ch4err[ch4 > 0] = np.sqrt(((21.58 - 2.5 * np.log10(ch4[ch4 > 0]-ch4err[ch4 > 0]) - (21.58 - 2.5 * np.log10(ch4[ch4 > 0]+ch4err[ch4 > 0]))) / 2.0) ** 2 + temperrfunc4 ** 2)
ch4[ch4 > 0] = 21.58 - 2.5 * np.log10(ch4[ch4 > 0])
ch4[ch4err >= 1] = -99.0
ch4err[ch4err >= 1] = -99.0


fileout = "IRAC.cat"
str = "x y ch1 ch1err ch2 ch2err ch3 ch3err ch4 ch4err"
dataout = np.c_[x,y,ch1,ch1err,ch2,ch2err,ch3,ch3err,ch4,ch4err]
np.savetxt(fileout,dataout,header=str,fmt='%.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f')

