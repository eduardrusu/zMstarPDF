# Plots the equivalent of Figure 13 in Rusu et al. 2017

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import scipy as sp
from scipy.stats import norm
from scipy.optimize import curve_fit
import numpy as np
import os
import glob

# I used inferkappa_unbiasedwithshear45and120FITSio_customzeta.py on J1206 24 with pure number counts, zeta=0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0(,2.6)
# with width +/-0.05 to explore the relation between kappamed and sigma:
#medpdf_45=np.array([ -0.036,-0.025,-0.014,-0.000,0.017,0.036,0.060,0.087,0.180])
#stdpdf_45=np.array([0.018 ,0.021,0.024,0.029,0.035,0.042,0.051,0.061,0.092])
#medpdf_120=np.array([-0.045,-0.030,-0.012,0.009,0.034,0.064,0.097,0.118])
#stdpdf_120=np.array([0.017,0.021,0.026,0.033,0.042,0.052,0.068,0.076])

# I used inferkappa_unbiasedwithshear45and120FITSio_customzeta.py on WFI2033 22.5 with pure number counts, zeta=0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,(2.6,3.0,3.4)
# with width +/-0.1 and E=1 to explore the relation between kappamed and sigma:
medpdf_45=np.array([ -0.027,-0.022,-0.016,-0.010,-0.004,0.003,0.011,0.020,0.031,0.043,0.056,0.071,0.086,0.101,0.118,0.135,0.155,0.181,0.202,0.220])
stdpdf_45=np.array([0.026,  0.026, 0.028, 0.031, 0.033, 0.037,0.042,0.044,0.050,0.061,0.063,0.067,0.075,0.085,0.092,0.105,0.116,0.132,0.137,0.132])
medpdf_120=np.array([-0.034,-0.025,-0.016,-0.005,0.006,0.020,0.036,0.054,0.074,0.094,0.118,0.139,0.159])
stdpdf_120=np.array([0.024, 0.027,  0.032,0.035, 0.040,0.047,0.053,0.063,0.076,0.090,0.094,0.104,0.121])

def func(x, a, b):
    return a * x + b

popt45, pcov45 = curve_fit(func, medpdf_45, stdpdf_45)
popt120, pcov120 = curve_fit(func, medpdf_120, stdpdf_120)
#plt.plot(medpdf_45, stdpdf_45, 'b-', label='data45')
#plt.plot(medpdf_120, stdpdf_120, 'r-', label='data120')
#plt.plot(medpdf_45, func(medpdf_45, *popt45), 'b-.',label='fit: a=%5.3f, b=%5.3f' % tuple(popt45))
#plt.plot(medpdf_120, func(medpdf_120, *popt120), 'r-.',label='fit: a=%5.3f, b=%5.3f' % tuple(popt120))
#plt.legend()
#plt.show()
popt = np.mean([popt45,popt120],axis=0)

# But the way inferkappa_unbiasedwithshear45and120FITSio_customzeta.py computes std is not a simple np.std, which is output by inferkappasimbias.py
# Normalizing std(truekappa-medkappa):
root = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/E2new5000_2/"
#x = np.loadtxt(root+'kappasim_WFI2033_measured_5innermask_nobeta_zgap-1.0_-1.0_45_gal_22.5_med_1.0_overdensities1.0.cat',unpack=True)
x = np.loadtxt(root+'kappasim_WFI2033_measured_5innermask_nobeta_zgap-1.0_-1.0_45_gal_22.5_med_overdensities1.44.cat',unpack=True)
#x = np.loadtxt(root+'kappasim_WFI2033_measured_5innermask_nobeta_zgap-1.0_-1.0_120_gal_22.5_med_overdensities1.55.cat',unpack=True)

list = glob.glob(root+'kappasim_WFI2033*.cat')
fout = root + 'medstdbias.dat'
os.system('rm -f %s' % fout)

for i in range(len(list)):
    file = list[i]
    data  = np.loadtxt(file)
    data = data[data[:,4] >= 3] # at least 3 LOS
    data = data.T
    if np.shape(data)[1] != 0:
        scaledstd = np.std(data[1]-data[0]) * (func(np.median(x[1]), *popt)/func(np.median(data[1]), *popt)) / np.std(x[1]-x[0])
        str = "%s %.3f %.3f %d %d %.3f \n" % (list[i],np.median(data[0]-data[1]),np.std(data[0]-data[1]),np.median(data[4]),len(data[0]),scaledstd)
        #str = "%s %.3f %.3f %d %d \n" % (list[i],np.median(data[1]),np.std(data[0]-data[1]),np.median(data[4]),len(data[0]))
    else: str = "%s 0.000 0.000 0 0 0.000 \n" % (list[i])
    file =open(fout,'a')
    file.write(str)
file.close()
