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
medpdf_45=np.array([ -0.036,-0.025,-0.014,-0.000,0.017,0.036,0.060,0.087,0.180])
stdpdf_45=np.array([0.018 ,0.021,0.024,0.029,0.035,0.042,0.051,0.061,0.092])
medpdf_120=np.array([-0.045,-0.030,-0.012,0.009,0.034,0.064,0.097,0.118])
stdpdf_120=np.array([0.017,0.021,0.026,0.033,0.042,0.052,0.068,0.076])

def func(x, a, b):
    return a * x + b

popt45, pcov45 = curve_fit(func, medpdf_45, stdpdf_45)
popt120, pcov120 = curve_fit(func, medpdf_120, stdpdf_120)
#plt.plot(medpdf_45, stdpdf_45, 'b-', label='data45')
#plt.plot(medpdf_120, stdpdf_120, 'b-', label='data120')
#plt.plot(medpdf_45, func(medpdf_45, *popt45), 'r-',label='fit: a=%5.3f, b=%5.3f' % tuple(popt45))
#plt.plot(medpdf_120, func(medpdf_120, *popt120), 'r-',label='fit: a=%5.3f, b=%5.3f' % tuple(popt120))
#plt.legend()
#plt.show()
popt = np.mean(popt45,popt120,axis=0)

# But the way inferkappa_unbiasedwithshear45and120FITSio_customzeta.py computes std is not a simple np.std, which is output by inferkappasimbias.py
# Normalizing std(truekappa-medkappa):
root = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/"
x = np.loadtxt(root+'calibs/kappasim_WFI2033_measured_5innermask_nobeta_zgap-1.0_-1.0_45_gal_22.5_med_1.8_overdensities1.8.cat',unpack=True)

func(0.01, *popt)







list = glob.glob(root+'kappasim_WFI2033*.cat')

#fout = root + 'medstd.dat'
#os.system('rm -f %s' % fout)

for i in range(len(list)):
    file = list[i]
    kappa  = np.loadtxt(file, usecols=[0], unpack=True, comments = '#')
    #median,stddev,kappa_values = statistics(kappa,bin_stat,min_kappa,max_kappa)
    winlen = 12
    #median_smooth,stddev_smooth,kappa_values_smooth = statistics(smooth(kappa,winlen,'flat')[(winlen/2-1):-(winlen/2)],bin_stat,min_kappa,max_kappa)
    #str = "%s %.3f %.3f %.3f %.3f \n" % (list[i],median,stddev,median_smooth,stddev_smooth)
    #file =open(fout,'a')
    #file.write(str)
#file.close()
