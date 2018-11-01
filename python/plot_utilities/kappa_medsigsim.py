# Plots the equivalent of Figure 13 in Rusu et al. 2017

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import scipy as sp
from scipy.stats import norm
import numpy as np
import os
import glob

min_kappa = -0.10
max_kappa = 1
bin_stat = 2000
halfwidth = (max_kappa - min_kappa) / (bin_stat * 2.0)

root = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/"
list = glob.glob(root+'kappasim_*.cat')

fout = root + 'medstd.dat'
#file =open('fout','a')
os.system('rm -f %s' % fout)

for i in range(len(list)):
    file = list[i]
    kappa  = np.loadtxt(file, usecols=[0], unpack=True, comments = '#')
    median,stddev,kappa_values = statistics(kappa,bin_stat,min_kappa,max_kappa)
    winlen = 12
    median_smooth,stddev_smooth,kappa_values_smooth = statistics(smooth(kappa,winlen,'flat')[(winlen/2-1):-(winlen/2)],bin_stat,min_kappa,max_kappa)
    str = "%s %.3f %.3f %.3f %.3f \n" % (list[i],median,stddev,median_smooth,stddev_smooth)
    file =open(fout,'a')
    file.write(str)
file.close()
