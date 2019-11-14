# Testing Figure 4 in Nierenberg et al. 2019

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import scipy as sp
#from scipy.stats import norm
import numpy as np

plt.clf()
x = np.random.normal(0, 1, 1000)
chi2 = ((0-x)**2)/(1**2)
argchi2sort = np.argsort(chi2)
f = (argchi2sort + 1)/1000.0
chi2sort = chi2[argchi2sort]
plt.xscale('log')
plt.xlim([0.01,10])
plt.scatter(chi2sort,np.sort(f))
plt.show()
