# Combines samples from multiple Gaussians with asymmetric error bars

import numpy as np
import pylab as plt

samples = 5000

med = np.array([0.316,0.371,0.349,0.332,0.275,0.291,0.360,0.275,0.360,0.349,0.337])
stdsup = np.array([+0.041,+0.033,+0.036,+0.028,+0.052,+0.034,+0.043,+0.031,+0.019,+0.038,+0.047])
stdinf = np.array([-0.007,-0.041,-0.019,-0.012,-0.002,-0.003,-0.003,-0.013,-0.02,-0.026,-0.013])

combined = np.array([])
for i in range(len(med)):
    left = med[i] - abs(med[i] - np.random.normal(med[i],abs(stdinf[i]),samples))
    right = med[i] + abs(med[i] - np.random.normal(med[i],stdsup[i],samples))
    combined = np.r_[combined,left,right]

print np.percentile(combined,[16,50,84])
plt.clf()
plt.hist(combined,bins=100)
plt.show()
