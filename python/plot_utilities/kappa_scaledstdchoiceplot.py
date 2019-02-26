# uses the output of kappa_medsigsim.py to decide the weight-based combination

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/scaledstdchoiceplot.dat"
names = np.genfromtxt(file,usecols=[0],dtype='S200')
data = np.loadtxt(file,usecols=[1,2,3,4,5,6,19,20],dtype='d,d,d,d,d,d,f,f',unpack=True)
medofstd,rmsofstd = data[6],data[7]
N = len(data[0])
los = np.zeros(N)

for i in range(N):
    los[i] = np.log10(np.median([data[0][i],data[1][i],data[2][i],data[3][i],data[4][i],data[5][i]]))
    trim = np.char.find(names[i], '_22.5')
    names[i] = names[i][:trim]

ind = 1 * np.arange(N)  # the x locations for the groups
width = 0.9      # the width of the bars

ax = plt.subplot(2,1,1)
col1 = los
rects1 = ax.bar(ind + width, col1, width)
ax.set_ylim([0,7])
ax.set_xlim([0,72])
ax.set_ylabel('log(No. of LOS)',fontsize=6)
ax.set_xticks(ind + width)
ax.set_xticklabels([])

ax = plt.subplot(2,1,2)
col2 = medofstd
col2err = rmsofstd
rects2 = ax.bar(ind + width, col2, width, yerr=rmsofstd)
ax.set_ylim([0,1.3])
ax.set_xlim([0,72])
ax.set_ylabel('scaled std($\kappa_{med} - \kappa_{true}$)',fontsize=6)
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(names, fontsize=5, rotation='vertical')
#ax.set_aspect(.0)
#w, h = figaspect(2.)
#fig = Figure(figsize=(w, h))
plt.subplots_adjust(left=0.08, bottom=0.5, right=0.99, top=0.95, wspace=0.7, hspace=0)
plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/scaledstdchoiceplot.png', dpi=250)
plt.clf()
