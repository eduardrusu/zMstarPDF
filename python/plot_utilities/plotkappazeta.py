#

import pylab as plt
import numpy as np

file = '/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappadependencezeta/medstd.dat'
data = np.loadtxt(file,usecols=[1])
mag = np.array([22.5,23,23.5,24])
kappa45zeta05 = np.array([data[1],data[19],data[20],data[13]])
kappa45zeta1 = np.array([data[5],data[6],data[7],data[18]])
kappa45zeta15 = np.array([data[3],data[12],data[14],data[21]])
kappa120zeta05 = np.array([data[4],data[9],data[11],data[15]])
kappa120zeta1 = np.array([data[2],data[22],data[17],data[10]])
kappa120zeta15 = np.array([data[0],data[16],data[23],data[8]])
plt.clf()
plt.plot(mag,kappa45zeta05,label='45 $\zeta=0.5$',color='b')
plt.plot(mag,kappa45zeta1,label='45 $\zeta=1.0$',color='k')
plt.plot(mag,kappa45zeta15,label='45 $\zeta=1.5$',color='r')
plt.plot(mag,kappa120zeta05,label='120 $\zeta=0.5$',linestyle='-.',color='b')
plt.plot(mag,kappa120zeta1,label='120 $\zeta=1.0$',linestyle='-.',color='k')
plt.plot(mag,kappa120zeta15,label='120 $\zeta=1.5$',linestyle='-.',color='r')

#plt.plot(kappa_values[:-1][::1],kappa_1[::1], linewidth=2, label ='$120: 1 + 1/r + \gamma$', linestyle='-.')
plt.xlabel(r'mag', fontsize=20)
plt.xlim([22.5,24])
#plt.set_xticklabels([22.5,23,23.5,24])
plt.ylabel(r'$\kappa_\mathrm{med}$', fontsize=20)
plt.legend(loc="lower right")
plt.xticks(np.linspace(22.5,24,4))
#plt.show()
plt.savefig('%s.png' % file[:-4], dpi=250, bbox_inches='tight')
