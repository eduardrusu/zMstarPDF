# Simple plot with custom labels

import matplotlib.pyplot as plt
import numpy as np

plt.clf()

SIS_b = 1.45
SIS_c = 1.54
SIS_d = 0.52
SIS_bc = 0.95
SIS_bd = 0.36
SIS_cd = 0.38
SIE_b = 1.47
SIE_c = 1.39
SIE_d = 0.58
SIE_bc = 0.95
SIE_bd = 0.39
SIE_cd = 0.42

data = np.loadtxt('fluxratioerr.cat', dtype={'names': ('filter', 'B/A', 'B/Ae', 'C/A', 'C/Ae', 'D/A', 'D/Ae', 'C/B', 'C/Be', 'D/B', 'D/Be', 'D/C', 'D/Ce'),'formats': ('S2', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
filt = [x[0] for x in data]
b = [x[1] for x in data]
b_e = [x[2] for x in data]
c = [x[3] for x in data]
c_e = [x[4] for x in data]
d = [x[5] for x in data]
d_e = [x[6] for x in data]
bc = [x[7] for x in data]
bc_e = [x[8] for x in data]
bd = [x[9] for x in data]
bd_e = [x[10] for x in data]
cd = [x[11] for x in data]
cd_e = [x[12] for x in data]

x=np.linspace(1,len(filt),len(filt))
plt.xticks(x, filt)
plt.plot(x, b, label=data.dtype.names[1], color='b')
plt.errorbar(x, b, yerr=b_e, color='b')
plt.plot(x, x*0+SIS_b, color='b', linestyle='--')
#plt.plot(x, x*0+SIE_b, color='b', linestyle=':')

plt.plot(x, c, label=data.dtype.names[3], color='g')
plt.errorbar(x, c, yerr=c_e, color='g')
plt.plot(x, x*0+SIS_c, color='g', linestyle='--')
#plt.plot(x, x*0+SIE_c, color='g', linestyle=':')

plt.plot(x, d, label=data.dtype.names[5], color='r')
plt.errorbar(x, d, yerr=d_e, color='r')
plt.plot(x, x*0+SIS_d, color='r', linestyle='--')
#plt.plot(x, x*0+SIE_d, color='r', linestyle=':')

plt.plot(x, bc, label=data.dtype.names[7], color='m')
plt.errorbar(x, bc, yerr=bc_e, color='m')
plt.plot(x, x*0+SIS_bc, color='m', linestyle='--')
#plt.plot(x, x*0+SIE_bc, color='m', linestyle=':')

plt.plot(x, bd, label=data.dtype.names[9], color='c')
plt.errorbar(x, bd, yerr=bd_e, color='c')
plt.plot(x, x*0+SIS_bd, color='c', linestyle='--')
#plt.plot(x, x*0+SIE_bd, color='c', linestyle=':')

plt.plot(x, cd, label=data.dtype.names[11], color='y')
plt.errorbar(x, cd, yerr=cd_e, color='y')
plt.plot(x, x*0+SIS_cd, color='y', linestyle='--')
#plt.plot(x, x*0+SIE_cd, color='y', linestyle=':')

plt.legend()



plt.savefig('fluxratioinversecomplete.eps', dpi=150, bbox_inches='tight')
