# Simple plot with custom labels

import matplotlib.pyplot as plt
import numpy as np

plt.clf()

data = np.loadtxt('fluxratio.cat', dtype={'names': ('filter', 'A/B', 'A/C', 'A/D', 'B/C'),'formats': ('S2', 'f4', 'f4', 'f4', 'f4')})
filt = [x[0] for x in data]
b = [x[1] for x in data]
c = [x[2] for x in data]
d = [x[3] for x in data]
bc = [x[4] for x in data]

x=np.linspace(1,len(filt),len(filt))
plt.xticks(x, filt)
plt.plot(x, b, label=data.dtype.names[1])
plt.plot(x, c, label=data.dtype.names[2])
plt.plot(x, d, label=data.dtype.names[3])
plt.plot(x, bc, label=data.dtype.names[4])
plt.legend()

plt.savefig('fluxratio.eps', dpi=150, bbox_inches='tight')

plt.clf()

data = np.loadtxt('fluxratio.cat', dtype={'names': ('filter', 'B/A', 'C/A', 'D/A', 'C/B'),'formats': ('S2', 'f4', 'f4', 'f4', 'f4')})
filt = [x[0] for x in data]
b = [1/x[1] for x in data]
c = [1/x[2] for x in data]
d = [1/x[3] for x in data]
bc = [1/x[4] for x in data]

x=np.linspace(1,len(filt),len(filt))
plt.xticks(x, filt)
plt.plot(x, b, label=data.dtype.names[1])
plt.plot(x, c, label=data.dtype.names[2])
plt.plot(x, d, label=data.dtype.names[3])
plt.plot(x, bc, label=data.dtype.names[4])
plt.legend()

plt.savefig('fluxratioinverse.eps', dpi=150, bbox_inches='tight')
