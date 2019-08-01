##########################
# Simple code for scatter plot with y error bars.
##########################

#RXJ Rein=1.64 Reff=1.85 -> Rein/Reff=0.89  H0=78.2+/-3.4
#J Rein=1.25 Reff=0.34 -> Rein/Reff=3.68  H0=68.9+/-5.2 (Simon confirmed it is circularized)
#HE Rein=1.18 Reff=1.33 -> Rein/Reff=0.89  H0=71.7+/-4.6 (Ken confirmed circularized)
#B Rein=0.85 (main galaxy) Reff=0.58 -> Rein/Reff=1.47 H0=71.0+/-3.1 (F814W. since the lens is spiral, filter is important)
#WFI Rein=0.93 Reff=1.31, 1.44 -> Rein/Reff=0.71-0.0.65 H0=71.6+/-4.4
#PG Rein=1.08 Reff=0.64  H0=81.1+/-7.5 (Chih-Fan)

import pylab as plt
import numpy as np
ax=plt.subplot(111)

fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)
#ax1.set(aspect=1)
ax1.set_aspect(1, adjustable='datalim')

ax = plt.subplot(1,1,1, sharex=ax1, sharey=ax1)
Bh = 71.0; Bhsup = 2.9; Bhinf = 3.3; Bein = 0.85; Beff = 0.58
Rh = 78.2; Rhsup = 3.4; Rhinf = 3.4; Rein = 1.64; Reff = 1.85
Hh = 71.7; Hhsup = 4.8; Hhinf = 4.5; Hein = 1.18; Heff = 1.33
Jh = 68.9; Jhsup = 5.4; Jhinf = 5.1; Jein = 1.25; Jeff = 0.34
Wh = 71.6; Whsup = 3.8; Whinf = 4.9; Wein = 0.93; Weff = 1.31
Ph = 81.1; Phsup = 8.0; Phinf = 7.1; Pein = 1.08; Peff = 0.64
x = np.array([Bein/Beff, Rein/Reff, Hein/Heff, Jein/Jeff, Wein/Weff, Pein/Peff])
y = np.array([Bh, Rh, Hh, Jh, Wh, Ph])
ysup = np.array([Bhsup, Rhsup, Hhsup, Jhsup, Whsup, Phsup])
yinf = np.array([Bhinf, Rhinf, Hhinf, Jhinf, Whinf, Phinf])

#x, y, yinf, ysup = np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz_specflag0_bpz.cat", usecols=(9, 1, 2, 3), unpack=True)
#x, y, yinf, ysup = np.loadtxt("bpzspeczeazy.cat", usecols=(3, 0, 1, 2), unpack=True)
#x = x[abs(y) <= zlim]
#y = y[abs(y) <= zlim]
#yinf = yinf[abs(y) <= zlim]
#ysup = ysup[abs(y) <= zlim]
#y = y[abs(x) <= zlim]
#x = x[abs(x) <= zlim]
#yinf = yinf[abs(x) <= zlim]
#ysup = ysup[abs(x) <= zlim]
ax.tick_params(labelsize=14)
#plt.scatter(x,y, color='k')
plt.errorbar(x[0], y[0], yerr=[[yinf[0]], [ysup[0]]], fmt='o', label='B1608')
plt.errorbar(x[1], y[1], yerr=[[yinf[1]], [ysup[1]]], fmt='o', label='RXJ1131')
plt.errorbar(x[2], y[2], yerr=[[yinf[2]], [ysup[2]]], fmt='o', label='HE0435')
plt.errorbar(x[3], y[3], yerr=[[yinf[3]], [ysup[3]]], fmt='o', label='J1206')
plt.errorbar(x[4], y[4], yerr=[[yinf[4]], [ysup[4]]], fmt='o', label='WFI2033')
plt.errorbar(x[5], y[5], yerr=[[yinf[5]], [ysup[5]]], fmt='o', label='PG1115')
#ax.plot(x, m*x + b, '--')
plt.xlabel('Rein/Reff')
#plt.ylabel('phot-z')
plt.ylabel('H_0')
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
#plt.tick_params(labelbottom='off')
#plt.title('HE0435 ugri specz-photz')
plt.legend()
plt.subplots_adjust(bottom=0.1, left =0.2, right=0.9, top=0.90, wspace=0, hspace=0)
#plt.tight_layout()
#fig.text(0.05, 0.5, 'photo-z', ha='center', va='center', size='20', rotation='vertical')
#plt.title('HE0435 ugri specz-photz')
plt.savefig('H0-Rein_over_Reff.png')
