##########################
# Given a list of spectroscopic redshifts, plots the histograms and marks the groups previously identified
##########################

#from matplotlib.colors import LogNorm
#import scipy.optimize as optimization
from pylab import *
import numpy as np
#ax=plt.subplot(111)
zlim = 1.8
#maglimit = 21
#outlim = 0.15

fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)
#ax1.set(aspect=1)
ax1.set_aspect(1, adjustable='datalim')
ax = plt.subplot(1,1,1, sharex=ax1, sharey=ax1)

z = np.loadtxt("/Users/cerusu/Dropbox/Davis_work/code/J1206/spec/J1206_specfromChrisandSDSSusethis.tab", usecols=[4], unpack=True)
plt.hist(z[z < zlim], bins = 72, color='k', histtype='step')

z1 = 0.052; z1_size = 10; z1_ra = -133; z1_raerr = 178; z1_dec = 368; z1_decerr = 98
z2 = 0.123; z2_size = 6; z2_ra = 59; z2_raerr = 203; z2_dec = -212; z2_decerr = 68
z3 = 0.482; z3_size = 14; z3_ra = -230; z3_raerr = 67; z3_dec = -41; z3_decerr = 93
z4 = 0.551; z4_size = 27; z4_ra = 127; z4_raerr = 98; z4_dec = -171; z4_decerr = 64
z5 = 0.688; z5_size = 16; z5_ra = 92; z5_raerr = 68; z5_dec = -112; z5_decerr = 47
z6 = 0.750; z6_size = 26; z6_ra = -33; z6_raerr = 54; z6_dec = 5; z6_decerr = 26
z7 = 0.778; z7_size = 4; z7_ra = -227; z7_raerr = 30; z7_dec = -94; z7_decerr = 34

plt.axvline(x=z1, color='k', linestyle='--')
plt.axvline(x=z2, color='k', linestyle='--')
plt.axvline(x=z3, color='k', linestyle='--')
plt.axvline(x=z4, color='k', linestyle='--')
plt.axvline(x=z5, color='k', linestyle='--')
plt.axvline(x=z6, color='k', linestyle='--')
plt.axvline(x=z7, color='k', linestyle='--')

#x = x[abs(y) <= zlim]
#y = y[abs(y) <= zlim]
#plt.scatter(x,y, color='k')
#plt.scatter(x, y)
txt1 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z1,z1_size,z1_ra,z1_raerr,z1_dec,z1_decerr)
txt2 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z2,z2_size,z2_ra,z2_raerr,z2_dec,z2_decerr)
txt3 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z3,z3_size,z3_ra,z3_raerr,z3_dec,z3_decerr)
txt4 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z4,z4_size,z4_ra,z4_raerr,z4_dec,z4_decerr)
txt5 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z5,z5_size,z5_ra,z5_raerr,z5_dec,z5_decerr)
txt6 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z6,z6_size,z6_ra,z6_raerr,z6_dec,z6_decerr)
txt7 = "%.3f %s %s$\pm$%s %s$\pm$%s" % (z7,z7_size,z7_ra,z7_raerr,z7_dec,z7_decerr)
#stdout = "scatter = %.3f" % std
plt.text(0.5, 0.9, txt1, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.8, txt2, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.7, txt3, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.6, txt4, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.5, txt5, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.4, txt6, fontsize=9, color='black', transform=ax.transAxes)
plt.text(0.5, 0.3, txt7, fontsize=9, color='black', transform=ax.transAxes)
plt.xlabel('spectroscopic redshift')
plt.ylabel('number of galaxies')
plt.xlim(0, zlim)
#plt.ylim(0, zlim)
#ax1.set_yticklabels(np.arange(0.0, zlim, 0.2))
#ax1.set_xticklabels(np.arange(0.0, zlim, 0.2))
#plt.subplots_adjust(bottom=0.1, left =0.2, right=0.9, top=0.90, wspace=0, hspace=0)
#plt.tight_layout()
#fig.text(0.05, 0.5, 'photo-z', ha='center', va='center', size='20', rotation='vertical')
#plt.title('HE0435 ugri specz-photz')
plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/J1206/speczhist.png')
