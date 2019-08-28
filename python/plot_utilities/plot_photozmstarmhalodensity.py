# produces density plot for catalogue vs. computed values of redshift, stellar masses and halo masses

from matplotlib.colors import LogNorm
import scipy.optimize as optimization
from pylab import *
import numpy as np

font = 10
ticksize = 10

plt.clf()
fig = plt.figure(figsize=(10,12))
#fig, axes = plt.subplots(nrows=2, ncols=2)

ax1 = fig.add_subplot(3,2,1)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images_forNAOJ.txt" % (i,j,k), usecols=(1, 8), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 1.66
outlim = 0.15
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
colorbar()
delta = (y-x)/(1+x)
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects; outliers" % (len(x))
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
xx = np.linspace(-0.15, 20, 1000)
xxx = np.linspace(0.15, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(xx, 0.85*xx-0.15, 'r--')
plt.plot(xxx, 1.15*xxx+0.15, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue', fontsize=font)
plt.ylabel('photoz ugrizJHK', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(0, zlim)
plt.ylim(0, zlim)
#plt.title(' ugriJHK catalogue - pdz')



ax1 = fig.add_subplot(3,2,2)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt" % (i,j,k), usecols=(1, 8), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 1.66
outlim = 0.15
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
colorbar()
delta = (y-x)/(1+x)
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % (len(x))
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
xx = np.linspace(-0.15, 20, 1000)
xxx = np.linspace(0.15, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(xx, 0.85*xx-0.15, 'r--')
plt.plot(xxx, 1.15*xxx+0.15, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue', fontsize=font)
plt.ylabel('photoz ugriz', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(0, zlim)
plt.ylim(0, zlim)
#plt.title(' ugriz catalogue - pdz')



ax1 = fig.add_subplot(3,2,3)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images_forNAOJ.txt" % (i,j,k), usecols=(5, 10), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 12
zlim_ = 7
outlim = 0.5
x = np.log10(x)
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
x = x[abs(y) > zlim_]
y = y[abs(y) > zlim_]
y = y[abs(x) > zlim_]
x = x[abs(x) > zlim_]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
colorbar()
delta = (y-x)
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % len(x)
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(x, x+0.5, 'r--')
plt.plot(x, x-0.5, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue $\log M_\star$ ugriJHK', fontsize=font)
plt.ylabel('measured $\log M_\star$ ugriJHK', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(zlim_, zlim)
plt.ylim(zlim_, zlim)
#plt.title(' ugriJHK catalogue - measured Mstar')



ax1=plt.subplot(3,2,4)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt" % (i,j,k), usecols=(5, 9), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim_ = 7
outlim = 0.5
x = np.log10(x)
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
x = x[abs(y) > zlim_]
y = y[abs(y) > zlim_]
y = y[abs(x) > zlim_]
x = x[abs(x) > zlim_]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
delta = y-x
colorbar()
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % len(x)
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(x, x+0.5, 'r--')
plt.plot(x, x-0.5, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue $\log M_\star$ ugriz', fontsize=font)
plt.ylabel('measured $\log M_\star$ ugriz', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(zlim_, zlim)
plt.ylim(zlim_, zlim)
#plt.title(' ugriz catalogue - measured Mstar')

ax1 = fig.add_subplot(3,2,5)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugrizJHK_WFI2033.images_forNAOJ.txt" % (i,j,k), usecols=(4, 12), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 14.5
zlim_ = 10
outlim = 0.5
x = np.log10(x)
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
x = x[abs(y) > zlim_]
y = y[abs(y) > zlim_]
y = y[abs(x) > zlim_]
x = x[abs(x) > zlim_]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
colorbar()
delta = (y-x)
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % len(x)
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(x, x+0.5, 'r--')
plt.plot(x, x-0.5, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue $\log M_{halo}$ ugriJHK', fontsize=font)
plt.ylabel('measured $\log M_{halo}$ ugriJHK', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(zlim_, zlim)
plt.ylim(zlim_, zlim)
#plt.title(' ugriJHK catalogue - measured Mstar')


ax1=plt.subplot(3,2,6)
ax1.set_aspect(1)
for i in range(7):
    for j in range(4):
        for k in range(4):
            x_, y_ = np.loadtxt("/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/GGL_los_8_0_%s_%s_%s_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt" % (i,j,k), usecols=(4, 10), unpack=True)
            if (i==0) & (j==0) & (k==0):
                x = x_
                y = y_
            else:
                x = np.append(x,x_)
                y = np.append(y,y_)
zlim = 14.5
zlim_ = 10
outlim = 0.5
x = np.log10(x)
x = x[abs(y) <= zlim]
y = y[abs(y) <= zlim]
y = y[abs(x) <= zlim]
x = x[abs(x) <= zlim]
x = x[abs(y) > zlim_]
y = y[abs(y) > zlim_]
y = y[abs(x) > zlim_]
x = x[abs(x) > zlim_]
hist2d(x, y, bins=[100, 100], norm=LogNorm())
plt.xticks(rotation='vertical',size = ticksize)
plt.yticks(size = ticksize)
#m, b = np.polyfit(x, y, 1)
def func(x, m):
    return m+x
delta = y-x
colorbar()
#delta = (y-(m+x))/(1+x)
#firstpass_std = np.std(delta)
firstpass_std = np.std(delta)
std = np.std(delta[abs(delta)<outlim])
stdoutobj = "%d objects" % len(x)
stdout = "scatter = %.3f" % std
#out = 100.0*len(delta[abs(delta)>(4 * std)])/len(delta)
outnr = len(delta[abs(delta)>outlim])
out = 100.0*len(delta[abs(delta)>outlim])/len(delta)
outlier = "outliers = %d (%.2f %%)" % (outnr,out)
#bestfit = "best-fit line = %.2f * x + %.2f" % (m,b)
#m, b = np.polyfit(x, y, 1)
x = x[np.where(abs(delta)<outlim)]
y = y[np.where(abs(delta)<outlim)]
def func(x, m):
    return m+x
fit=optimization.curve_fit(func, x, y, np.array([0])) # the array contains the initial guesses
m=fit[0][0]
bias = "bias = %.3f" % m
x = np.linspace(0, 20, 1000)
#plt.plot(x, m*x + b, '--')
plt.plot(x, m+x, '--')
plt.plot(x, x+0.5, 'r--')
plt.plot(x, x-0.5, 'r--')
plt.plot(x, x)
ax1.text(0.05, 0.95, stdoutobj, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.90, stdout, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.85, outlier, fontsize=7, color='black',transform=ax1.transAxes)
ax1.text(0.05, 0.80, bias, fontsize=7, color='black',transform=ax1.transAxes)
#ax1.text(0.05, 0.80, bestfit, fontsize=7, color='black',transform=ax1.transAxes)
plt.xlabel('catalogue $\log M_{halo}$ ugriz', fontsize=font)
plt.ylabel('measured $\log M_{halo}$ ugriz', fontsize=font)
#plt.xlim(0, 3.5)
#plt.ylim(0, 3.5)
plt.xlim(zlim_, zlim)
plt.ylim(zlim_, zlim)
#plt.title(' ugriz catalogue - measured Mstar')

#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
#cb = plt.colorbar(ax1, cax = cbaxes)

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(im, cax=cbar_ax)

#plt.subplots_adjust(left=0.1, bottom=0.1, right=0.80, top=0.90, wspace=0.4, hspace=0.4)
#plt.tight_layout()
plt.savefig('/Volumes/perseus_1/simulations/lensing_simulations/SA_galaxies/original/WFI2033/WFI2033_photoz_Mstar_Mhalo.png' , dpi=250)
