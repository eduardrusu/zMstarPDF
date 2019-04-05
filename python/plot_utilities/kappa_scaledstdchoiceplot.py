# uses the output of kappa_medsigsim.py to decide the weight-based combination

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/scaledstdchoiceplot.dat"
names = np.genfromtxt(file,usecols=[0],dtype='S200')
data = np.loadtxt(file,usecols=[1,2,3,4,5,6,19,20],dtype='d,d,d,d,d,d,f,f',unpack=True)
medofstd,rmsofstd = data[6],data[7]
N = len(data[0])
los = np.zeros(N)

for i in range(N):
    los[i] = np.log10(np.median([data[0][i],data[1][i],data[2][i],data[3][i],data[4][i],data[5][i]]))
    trim = np.char.find(names[i], '+22.5')
    names[i] = names[i][:trim]

ind = 1 * np.arange(N)  # the x locations for the groups
width = 0.9      # the width of the bars

ax = plt.subplot(2,1,1)
#ax.set(aspect=15)
col1 = los
rects1 = ax.bar(ind + width, col1, width,color='gray')
#ax.set_ylim([0,7])
ax.set_xlim([0.5,N+0.5])
ax.set_ylabel('log(No. of LOS)',fontsize=14)
ax.set_xticks(ind + width)
ax.set_xticklabels([])
plt.title('Relative scaled widths of selected $\kappa_{ext}^{med}-\kappa_{true}$ distributions',fontsize=18)

ax = plt.subplot(2,1,2)
#ax.set(aspect=15)
col2 = medofstd
col2err = rmsofstd
rects2 = ax.bar(ind + width, col2, width, yerr=rmsofstd,color='gray')
#ax.set_ylim([0,1.3])
ax.set_xlim([0.5,N+0.5])
ax.set_ylabel('scaled std($\kappa_{ext}^{med} - \kappa_{true}$)',fontsize=14)
#ax.set_ylabel('$\sigma_{\kappa_\mathrm{med} - \kappa_\mathrm{true}}$')
ax.set_xticks(ind + width)
ax.set_xticklabels(names, fontsize=8, rotation='vertical')
#ax.set_aspect(.0)
#w, h = figaspect(2.)
#fig = Figure(figsize=(w, h))
fsannotate=14
plt.plot([0,N],[1,1],linewidth=1, color='k', linestyle='--')
width = 0.35*2
ax.annotate('1', xy=(1./N-0.0025, -0.82), xytext=(1./N-0.0025, -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*(5-0.5)
ax.annotate('2', xy=(1./N-0.0025+(1./N+5.0/(2*N))+0.0015, -0.82), xytext=(1./N-0.0025+(1./N+5.0/(2*N))+0.0015, -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*(8-0.7)
ax.annotate('3', xy=(1./N-0.0025+(1./N+5.0/(1*N))+8.0/(2*N)+0.000, -0.82), xytext=(1./N-0.0025+(1./N+5.0/(1*N))+8.0/(2*N)+0.000, -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*4-0.1
ax.annotate('4', xy=(1./N-0.0025+(1./N+13.0/N)+4.0/(2*N), -0.82), xytext=(1./N-0.0025+(1./N+13.0/N)+4.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*1-0.1
ax.annotate('5', xy=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(2*N), -0.82), xytext=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*11-0.2
ax.annotate('6', xy=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(2*N), -0.82), xytext=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*6-0.2
ax.annotate('7', xy=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(2*N), -0.82), xytext=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*16-0.3
ax.annotate('8', xy=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(1*N)+16.0/(2*N), -0.82), xytext=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(1*N)+16.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
width = 0.35*16
ax.annotate('9', xy=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(1*N)+16.0/(1*N)+17.0/(2*N), -0.82), xytext=(1./(1*N)-0.0025+13.0/(1*N)+5.0/(1*N)+1.0/(1*N)+11.0/(1*N)+6.0/(1*N)+16.0/(1*N)+17.0/(2*N), -0.92), xycoords='axes fraction',
            fontsize=fsannotate, ha='center', va='bottom'
            #,bbox=dict(boxstyle='square', fc='white')
            ,arrowprops=dict(arrowstyle='-[, widthB=%f, lengthB=13.2' %width, lw=1.0)
            )
plt.subplots_adjust(left=0.08, bottom=0.3, right=0.99, top=0.95, wspace=0.7, hspace=0)
plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/WFI2033/kappasim/scaledstdchoiceplot.png', dpi=250)
plt.clf()
