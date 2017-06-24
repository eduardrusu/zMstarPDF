##########################
# The code finds more realistic error bars for the magnitudes in a catalogue produced with Sextractor
# It does this by finding a beter estimate of the background variance, which accounts for pixel correlations
# The empirical method is based on Liu et al. 2016 (arXiv:1612.01101v1) and references within
##########################

# IF THE CODE FAILS WITH ERROR "ValueError: zero-size array to reduction operation maximum which has no identity", EITHER RUN REPEATEDLY UNTIL IT WORKS, OR REDUCE THE VALUE OF success

# VERY IMPORTANT TO REMEMBER THAT IN DUAL MODE SEXTRACTOR PRODUCES STD AND BACKGROUND MAPS ONLY FOR THE DETECTION IMAGE, SO I FIRST NEED TO RUN SEXTRACTOR IN SINGLE IMAGE MODE AND SAVE THOSE STD AND BACKGROUND MAPS INSTEAD. ALSO, SEXTRACTOR NEEDS TO BE RUN WITH A CUSTOM DEFAULT.PARAM

import numpy as np
import pylab as plt
from astropy.io import fits
from scipy.optimize import leastsq

#image = "FINALweighted_r_covernolens.fits"
#wht = "FINALweighted_r_covernolens_matchg_wht.fits"
#segmentation = "r_noconv_segm.fits"
#back = fits.open("r_noconv_backUSE.fits")[0].data
#std_im = fits.open("r_noconv_stdUSE.fits")[0].data
#savecat = "r_detectin_ir_noconv_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("r_detectin_ir_noconv.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_i_covernolens.fits"
#wht = "FINALweighted_i_covernolens_matchg_wht.fits"
#segmentation = "i_noconv_segm_detect_i.fits"
#back = fits.open("i_noconv_backUSE_detect_i.fits")[0].data
#std_im = fits.open("i_noconv_stdUSE_detect_i.fits")[0].data
#savecat = "i_detectin_i_noconv_newisoerr.cat"
#xpos,ypos,F,err,area = np.loadtxt("i_detectin_i_noconv.cat",usecols=[0,1,10,15,9],unpack=True)

image = "FINALweighted_i_covernolens.fits"
wht = "FINALweighted_i_covernolens_matchg_wht.fits"
segmentation = "i_noconv_segm_detect_ir.fits"
back = fits.open("i_noconv_backUSE_detect_i.fits")[0].data
std_im = fits.open("i_noconv_stdUSE_detect_i.fits")[0].data
savecat = "i_detectin_ir_noconv_newisoerr.cat"
xpos,ypos,F,err,area = np.loadtxt("i_detectin_ir_noconv.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_u_covernolens_matchg_bkgsubt.fits"
#wht = "FINALweighted_u_covernolens_matchg_wht.fits"
#segmentation = "u_segm.fits"
#back = fits.open("u_backUSE.fits")[0].data
#std_im = fits.open("u_stdUSE.fits")[0].data
#savecat = "u_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("u_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_g_covernolens.fits"
#wht = "FINALweighted_g_covernolens_wht.fits"
#segmentation = "g_segm.fits"
#back = fits.open("g_backUSE.fits")[0].data
#std_im = fits.open("g_stdUSE.fits")[0].data
#savecat = "g_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("g_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_r_covernolens_matchg.fits"
#wht = "FINALweighted_r_covernolens_matchg_wht.fits"
#segmentation = "r_segm.fits"
#back = fits.open("r_backUSE.fits")[0].data
#std_im = fits.open("r_stdUSE.fits")[0].data
#savecat = "r_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("r_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_i_covernolens_matchg.fits"
#wht = "FINALweighted_i_covernolens_matchg_wht.fits"
#segmentation = "i_segm.fits"
#back = fits.open("i_backUSE.fits")[0].data
#std_im = fits.open("i_stdUSE.fits")[0].data
#savecat = "i_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("i_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_z_covernolens_matchg.fits"
#wht = "FINALweighted_z_covernolens_matchg_wht.fits"
#segmentation = "z_segm.fits"
#back = fits.open("z_backUSE.fits")[0].data
#std_im = fits.open("z_stdUSE.fits")[0].data
#savecat = "z_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("z_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALweighted_Y_covernolens_matchg.fits"
#wht = "FINALweighted_Y_covernolens_matchg_wht.fits"
#segmentation = "Y_segm.fits"
#back = fits.open("Y_backUSE.fits")[0].data
#std_im = fits.open("Y_stdUSE.fits")[0].data
#savecat = "Y_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("Y_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALHEADmedian_J_covernolens_matchg_bkgsubt.fits"
#wht = "FINALweighted_J_covernolens_matchg_wht.fits"
#segmentation = "J_segm.fits"
#back = fits.open("J_backUSE.fits")[0].data
#std_im = fits.open("J_stdUSE.fits")[0].data
#savecat = "J_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("J_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALHEADmedian_H_covernolens_matchg_bkgsubt.fits"
#wht = "FINALweighted_H_covernolens_matchg_wht.fits"
#segmentation = "H_segm.fits"
#back = fits.open("H_backUSE.fits")[0].data
#std_im = fits.open("H_stdUSE.fits")[0].data
#savecat = "H_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("H_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)

#image = "FINALHEADmedian_Ks_covernolens_matchg_bkgsubt.fits"
#wht = "FINALweighted_Ks_covernolens_matchg_wht.fits"
#segmentation = "Ks_segm.fits"
#back = fits.open("Ks_backUSE.fits")[0].data
#std_im = fits.open("Ks_stdUSE.fits")[0].data
#savecat = "Ks_detectin_ir_newisoerror.cat"
#xpos,ypos,F,err,area = np.loadtxt("Ks_detectin_ir.cat",usecols=[0,1,10,15,9],unpack=True)



img = fits.open(image)[0].data
gain = fits.open(image)[0].header['gain']
segm = fits.open(segmentation)[0].data
aper_img = fits.open(segmentation) # will be used to register the apertures
med_std = np.median(std_im) # median value of the std computed by Sextractor across the image
aper = aper_img[0].data
weight = fits.open(wht)[0].data
Nmax = 20
apertures_nr = 300
gauss = np.zeros([Nmax,3]) # scale, mean and std for the fitted gaussians for N in range (1,20+1)
plt.clf()

for N in range(1,Nmax+1): # sqrt(Npix)
    i = 1 # aperture number, including bad apertures
    success = 0 # number of successful apertures
    aper[aper != 0] = 0 # initialize with no apertures
    flux = np.array([]) # the flux measured in each aperture
    while i <= 3000 and success < apertures_nr:
        ok = True
        xcenter,ycenter = np.random.randint(low=0, high=img.shape[0], size=2)
        xlow = int(xcenter - N/2.0)
        xhigh = int(xcenter + N/2.0)
        ylow = int(ycenter - N/2.0)
        yhigh = int(ycenter + N/2.0)
        # for simplicity, make all apertures squares
        # do not use the aperture if it is not contained inside the image borders, or if it corresponds to non-zero segmentation map values, or if it overlaps with a previous aperture, or if it falls in a region with zero weight
        if xlow < 0 or ylow < 0 or xhigh > img.shape[0] or yhigh > img.shape[0] or np.max(segm[xlow:xhigh,ylow:yhigh] > 0) == True or np.max(aper[xlow:xhigh,ylow:yhigh] > 0) == True or np.max(weight[xlow:xhigh,ylow:yhigh] == 0) == True:
            ok = False
        if ok == True:
            aper[xlow:xhigh,ylow:yhigh] = 1
            # here I account for the fact that the noise is not uniform across the image, and I read it from the Sextractor noise map
            scaledflux = np.sum(img[xlow:xhigh,ylow:yhigh] - back[xlow:xhigh,ylow:yhigh]) * std_im[xcenter][ycenter] / med_std
            flux = np.append(flux,np.array(scaledflux))
            success += 1
        i += 1
    # fitting gaussian to the data
    fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
    errfunc  = lambda p, x, y: (y - fitfunc(p, x))
    init = [100.0, -1, 5]
    ydata,xdatabin = np.histogram(flux)
    bin = xdatabin[1]-xdatabin[0]
    xdata = np.linspace(xdatabin[0]+bin/2.0,xdatabin[len(xdatabin)-1]-bin/2.0,len(xdatabin)-1)
    out = leastsq( errfunc, init, args=(xdata, ydata))
    c = out[0]
    c[2] = np.abs(c[2]) # because by definition, c[2] may be negative
    gauss[N - 1] = c
    #plt.plot(xdata, fitfunc(c, xdata))
    #plt.hist(flux)
    #plt.title(r'$N = %.0f\  A = %.3f\  \mu = %.3f\  \sigma = %.3f\ $' %(N,c[0],c[1],c[2]));
    #plt.show()
#aper_img.writeto("check_aper.fits",clobber=True)

# fitting the power law
std_0 = gauss[0][2]
fitfunc  = lambda p, x: std_0 * (x ** p[0])
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
init = [1.5]
ydata,xdata = gauss[:,2], range(1,Nmax+1)
out = leastsq( errfunc, init, args=(xdata, ydata))
c = out[0]
plt.yscale('log')
plt.plot(xdata, fitfunc(c, xdata))
plt.plot(xdata, fitfunc(np.array([1]), xdata),'r--')
plt.plot(xdata, fitfunc(np.array([2]), xdata),'r--')
plt.scatter(xdata,ydata)
plt.xlabel('N')
plt.ylabel('$\sigma$ [counts/s]')
plt.xlim(0,Nmax + 1)
plt.title(r'$\sigma = %.3f\ \beta = %.3f\ $' %(std_0,c[0]))
plt.show()

# correct the errors
std = np.zeros(len(xpos))
err_recover = np.zeros(len(xpos))
err_new = np.zeros(len(xpos))
for i in range(len(xpos)):
    std[i] = std_im[int(ypos[i])][int(xpos[i])] # switched xpos and ypos because this is how fits images are read to the array
    err_recover[i] = 1.0857 * np.sqrt(area[i] * (std[i] ** 2) + F[i]/gain) / F[i]
    err_new[i] = 1.0857 * np.sqrt((std[i] ** 2) * (area[i] ** c[0]) + F[i]/gain) / F[i]

z=np.linspace(0,0.3,1000)
plt.plot(z,z)
plt.xlim(0,0.5)
plt.ylim(0,0.5)
plt.scatter(err,err_recover,label='recovered')
plt.scatter(err,err_new,color="red",label='corrected')
plt.xlabel('sextractor error')
plt.ylabel('computed error')
plt.legend()
plt.show()

# use the maximum errors between original and recalculated, and replace the NaN values
err_recover[err_recover > 99] = 99.0000
err_new[err_new > 99] = 99.0000
err_recover = np.maximum(err,err_recover)
err_new = np.maximum(err,err_new)
err_recover[np.isnan(err_recover)] = 99.0000
err_new[np.isnan(err_new)] = 99.0000
np.savetxt(savecat,err_new,fmt='%.4f')

# below is for fitting power law with 2 parameters
'''
# fitting the power law
std_0 = gauss[0][2]
fitfunc  = lambda p, x: std_0 * p[0] * (x ** p[1])
errfunc  = lambda p, x, y: (y - fitfunc(p, x))
init = [1,1.5]
ydata,xdata = gauss[:,2], range(1,Nmax+1)
out = leastsq( errfunc, init, args=(xdata, ydata))
c = out[0]
plt.yscale('log')
plt.plot(xdata, fitfunc(c, xdata))
plt.plot(xdata, fitfunc(np.array([1,1]), xdata),'r--')
plt.plot(xdata, fitfunc(np.array([1,2]), xdata),'r--')
plt.scatter(xdata,ydata)
plt.xlabel('N')
plt.ylabel('RMS [counts/s]')
plt.xlim(0,Nmax + 1)
plt.title(r'$\sigma^2 = %.3f\ \alpha = %.3f\ \beta = %.3f\ $' %(std_0**2,c[0],c[1]))
plt.show()

# read the Sextractor catalogue and correct the errors
ID,xpos,ypos,F,err,area = np.loadtxt('r.cat',usecols=[0,1,2,4,7,8],unpack=True)
var = np.zeros(len(xpos))
err_recover = np.zeros(len(xpos))
err_new = np.zeros(len(xpos))
for i in range(len(xpos)):
    var[i] = rms[int(ypos[i])][int(xpos[i])] # switched xpos and ypos because this is how fits images are read to the array
    err_recover[i] = 1.0857 * np.sqrt(area[i] * var[i] + F[i]/gain)/F[i]
    err_new[i] = 1.0857 * np.sqrt(var[i] * (c[0] ** 2) * (area[i]**c[1]) + F[i]/gain)/F[i]
    '''
