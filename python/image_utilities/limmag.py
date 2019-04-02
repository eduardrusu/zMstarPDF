# Calculates limiting mag from a given region in an image

from astropy.io import fits
import numpy as np

def lim_mag(image,zpt,scale,xpix,ypix):
        i=fits.open(image)
        y=i[0].data
        listx=[]
        listy=[]
        for j in range(y.shape[0]):
            for k in range(y.shape[0]):
                if ((j-xpix)*scale)**2+((k-ypix)*scale)**2<=4:
                    listx=listx+[j]
                    listy=listy+[k]
        list=[]
        for x in range(len(listx)):
            list=list+[y[listx[x],listy[x]]]
        lim=zpt-2.5*np.log10(5*np.sqrt(len(list))*np.std(np.asarray(list))) # 5 sigma
        print lim
        return lim

zeropt = 23.57 #+ 2.5*np.log10(1020) # do I use the exposure time?
#image ="J1206_NIRI_nativescale_weightedmedian.fits"
#scale = 0.1164 # arcsec
#image ="J1206_NIRI_GMOSscale_weightedmedian.fits"
#scale = 0.1618 # arcsec
image ="1016cut.fits"
scale = 0.052 # arcsec

limmag = np.array([])
x = lim_mag(image,zeropt,scale,346,33)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,280,236)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,490,110)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,340,160)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,480,70)
limmag = np.append(limmag,x)
#x = lim_mag(image,zeropt,scale,980,420)
#limmag = np.append(limmag,x)
#x = lim_mag(image,zeropt,scale,130,250)
#limmag = np.append(limmag,x)
#x = lim_mag(image,zeropt,scale,470,120)
#limmag = np.append(limmag,x)
#x = lim_mag(image,zeropt,scale,760,470)
#limmag = np.append(limmag,x)
#x = lim_mag(image,zeropt,scale,530,620)
#limmag = np.append(limmag,x)

print np.median(limmag),np.std(limmag)
