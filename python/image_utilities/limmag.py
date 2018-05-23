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

zeropt = 23.43 #+ 2.5*np.log10(1020) # do I use the exposure time?
#image ="J1206_NIRI_nativescale_weightedmedian.fits"
#scale = 0.1164 # arcsec
#image ="J1206_NIRI_GMOSscale_weightedmedian.fits"
#scale = 0.1618 # arcsec
image ="J1206_NIRI_CFHTLSscale_weightedmedian.fits"
scale = 0.187 # arcsec

limmag = np.array([])
x = lim_mag(image,zeropt,scale,160,830)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,390,960)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,810,830)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,1190,1080)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,880,80)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,980,420)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,130,250)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,470,120)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,760,470)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,530,620)
limmag = np.append(limmag,x)

print np.median(limmag),np.std(limmag)

