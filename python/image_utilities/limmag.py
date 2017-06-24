# Calculates limiting mag from a given region in an image

from astropy.io import fits
import numpy as np

def lim_mag(image,zpt,xpix,ypix):
        i=fits.open(image)
        y=i[0].data
        listx=[]
        listy=[]
        for j in range(y.shape[0]):
            for k in range(y.shape[0]):
                if ((j-xpix)*0.26)**2+((k-ypix)*0.26)**2<=4:
                    listx=listx+[j]
                    listy=listy+[k]
        list=[]
        for x in range(len(listx)):
            list=list+[y[listx[x],listy[x]]]
        lim=zpt-2.5*np.log10(5*np.sqrt(len(list))*np.std(np.asarray(list))) # 5 sigma
        return lim

lim_mag("FINALHEADmedian_Ks.fits",22.90,140,740)
lim_mag("FINALHEADmedian_Ks.fits",22.90,330,770)
lim_mag("FINALHEADmedian_Ks.fits",22.90,610,740)
lim_mag("FINALHEADmedian_Ks.fits",22.90,400,530)
lim_mag("FINALHEADmedian_Ks.fits",22.90,770,560)
lim_mag("FINALHEADmedian_Ks.fits",22.90,320,470)
lim_mag("FINALHEADmedian_Ks.fits",22.90,570,400)
lim_mag("FINALHEADmedian_Ks.fits",22.90,230,360)
lim_mag("FINALHEADmedian_Ks.fits",22.90,120,110)
lim_mag("FINALHEADmedian_Ks.fits",22.90,640,110)
