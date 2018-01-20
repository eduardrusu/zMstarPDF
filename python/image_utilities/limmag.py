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

zeropt = 24.74 + 2.5*np.log10(1020)
image ="g.fits"
scale = 0.25

limmag = np.array([])
x = lim_mag(image,zeropt,scale,390,480)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,290,420)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,330,270)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,230,170)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,200,350)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,210,490)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,300,550)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,100,220)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,430,280)
limmag = np.append(limmag,x)
x = lim_mag(image,zeropt,scale,650,90)
limmag = np.append(limmag,x)

print np.median(limmag),np.std(limmag)

