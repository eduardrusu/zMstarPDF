# Given a folder with .fits RGB images and the input list for DAS Cutout, it produces a color-combined RGB using Humvi (https://github.com/drphilmarshall/HumVI/blob/master/examples/Examples.ipynb )
# run as python HSCcolorcutouts.py /Volumes/LaCieSubaru/Gaia/cutout0/ /Volumes/LaCieSubaru/Gaia/DelchambreHSC5arcsecunique_cutout0.cat

#from astropy.io import fits
import numpy as np
import os
import sys
import glob
import humvi
from IPython.display import Image
from astropy.io import fits

imgpath = str(sys.argv[1])
cat = str(sys.argv[2])
os.chdir(imgpath)
files = glob.glob('*.fits')

'''
    The input DAS Cutout catalogue is designed to ask for grizy for each object, in case some of the objects are missing some of the gri frames. The first part of the name of the output images produced by DAS Cutout is the line number in the input list.
'''

list = np.zeros(len(files))
for i in range(len(files)):
    list[i] = int(files[i].split('-')[0])
list=np.sort(list)

imcat = np.genfromtxt(cat,usecols=[1,2,3,4],unpack=True,dtype='string')
imgnr = len(imcat[0]) + 1 # image ID starts at 2

pos_g = 2
pos_r = 3
pos_i = 4
pos_z = 5
pos_y = 6
while pos_y <= imgnr:
    print pos_y,'/',imgnr
    pos = [pos_g,pos_r,pos_i,pos_z,pos_y]
    #print pos
    inside = np.in1d(pos, list)
    color = ''
    colorcount = 0
    if inside[0] == True:
        color += 'g'
        if colorcount < 3: colorcount += 1
    if inside[1] == True:
        color += 'r'
        if colorcount < 3: colorcount += 1
    if inside[2] == True:
        color += 'i'
        if colorcount < 3: colorcount += 1
    if inside[3] == True:
        if colorcount < 3:
            color += 'z'
            colorcount += 1
    if inside[4] == True:
        if colorcount < 3:
            color += 'y'
            colorcount += 1
    out = imgpath + imcat[1][pos_g - 2].replace(':','') + imcat[2][pos_g - 2].replace(':','') + '_' + imcat[3][pos_g - 2] + '_' + color + '.png'
    if colorcount == 1:
        bfile = ''
        gfile = ''
        rfile = ''
        here = False
        i = 0
        while here == False:
            if int(files[i].split('-')[0]) in pos:
                bfile, gfile, rfile, outfile = files[i], files[i], files[i], out
                scales, offset, Q, alpha, masklevel, saturation = (1.0,1.0,1.0), 0.5, 1.0, 0.1, -1.0, 'white'
                humvi.compose(rfile, gfile, bfile, scales=scales, Q=Q, alpha=alpha, masklevel=masklevel, saturation=saturation, offset=offset, backsub=False, vb=True, outfile=outfile)
                Image(filename=outfile,width=400)
                here = True
            i += 1
    if colorcount == 2:
        bfile = ''
        gfile = ''
        rfile = ''
        here1 = False
        here2 = False
        i = 0
        while here1 == False or here2 == False:
            if int(files[i].split('-')[0]) in pos:
                if here1 == False:
                    bfile = files[i]
                    here1 = True
                else:
                    here2 = True
                    rfile = files[i]
            i += 1
        im1 = fits.open(bfile)
        im2 = fits.open(rfile)
        data = im1[1].data
        data = (im1[1].data + im2[1].data)/2.0
        im1.writeto('tmp.fits',clobber=True)
        gfile, outfile = 'tmp.fits', out
        scales, offset, Q, alpha, masklevel, saturation = (1.0,1.0,1.0), 0.0, 2.0, 0.1, -1.0, 'white'
        humvi.compose(rfile, gfile, bfile, scales=scales, Q=Q, alpha=alpha, masklevel=masklevel, saturation=saturation, offset=offset, backsub=False, vb=True, outfile=outfile)
        Image(filename=outfile,width=400)
    if colorcount == 3:
        bfile = ''
        gfile = ''
        rfile = ''
        here1 = False
        here2 = False
        here3 = False
        i = 0
        while here1 == False or here2 == False or here3 == False:
            if int(files[i].split('-')[0]) in pos:
                if here1 == False:
                    bfile = files[i]
                    here1 = True
                else:
                    if here1 == True and here2 == False:
                        here2 = True
                        gfile = files[i]
                    else:
                        if here1 == True and here2 == True:
                            here3 = True
                            rfile = files[i]
            i += 1
        outfile = out
        scales, offset, Q, alpha, masklevel, saturation = (1.0,1.0,1.0), 0.0, 2.0, 0.1, -1.0, 'white'
        humvi.compose(rfile, gfile, bfile, scales=scales, Q=Q, alpha=alpha, masklevel=masklevel, saturation=saturation, offset=offset, backsub=False, vb=True, outfile=outfile)
        Image(filename=outfile,width=400)
    pos_g += 5
    pos_r += 5
    pos_i += 5
    pos_z += 5
    pos_y += 5
    print colorcount, out
os.system("rm -f tmp.fits")
