# run as: source activate iraf27; python /Users/cerusu/GITHUB/zMstarPDF/python/reduction_utilities/IRCSquickreduce.py /Users/cerusu/Desktop/IRCS/ APM08279+5255

import os
import numpy as np
from numpy import inf
from astropy.io import fits
import numpy as np
import sys
import glob

os.system('/Applications/ds9.darwinsierra.7.5/ds9 &')
path = str(sys.argv[1])
target = str(sys.argv[2])
try: os.mkdir(path+target)
except: pass
os.chdir(path+target)
files = glob.glob('../*.fits')
filesuse = []
for i in range(len(files)):
    frame = fits.open(files[i])
    if (frame[0].header['OBJECT'] == target) and (files[i] != '../flat.fits'): filesuse = np.append(filesuse,files[i])
np.savetxt(target+'.cat',filesuse,fmt='%s')
os.system('python /Users/cerusu/GITHUB/pyircs_imgred/frcheck.py %s' % target+'.cat')
createflat = False
if createflat == True:
    os.system('python /Users/cerusu/GITHUB/pyircs_imgred/imgred_all.py %s %s.fits --combine=median --skyflat --flat=flat.fits --bpm=/Users/cerusu/GITHUB/pyircs_imgred/DATA/ircs_bpmask.fits --start=0 --end=8' %(target+'.cat',target))
    os.system('cp flat.fits ../')
else:
    os.system('python /Users/cerusu/GITHUB/pyircs_imgred/imgred_all.py %s %s.fits --combine=median --flat=../flat.fits --bpm=/Users/cerusu/GITHUB/pyircs_imgred/DATA/ircs_bpmask.fits --start=0 --end=8' %(target+'.cat',target))
