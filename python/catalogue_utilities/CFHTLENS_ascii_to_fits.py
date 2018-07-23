# CE Rusu, July 22 2018
# Reads the necessary columns form files like W4p2m2_24galphotmstar.cat and converts them to a FITS file of same data type for each column
# run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/CFHTLENS_ascii_to_fits.py /Users/cerusu/Dropbox/Davis_work/code/J1206/W4p2m2_24galphotmstar.cat

import sys
import numpy as np
import fitsio # https://github.com/esheldon/fitsio
import astropy.table as table

file = str(sys.argv[1])
data = np.genfromtxt(file, dtype=None, usecols=[1,2,4,14,16,20,22], unpack=True)
t = table.Table(data.T, names=('RA', 'DEC', 'photoz', 'i', 'y', 'mass_BEST', 'mass_MED'), dtype=(np.float64,np.float64,np.float32,np.float32,np.int32,np.float32,np.float32))
t.write('%s.fits' %file[:-4], overwrite=True)
    
#end = False
#i = 0
#while end == False:
#    try:
#        x = np.genfromtxt(file, dtype=None, usecols=[i])
#        if i == 0: data = x
#        else: data = np.c_[data,x]
#        i += 1
#    except:
#        end = True
#print i,end
#t = table.Table(data.T)
#t.write('%s.fits' %file[:-4], overwrite=True)
