# CE Rusu, July 22 2018
# Reads the necessary columns form files like W4p2m2_24galphotmstar.cat and converts them to a FITS file of same data type for each column
# run as python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/CFHTLENS_ascii_to_fits.py /Users/cerusu/Dropbox/Davis_work/code/J1206/W4p2m2_24galphotmstar.cat

import fitsio # https://github.com/esheldon/fitsio
import astropy.table as table

file = str(sys.argv[1])
data = np.genfromtxt(file, dtype=None, usecols=[1,2,4,14,16,20,22,24,25], unpack=True)
t = table.Table(data, names=('RA', 'DEC', 'gamma1', 'gamma2', 'w_gal_%s' % limmag, 'w_zweight_%s' % limmag, 'w_oneoverr_%s' % limmag, 'w_zoverr_%s' % limmag, 'galinner_%s' % limmag), dtype=(np.int32,np.float32,np.float32,np.float32,np.int32,np.float32,np.float32,np.float32,np.int32))
t.write('%s.fits' %file[:-4], overwrite=True)


tableout = table.Table(cellkappagammafinal, names=('ID', 'kappa', 'gamma1', 'gamma2', 'w_gal_%s' % limmag, 'w_zweight_%s' % limmag, 'w_oneoverr_%s' % limmag, 'w_zoverr_%s' % limmag, 'galinner_%s' % limmag), dtype=(np.int32,np.float32,np.float32,np.float32,np.int32,np.float32,np.float32,np.float32,np.int32))
    
    
#end = False
#i = 0
#while end == False:
#    try:
#        x = np.genfromtxt(file, dtype=None, usecols=[i])
#        if i == 0: data = x
#        else: data = np.c_[data,x]
#        i += 1
#    except:
        end = True
#print i,end
#t = table.Table(data.T)
#t.write('%s.fits' %file[:-4], overwrite=True)
