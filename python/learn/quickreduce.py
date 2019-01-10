#from pyraf import iraf
#import pyds9
import imexam
#from os import system
#system("/Applications/ds9.darwinsierra.7.5/ds9 &")
# if this does not work, consult https://imexam.readthedocs.io/en/latest/imexam/description.html#how-to-install or perhaps do from outside source activate iraf27
viewer=imexam.connect()
# viewer=imexam.connect("/tmp/xpa/DS9_ds9.1979") # this myight work if the above does not, where the address is from ds9 > File > XPA > Information... > XPA_METHOD
viewer.load_fits('/Users/cerusu/Desktop/IRCA00335317.fits')
# this might work for croscorrelation, especially if I choose image subsection:
from astropy.io import fits
from image_registration import chi2_shift
from image_registration.fft_tools import shift
data1 = fits.open('/Users/cerusu/Desktop/IRCA00335317.fits')[0].data
data2 = fits.open('/Users/cerusu/Desktop/IRCA00335303.fits')[0].data
xoff, yoff, exoff, eyoff = chi2_shift(data1, data2, return_error=True, upsample_factor='auto')

#from stellarpop import distances
#Distance=distances.Distance()
#Distance.reset
#Distance.age(1)
