# CE Rusu, July 22 2018
# Reads an ascii file with multiple columns and converts it into a FITS file with same data types and header

import sys
import numpy as np
import fitsio # https://github.com/esheldon/fitsio
import astropy.table as table
import glob
import time
import collections

root = '/lfs08/rusucs/HE0435/'
list = glob.glob(root+'nobetaave3435NEWMEASUREDmedinject_ugriJHK_HE0435_GGL_los_8_*_45.cat')
start_time = time.time()
#for i in range(1):
for i in range(len(list)):
    print list[i]
    headfile = open(list[i], 'r')
    head1 = headfile.readline()
    head2 = headfile.readline()
    head = False
    if head1.split()[0] == '#' and len(head1.split()[1:]) == len(head2.split()): # if there is a header
        head1 = head1.split()[1:] # ignore the # character
        head = True
    dict = collections.OrderedDict()
    for j in range(len(head2.split())):
        data = np.loadtxt(list[i],usecols=[j],unpack=True)
        if data.dtype == np.float64: type = 'float32'
        else: type = data.dtype
        if head == True:
            dict[head1[j]] = np.array(data, dtype=type)
        else:
            dict['col%s' %j] = np.array(data, dtype=type)
    del data
    t = table.Table(dict)
    del dict
    t.write('%s.fits' %list[i][:-4], overwrite = True)
    del t
print(" Total time --- %s seconds ---" % (time.time() - start_time))
