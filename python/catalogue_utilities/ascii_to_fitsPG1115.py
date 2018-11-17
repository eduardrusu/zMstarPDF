# CE Rusu, July 22 2018
# Reads an ascii file with multiple columns and converts it into a FITS file with same data types and header. Optionally outut only selected columns and conditions.

import sys
import numpy as np
#import fitsio # https://github.com/esheldon/fitsio
import astropy.table as table
import glob
import time
import collections

#root = '/lfs08/rusucs/HE0435/'
root = '/Volumes/LaCieDavis/CFHTlens/'
rootout = '/Volumes/LaCieDavis/CFHTcatalogues/'
#list = glob.glob(root+'nobetaave3435NEWMEASUREDmedinject_ugriz_HE0435_GGL_los_8_*_120.cat')
listfile = '/Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/fields.cat'
list = np.genfromtxt(listfile,dtype='S')
start_time = time.time()
#for i in range(1):
for i in range(len(list)):
    print list[i]
    #list=['W3p2p3']
    headfile = open(root+list[i]+'.cat', 'r')
    head1 = headfile.readline() # first line
    head2 = headfile.readline() # second line
    head = False
    if (head1.split()[0][0] == '#' and len(head1.split()[0:]) == len(head2.split())) or (head1.split()[0] == '#' and len(head1.split()[1:]) == len(head2.split())): # if there is a header
        if head1.split()[0][0] == '#' and len(head1.split()[0:]) == len(head2.split()):
            head1 = head1.split()[0:] # ignore the # character
            head1[0] = head1[0][1:]
        else: head1 = head1.split()[1:]
        head = True
    dict = collections.OrderedDict()
    for j in range(len(head2.split())):
        if head1[j] in ['MASK','star_flag','MAG_r','MAGERR_r','id','ALPHA_J2000','DELTA_J2000','Flag']: # select desired columns
            if head1[j] == 'id': data = np.genfromtxt(root+list[i]+'.cat',usecols=[j],unpack=True,dtype='S')
            else: data = np.genfromtxt(root+list[i]+'.cat',usecols=[j],unpack=True)
            #if data.dtype == np.float64: type = 'float32'
            #else: type = data.dtype
            type = data.dtype
            if head == True:
                dict[head1[j]] = np.array(data, dtype=type)
            else:
                dict['col%s' %j] = np.array(data, dtype=type)
    del data
    dict['Flag'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data1 = dict['id'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data2 = dict['ALPHA_J2000'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data3 = dict['DELTA_J2000'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data4 = dict['MAG_r'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data5 = dict['MAGERR_r'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    data6 = dict['Flag'][(dict['MAG_r']>0) & (dict['MAG_r']<=24) & (dict['MASK']==0) & (dict['star_flag']==0)]
    dict['id'] = data1; dict['ALPHA_J2000'] = data2; dict['DELTA_J2000'] = data3; dict['MAG_r'] = data4; dict['MAGERR_r'] = data5; dict['Flag'] = data6;
    del dict['MASK']; del dict['star_flag']
    del data1; del data2; del data3; del data4; del data5; del data6;
    t = table.Table(dict)
    del dict
    t.write(rootout+list[i]+'_r24galphot.fits', overwrite = True)
    del t
print(" Total time --- %s seconds ---" % (time.time() - start_time))
