# Given a catalogue of RA and DEC, produces the input lists for DAS Cutout and for the PSF cutout
# run as python HSClistforcutouts.py /Volumes/LaCieSubaru/Gaia/James/SecrestHSC5arcsecS17A.fits # before running I crossmatched the catalogue with itself to check that there are no multiple detections inside 5 arcsec

import os
import sys
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

cat = str(sys.argv[1])
s = cat.split('/')[:-1]
path = ''
for i in range(len(s)):
    path = path + s[i] + '/'

from astropy.table import Table
t = Table.read(cat)
coord = np.c_[t['ra'],t['dec']].T
head = '#? rerun      filter    ra       dec       sw     sh  # column descriptor\n'
headpsf = '#? ra          dec      filter  type    rerun centered\n'
rerun = 's18a_wide'
filters = ['HSC-G','HSC-R','HSC-I','HSC-Z','HSC-Y']
filterspsf = ['g','r','i','z','y']
sw = '3.5asec'
sh = '3.5asec'
descriptor = '#'

len1 = np.shape(coord[0])[0] * len(filters) / 1000
len2 = np.shape(coord[0])[0] * len(filters) % 1000
if len2 > 0: len1 +=1
len3 = np.shape(coord[0])[0] * len(filters) / len1
len4 = np.shape(coord[0])[0] * len(filters) % len1
len4 = len3 + len4

fout = []
foutpsf = []
for i in range(len1):
    fout.append(path + cat.split('/')[-1][:-5] + '_cutout' + str(i) +'.cat')
    foutpsf.append(path + cat.split('/')[-1][:-5] + '_psfcutout' + str(i) +'.cat')
    os.system("rm -f %s" % fout[i])
    os.system("rm -f %s" % foutpsf[i])

strcoord = []
strcoordpsf = []
for i in range(np.shape(coord[0])[0]):
    x = SkyCoord(np.float(coord[0][i]), np.float(coord[1][i]), unit='deg')
    strcoord.append('{0} {1}'.format(x.ra.to_string(unit=u.hourangle, sep=':', precision=2, pad=True), x.dec.to_string(sep=':', precision=2, alwayssign=True, pad=True)))
    strcoordpsf.append('{0} {1}'.format(x.ra.deg, x.dec.deg))


pos = 0
for i in range(len1):
    f = open(fout[i],'a')
    g = open(foutpsf[i],'a')
    f.write(head)
    g.write(headpsf)
    if i == 0:
        for j in range(len4):
            f.write(rerun + ' ' + filters[pos % len(filters)] + ' ' + strcoord[pos / len(filters)] + ' ' + sw + ' ' + sh + ' ' + descriptor + '\n')
            g.write(strcoordpsf[pos / len(filters)] + ' ' + filterspsf[pos % len(filters)] + ' coadd ' + rerun + ' true \n')
            pos += 1
    if i != 0:
        for j in range(len3):
            f.write(rerun + ' ' + filters[pos % len(filters)] + ' ' + strcoord[pos / len(filters)] + ' ' + sw + ' ' + sh + ' ' + descriptor + '\n')
            g.write(strcoordpsf[pos / len(filters)] + ' ' + filterspsf[pos % len(filters)] + ' coadd ' + rerun + ' true \n')
            pos += 1
    f.close()
    g.close()
