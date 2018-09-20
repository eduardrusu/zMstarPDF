# Given folders containing image cutouts and corresponding PSFs, match them by accounting for duplicates and missing files
# run as python HSCmatchcutoutpsf.py

#from astropy.io import fits
import numpy as np
import os
import glob
from astropy.wcs import WCS
from astropy.io import fits

path = "/Volumes/LaCieSubaru/Gaia/James/"
os.chdir(path)
folder_cutout = glob.glob('arch*')
folder_psf = glob.glob('psf-*')

list_cutout = np.array([])
for i in range(len(folder_cutout)):
    files = glob.glob('%s/*' % folder_cutout[i])
    files_array = np.asarray(files)
    list_cutout = np.r_[list_cutout,files_array]

list_psf = np.array([])
for i in range(len(folder_psf)):
    files = glob.glob('%s/*' % folder_psf[i])
    files_array = np.asarray(files)
    list_psf = np.r_[list_psf,files_array]

list_psf_filter = np.empty(len(list_psf),dtype="string")
list_psf_ra = np.zeros(len(list_psf))
list_psf_dec = np.zeros(len(list_psf))
for i in range(len(list_psf)):
    str = list_psf[i].split("-")
    ra = float(str[10])
    if str[-2] == '': str[-1] = '-' + str[-1]
    dec = float(str[-1][:-5])
    filter = str[7]
    list_psf_filter[i] = filter
    list_psf_ra[i] = ra
    list_psf_dec[i] = dec

list_cutout_filter = np.empty(len(list_cutout),dtype="string")
list_cutout_ra = np.zeros(len(list_cutout))
list_cutout_dec = np.zeros(len(list_cutout)) 
for i in range(len(list_cutout)):
    if i % 1000 == 0: print i
    file = fits.open(list_cutout[i])
    w = WCS(file[1].header)
    coord = w.wcs_pix2world(21,21,1)
    list_cutout_ra[i] = coord[0]
    list_cutout_dec[i] = coord[1]
    file.close()
    str = list_cutout[i].split("-")
    filter = str[7]
    list_cutout_filter[i] = filter

out_cutout = np.c_[list_cutout,list_cutout_filter,list_cutout_ra,list_cutout_dec]
out_psf = np.c_[list_psf,list_psf_filter,list_psf_ra,list_psf_dec]
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/cutouts_g.cat",out_cutout[out_cutout[:,1]=='G'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/cutouts_r.cat",out_cutout[out_cutout[:,1]=='R'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/cutouts_i.cat",out_cutout[out_cutout[:,1]=='I'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/cutouts_z.cat",out_cutout[out_cutout[:,1]=='Z'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/cutouts_y.cat",out_cutout[out_cutout[:,1]=='Y'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/psf_g.cat",out_psf[out_psf[:,1]=='G'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/psf_r.cat",out_psf[out_psf[:,1]=='R'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/psf_i.cat",out_psf[out_psf[:,1]=='I'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/psf_z.cat",out_psf[out_psf[:,1]=='Z'],'%s \t %s \t %s \t %s')
np.savetxt("/Volumes/LaCieSubaru/Gaia/James/psf_y.cat",out_psf[out_psf[:,1]=='Y'],'%s \t %s \t %s \t %s')
os.chdir("/Users/cerusu/GITHUB/zMstarPDF/python/image_utilities")
#os.chdir(imgpath)
#os.system("rm -f tmp.fits")
