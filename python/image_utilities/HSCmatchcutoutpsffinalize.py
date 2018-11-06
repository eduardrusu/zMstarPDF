# Run after running Topcat on the outputs from HSCmatchcutoutpsf.py. This code renames the files which are suitable for CHITAH.
# run as python HSCmatchcutoutpsffinalize.py

#from astropy.io import fits
import numpy as np
import os
import glob
from astropy.wcs import WCS
from astropy.io import fits

path = "/Volumes/LaCieSubaru/Gaia/James/"
os.chdir(path)
# read the outputs from Topcat
# It seems when I have columns of strings I need to read one column at a time, or they are not searchable...
# ... so I cannot do cutoutspsf_g = np.genfromtxt('cutoutspsf_g.cat',usecols=[0,4,6,7],dtype="S90,S90,float,float")
cutoutspsf_g_cutout = np.genfromtxt('cutoutspsf_g.cat',usecols=[0],dtype="S90")
cutoutspsf_g_psf = np.genfromtxt('cutoutspsf_g.cat',usecols=[4],dtype="S90")
cutoutspsf_g_coord = np.loadtxt('cutoutspsf_g.cat',usecols=[6,7])

cutoutspsf_r_cutout = np.genfromtxt('cutoutspsf_r.cat',usecols=[0],dtype="S90")
cutoutspsf_r_psf = np.genfromtxt('cutoutspsf_r.cat',usecols=[4],dtype="S90")
cutoutspsf_r_coord = np.loadtxt('cutoutspsf_r.cat',usecols=[6,7])

cutoutspsf_i_cutout = np.genfromtxt('cutoutspsf_i.cat',usecols=[0],dtype="S90")
cutoutspsf_i_psf = np.genfromtxt('cutoutspsf_i.cat',usecols=[4],dtype="S90")
cutoutspsf_i_coord = np.loadtxt('cutoutspsf_i.cat',usecols=[6,7])

cutoutspsf_z_cutout = np.genfromtxt('cutoutspsf_z.cat',usecols=[0],dtype="S90")
cutoutspsf_z_psf = np.genfromtxt('cutoutspsf_z.cat',usecols=[4],dtype="S90")
cutoutspsf_z_coord = np.loadtxt('cutoutspsf_z.cat',usecols=[6,7])

cutoutspsf_y_cutout = np.genfromtxt('cutoutspsf_y.cat',usecols=[0],dtype="S90")
cutoutspsf_y_psf = np.genfromtxt('cutoutspsf_y.cat',usecols=[4],dtype="S90")
cutoutspsf_y_coord = np.loadtxt('cutoutspsf_y.cat',usecols=[6,7])

psf_onlybandg = np.genfromtxt('psf_onlybandg.cat',usecols=[0],dtype="S90")
psf_onlybandr = np.genfromtxt('psf_onlybandr.cat',usecols=[0],dtype="S90")
psf_onlybandi = np.genfromtxt('psf_onlybandi.cat',usecols=[0],dtype="S90")
psf_onlybandz = np.genfromtxt('psf_onlybandz.cat',usecols=[0],dtype="S90")
psf_onlybandy = np.genfromtxt('psf_onlybandy.cat',usecols=[0],dtype="S90")

folder_cutout = glob.glob('arch*')
list_cutout = np.array([])
for i in range(len(folder_cutout)):
    files = glob.glob('%s/*' % folder_cutout[i])
    files_array = np.asarray(files)
    list_cutout = np.r_[list_cutout,files_array]

os.system("mkdir S17AforJames")
os.system("mkdir S17Aonebandpsf")
os.system("mkdir S17Anopsf")

for i in range(len(list_cutout)):
    if '-G-' in list_cutout[i]:
        if list_cutout[i] in cutoutspsf_g_cutout:
            psf = cutoutspsf_g_psf[np.where(cutoutspsf_g_cutout == list_cutout[i])[0][0]]
            coord = cutoutspsf_g_coord[np.where(cutoutspsf_g_cutout == list_cutout[i])[0][0]]
            if psf not in psf_onlybandg:
                os.system("cp %s S17AforJames/cutout_%.5f_%.5f_G.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17AforJames/psf_%.5f_%.5f_G.fits" % (psf,coord[0],coord[1]))
            else:
                os.system("cp %s S17Aonebandpsf/cutout_%.5f_%.5f_G.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17Aonebandpsf/psf_%.5f_%.5f_G.fits" % (psf,coord[0],coord[1]))
        else: os.system("cp %s S17Anopsf/cutout_%.5f_%.5f_G_%d.fits" % (list_cutout[i],coord[0],coord[1],i)) # adding the number to the name so I don't lose files with identical input
    if '-R-' in list_cutout[i]:
        if list_cutout[i] in cutoutspsf_r_cutout:
            psf = cutoutspsf_r_psf[np.where(cutoutspsf_r_cutout == list_cutout[i])[0][0]]
            coord = cutoutspsf_r_coord[np.where(cutoutspsf_r_cutout == list_cutout[i])[0][0]]
            if psf not in psf_onlybandr:
                os.system("cp %s S17AforJames/cutout_%.5f_%.5f_R.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17AforJames/psf_%.5f_%.5f_R.fits" % (psf,coord[0],coord[1]))
            else:
                os.system("cp %s S17Aonebandpsf/cutout_%.5f_%.5f_R.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17Aonebandpsf/psf_%.5f_%.5f_R.fits" % (psf,coord[0],coord[1]))
        else: os.system("cp %s S17Anopsf/cutout_%.5f_%.5f_R_%d.fits" % (list_cutout[i],coord[0],coord[1],i))
    if '-I-' in list_cutout[i]:
        if list_cutout[i] in cutoutspsf_i_cutout:
            psf = cutoutspsf_i_psf[np.where(cutoutspsf_i_cutout == list_cutout[i])[0][0]]
            coord = cutoutspsf_i_coord[np.where(cutoutspsf_i_cutout == list_cutout[i])[0][0]]
            if psf not in psf_onlybandi:
                os.system("cp %s S17AforJames/cutout_%.5f_%.5f_I.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17AforJames/psf_%.5f_%.5f_I.fits" % (psf,coord[0],coord[1]))
            else:
                os.system("cp %s S17Aonebandpsf/cutout_%.5f_%.5f_I.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17Aonebandpsf/psf_%.5f_%.5f_I.fits" % (psf,coord[0],coord[1]))
        else: os.system("cp %s S17Anopsf/cutout_%.5f_%.5f_I_%d.fits" % (list_cutout[i],coord[0],coord[1],i))
    if '-Z-' in list_cutout[i]:
        if list_cutout[i] in cutoutspsf_z_cutout:
            psf = cutoutspsf_z_psf[np.where(cutoutspsf_z_cutout == list_cutout[i])[0][0]]
            coord = cutoutspsf_z_coord[np.where(cutoutspsf_z_cutout == list_cutout[i])[0][0]]
            if psf not in psf_onlybandz:
                os.system("cp %s S17AforJames/cutout_%.5f_%.5f_Z.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17AforJames/psf_%.5f_%.5f_Z.fits" % (psf,coord[0],coord[1]))
            else:
                os.system("cp %s S17Aonebandpsf/cutout_%.5f_%.5f_Z.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17Aonebandpsf/psf_%.5f_%.5f_Z.fits" % (psf,coord[0],coord[1]))
        else: os.system("cp %s S17Anopsf/cutout_%.5f_%.5f_Z_%d.fits" % (list_cutout[i],coord[0],coord[1],i))
    if '-Y-' in list_cutout[i]:
        if list_cutout[i] in cutoutspsf_y_cutout:
            psf = cutoutspsf_y_psf[np.where(cutoutspsf_y_cutout == list_cutout[i])[0][0]]
            coord = cutoutspsf_y_coord[np.where(cutoutspsf_y_cutout == list_cutout[i])[0][0]]
            if psf not in psf_onlybandy:
                os.system("cp %s S17AforJames/cutout_%.5f_%.5f_Y.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17AforJames/psf_%.5f_%.5f_Y.fits" % (psf,coord[0],coord[1]))
            else:
                os.system("cp %s S17Aonebandpsf/cutout_%.5f_%.5f_Y.fits" % (list_cutout[i],coord[0],coord[1]))
                #os.system("cp %s S17Aonebandpsf/psf_%.5f_%.5f_Y.fits" % (psf,coord[0],coord[1]))
        else: os.system("cp %s S17Anopsf/cutout_%.5f_%.5f_Y_%d.fits" % (list_cutout[i],coord[0],coord[1],i))

os.chdir("/Users/cerusu/GITHUB/zMstarPDF/python/image_utilities")
