# CE Rusu, Feb 14 2018
# This code makes use of the CFHTLens *galphotmstar.cat files and the lens photometric+Mstar+photoz catalogue; it computes weighted ratios (lens/field) with proper masking, for various radii, limiting mag, number of samples and classification scheme.
# run as: python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_overlap_sampling_nobeta.py WFI2033 /Volumes/LaCieDavis/CFHTcatalogues/W1m0m0_24galphotmstar.cat /Volumes/LaCieDavis/CFHTLenSmasks/W1m0m0_izrgu_finalmask_mosaic.fits /Users/cerusu/Dropbox/Davis_work/code /Volumes/LaCieSubaru/weightedcounts/WFI2033 45 5 IRAC deti meds removegrouphandpicked -1 -1
# the code produces output for two different limiting mags: 23 and 24
# the scripts to run en masse are in /Users/cerusu/GITHUB/zMstarPDF/python/scripts/DESKTOP/
# the code is optimized for speed, but may be memory intensive because it stores the input catalogue in memory
# definitions:
#field refers to one of the 171 CHFT fields
#lens refers to the name of the H0licow lens
#field mask refers to the mask .fits file of the corresponding CHFT field
#cell refers to a rectangular surface of length on a side 2*radius
#radius is the aperture radius around the actual lens object around which the environmant is explored
#irac shows whether IRAC bands are used or not; use either "IRAC" or "noIRAC"
#det characterizes the object detection type; for example, for WFI2033 I detect in i or i+r; use either "detir" or "deti"
#mode (sum or meds): whether the weighted counts are summed or Median * counts is considered
#remove=[filename w/o extension]/no: refers to whether or not some objects specified in the code should be removed from the lens catalogue
#zinf, zsup: the lower and upper redshift gap in case I remove redshift slices; just use a negative range if not removing anything
#overlap: the CFHTLenS field is covered with a grid of cells that can overlap; overlap=2 means that the field is covered with two grids: a first grid of non-overlapping cells, and then another of non-overlapping cells, and the two grids are offset relatively by 1/2 of the width of a cell; overlap is used to increase the number of data points; overlap of 2 increases the data points by a factor of 4
# convergence refers to the effective convergence from Momcheva et al. 2006

import numpy as np
import scipy
import sys
import os
from scipy import ndimage
from scipy import special
from scipy import stats
from scipy.interpolate import griddata
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
#import distances
import time
import pandas as pd

lensID = str(sys.argv[1])
fieldID = str(sys.argv[2])
fmask = str(sys.argv[3])
rootlens = str(sys.argv[4])
output = str(sys.argv[5])
radius = int(str(sys.argv[6]))
inner = int(str(sys.argv[7]))
irac = str(sys.argv[8])
det = str(sys.argv[9])
type = str(sys.argv[10])
remove = str(sys.argv[11])
zinf = float(str(sys.argv[12]))
zsup = float(str(sys.argv[13]))

suffix = '' # added to the name of the output catalogues
samples = 10

#print("Arguments: \n Lens: %s \n Field file: %s \n Field mask file: %s \n Limiting radius: %s \n" % (lensID, fieldID, fmask, radius))

def lensinit(lensbpz_masked_,lens_gal_,lens_zweight_,lens_mass_,lens_mass2_,lens_mass3_,lens_oneoverr_,lens_zoverr_,lens_massoverr_,lens_mass2overr_,lens_mass3overr_,lens_flexion_,lens_tidal_,lens_convergence_,lens_convergencehalo_,lens_mass2rms_,lens_mass3rms_,lens_mass2overrms_,lens_mass3overrms_):
    if (lensbpz_masked_.shape[1] != 0):
        lens_gal_[k][l][i][j][n] = lensbpz_masked_.shape[1]
        #print lens_gal_[k][l][i][j][n] # test to check if the function actually returns the result globally
        lens_zweight = (z_s * lensbpz_masked_[z_lens + n * 3]) - lensbpz_masked_[z_lens + n * 3] ** 2
        lens_mass = 10 ** lensbpz_masked_[Mlens + n * 3]
        lens_mass2 = (10 ** lensbpz_masked_[Mlens + n * 3]) ** 2
        lens_mass3 = (10 ** lensbpz_masked_[Mlens + n * 3]) ** 3
        lens_oneoverr = 1.0/lensbpz_masked_[sep]
        lens_zoverr = lens_zweight * 1.0 / lensbpz_masked_[sep]
        lens_massoverr = lens_mass * 1.0 / lensbpz_masked_[sep]
        lens_mass2overr = lens_mass2 * 1.0 / lensbpz_masked_[sep]
        lens_mass3overr = lens_mass3 * 1.0 / lensbpz_masked_[sep]
        lens_flexion = lens_mass / (lensbpz_masked_[sep] ** 3) # ignoring const 4GM/c^2 and unit conversion and D_PS / (D_P * D_S): M / R^3
        lens_tidal = lens_mass / (lensbpz_masked_[sep] ** 2) # (M * (D_PS / (D_P * D_S)) / R^2) * (1 - D_12 * D_S / (D_P * D_LS))
        convergence = np.sqrt(lens_mass) / lensbpz_masked_[sep] # sqrt(const * M * (D_PS / (D_P * D_S))) / R
        lens_masshalo = 10 ** lensbpz_masked_[Mhalo + n * 3]
        convergencehalo = np.sqrt(lens_masshalo) / lensbpz_masked_[sep] # sqrt(const * M * (D_PS / (D_P * D_S))) / R
        if type == "meds":
            lens_zweight_[k][l][i][j][n] = np.median(lens_zweight) * lens_gal_[k][l][i][j][n]
            lens_mass_[k][l][i][j][n] = np.median(lens_mass) * lens_gal_[k][l][i][j][n]
            lens_mass2_[k][l][i][j][n] = np.median(lens_mass2) * lens_gal_[k][l][i][j][n]
            lens_mass3_[k][l][i][j][n] = np.median(lens_mass3) * lens_gal_[k][l][i][j][n]
            lens_oneoverr_[k][l][i][j][n] = np.median(lens_oneoverr) * lens_gal_[k][l][i][j][n]
            lens_zoverr_[k][l][i][j][n] = np.median(lens_zoverr) * lens_gal_[k][l][i][j][n]
            lens_massoverr_[k][l][i][j][n] = np.median(lens_massoverr) * lens_gal_[k][l][i][j][n]
            lens_mass2overr_[k][l][i][j][n] = np.median(lens_mass2overr) * lens_gal_[k][l][i][j][n]
            lens_mass3overr_[k][l][i][j][n] = np.median(lens_mass3overr) * lens_gal_[k][l][i][j][n]
            lens_flexion_[k][l][i][j][n] = np.median(lens_flexion) * lens_gal_[k][l][i][j][n]
            lens_tidal_[k][l][i][j][n] = np.median(lens_tidal) * lens_gal_[k][l][i][j][n]
            lens_convergence_[k][l][i][j][n] = np.median(convergence) * lens_gal_[k][l][i][j][n]
            lens_convergencehalo_[k][l][i][j][n] = np.median(convergencehalo) * lens_gal_[k][l][i][j][n]
        if type == "sum":
            lens_zweight_[k][l][i][j][n] = np.sum(lens_zweight)
            lens_mass_[k][l][i][j][n] = np.sum(lens_mass)
            lens_mass2_[k][l][i][j][n] = np.sum(lens_mass2)
            lens_mass3_[k][l][i][j][n] = np.sum(lens_mass3)
            lens_oneoverr_[k][l][i][j][n] = np.sum(lens_oneoverr)
            lens_zoverr_[k][l][i][j][n] = np.sum(lens_zoverr)
            lens_massoverr_[k][l][i][j][n] = np.sum(lens_massoverr)
            lens_mass2overr_[k][l][i][j][n] = np.sum(lens_mass2overr)
            lens_mass3overr_[k][l][i][j][n] = np.sum(lens_mass3overr)
            lens_flexion_[k][l][i][j][n] = np.sum(lens_flexion)
            lens_tidal_[k][l][i][j][n] = np.sum(lens_tidal)
            lens_convergence_[k][l][i][j][n] = np.sum(convergence)
            lens_convergencehalo_[k][l][i][j][n] = np.sum(convergencehalo)
        lens_mass2rms_[k][l][i][j][n] = np.sqrt(lens_mass2_[k][l][i][j][n])
        lens_mass3rms_[k][l][i][j][n] = scipy.special.cbrt(lens_mass3_[k][l][i][j][n])
        lens_mass2overrms_[k][l][i][j][n] = np.sqrt(lens_mass2overr_[k][l][i][j][n])
        lens_mass3overrms_[k][l][i][j][n] = scipy.special.cbrt(lens_mass3overr_[k][l][i][j][n])
    return lens_gal_,lens_zweight_,lens_mass_,lens_mass2_,lens_mass3_,lens_oneoverr_,lens_zoverr_,lens_massoverr_,lens_mass2overr_,lens_mass3overr_,lens_flexion_,lens_tidal_,lens_convergence_,lens_convergencehalo_,lens_mass2rms_,lens_mass3rms_,lens_mass2overrms_,lens_mass3overrms_

def fieldinit(field_masked_,w_gal_,field_gal_,field_zweight_,field_mass_,field_mass2_,field_mass3_,field_oneoverr_,field_zoverr_,field_massoverr_,field_mass2overr_,field_mass3overr_,field_mass2rms_,field_mass3rms_,field_mass2overrms_,field_mass3overrms_,field_flexion_,field_tidal_,field_convergence_,field_convergencehalo_):
    if (w_gal_ != 0):
        field_gal_[k][l][i][j] = w_gal_
        field_zweight = (z_s * field_masked_[z_field]) - field_masked_[z_field] ** 2
        field_mass = 10 ** field_masked_[mass_MED_field]
        field_mass2 = (10 ** field_masked_[mass_MED_field]) ** 2
        field_mass3 = (10 ** field_masked_[mass_MED_field]) ** 3
        field_oneoverr = 1.0/field_masked_[cell_sep]
        field_zoverr = field_zweight * 1.0 / field_masked_[cell_sep]
        field_massoverr = field_mass * 1.0 / field_masked_[cell_sep]
        field_mass2overr = field_mass2 * 1.0 / field_masked_[cell_sep]
        field_mass3overr = field_mass3 * 1.0 / field_masked_[cell_sep]
        field_flexion = field_mass / (field_masked_[cell_sep] ** 3) # ignoring const 4GM/c^2 and unit conversion and D_PS / (D_P * D_S): M / R^3
        field_tidal = field_mass / (field_masked_[cell_sep] ** 2) # (M * (D_PS / (D_P * D_S)) / R^2) * (1 - D_12 * D_S / (D_P * D_LS))
        convergence = np.sqrt(field_mass) / field_masked_[cell_sep] # sqrt(const * M * (D_PS / (D_P * D_S))) / R
        field_masshalo = 10 ** field_masked_[Mhalo_field]
        convergencehalo = np.sqrt(field_masshalo) / field_masked_[cell_sep] # sqrt(const * M * (D_PS / (D_P * D_S))) / R
        if type == "meds":
            field_zweight_[k][l][i][j] = np.median(field_zweight) * w_gal_
            field_mass_[k][l][i][j] = np.median(field_mass) * w_gal_
            field_mass2_[k][l][i][j] = np.median(field_mass2) * w_gal_
            field_mass3_[k][l][i][j] = np.median(field_mass3) * w_gal_
            field_oneoverr_[k][l][i][j] = np.median(field_oneoverr) * w_gal_
            field_zoverr_[k][l][i][j] = np.median(field_zoverr) * w_gal_
            field_massoverr_[k][l][i][j] = np.median(field_massoverr) * w_gal_
            field_mass2overr_[k][l][i][j] = np.median(field_mass2overr) * w_gal_
            field_mass3overr_[k][l][i][j] = np.median(field_mass3overr) * w_gal_
            field_flexion_[k][l][i][j] = np.median(field_flexion) * w_gal_
            field_tidal_[k][l][i][j] = np.median(field_tidal) * w_gal_
            field_convergence_[k][l][i][j] = np.median(convergence) * w_gal_
            field_convergencehalo_[k][l][i][j] = np.median(convergencehalo) * w_gal_
        if type == "sum":
            field_zweight_[k][l][i][j] = np.sum(field_zweight)
            field_mass_[k][l][i][j] = np.sum(field_mass)
            field_mass2_[k][l][i][j] = np.sum(field_mass2)
            field_mass3_[k][l][i][j] = np.sum(field_mass3)
            field_oneoverr_[k][l][i][j] = np.sum(field_oneoverr)
            field_zoverr_[k][l][i][j] = np.sum(field_zoverr)
            field_massoverr_[k][l][i][j] = np.sum(field_massoverr)
            field_mass2overr_[k][l][i][j] = np.sum(field_mass2overr)
            field_mass3overr_[k][l][i][j] = np.sum(field_mass3overr)
            field_flexion_[k][l][i][j] = np.sum(field_flexion)
            field_tidal_[k][l][i][j] = np.sum(field_tidal)
            field_convergence_[k][l][i][j] = np.sum(convergence)
            field_convergencehalo_[k][l][i][j] = np.sum(convergencehalo)
        field_mass2rms_[k][l][i][j] = np.sqrt(field_mass2_[k][l][i][j])
        field_mass3rms_[k][l][i][j] = scipy.special.cbrt(field_mass3_[k][l][i][j])
        field_mass2overrms_[k][l][i][j] = np.sqrt(field_mass2overr_[k][l][i][j])
        field_mass3overrms_[k][l][i][j] = scipy.special.cbrt(field_mass3overr_[k][l][i][j])
    return field_gal_,field_zweight_,field_mass_,field_mass2_,field_mass3_,field_oneoverr_,field_zoverr_,field_massoverr_,field_mass2overr_,field_mass3overr_,field_mass2rms_,field_mass3rms_,field_mass2overrms_,field_mass3overrms_,field_flexion_,field_tidal_,field_convergence_,field_convergencehalo_

start_time = time.time()

if lensID == "B1608":
    z_s = 1.39
    z_l = 0.63
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "HE0435":
    z_s = 1.69
    z_l = 0.455
    brightmag = 17.48
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "HE1104":
    z_s = 2.32
    z_l = 0.73
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "RX1131":
    z_s = 0.66
    z_l = 0.295
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "WFI2033":
    z_s = 1.66
    z_l = 0.66
    brightmag = 16.90
    pixnr = 915
    pixlens = 0.2625 * u.arcsec

#dist = distances.Distance()
#dist.OMEGA_M = 0.25 # change to the cosmology used by the Millennium Simulation
#dist.OMEGA_L = 0.75
#dist.h = 0.73
#D_S = dist.angular_diameter_distance(0,z_s) # in Mpc
#D_LS = dist.angular_diameter_distance(z_l,z_s)

# deginrad = 0.0174532925
const = 9.626*(10**(-20)) # 4GM_sol/c^2 in Mpc
radinsec = 206265  # radian in arcsec

#rootlensmask = "/Users/cerusu/Dropbox/Davis_work/code/CODETOCHECK"
rootlenscat = '%s/%s' % (rootlens,lensID) # where the lens catalogue and its masks are

# Behroozi et al 2010 parameters for z < 1:
M10_ = 12.35
M1a_ = 0.28
Ms00_ = 10.72
Ms0a_ = 0.55
b0_ = 0.44
ba_ = 0.18
d0_ = 0.57
da_ = 0.17
g0_ = 1.56
ga_ = 2.51
# z >= 1:
M10 = 12.27
M1a = -0.84
Ms00 = 11.09
Ms0a = 0.56
b0 = 0.65
ba = 0.31
d0 = 0.55
da = -0.12
g0 = 1.12
ga = -0.53

######################################################################
# WORKING ON THE LENS CATALOGUE
######################################################################

print "Initializing the lens catalogue..."

if lensID == "B1608":
    center_lens = SkyCoord('16:09:13.956 +65:32:28.00', frame='fk5', unit=(u.hourangle, u.deg))
if lensID == "HE0435":
    center_lens = SkyCoord('04:38:14.871 -12:17:14.96', frame='fk5', unit=(u.hourangle, u.deg))
if lensID == "HE1104":
    center_lens = SkyCoord('11:06:33.450 -18:21:24.20', frame='fk5', unit=(u.hourangle, u.deg))
if lensID == "RX1131":
    center_lens = SkyCoord('11:31:51.435 -12:31:58.24', frame='fk5', unit=(u.hourangle, u.deg))
if lensID == "WFI2033":
    center_lens = SkyCoord('20:33:42.080 -47:23:43.00', frame='fk5', unit=(u.hourangle, u.deg))

lensbpz = np.loadtxt('%s/%s%sbpz_nobeta_%s_testothersamples.cat' % (rootlenscat,lensID,irac,det), unpack=True)
lenseazy = np.loadtxt('%s/%s%seazy_nobeta_%s_testothersamples.cat' % (rootlenscat,lensID,irac,det), unpack=True)
#msk_lens = fits.open('%s/testmask45.fits' % (rootlenscat)) # testing
msk_lens = fits.open('%s/msk%sarcsecrad%sarcsecgap.fits' % (rootlenscat,radius,inner))
# defining the columns:
x_lens = 0
y_lens = 1
RA_lens = 2
DEC_lens = 3
i_lens = 4
i_err_lens = 5
classify = 7
z_lens = 8 # first instance
Mlens = 9 # first instance
Mhalo = 10 # first instance
# soon to be defined:
sep = 38 # distance in arcsec between the lens center and the galaxies

def lensprep(lenscat):
    coord_lensinit = SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5')
    sep_lens = coord_lensinit.separation(center_lens).arcsec
    sep_lens[sep_lens < 10] = 10 # limiting the minimum distance to the lens to 10 arcsec; this does not mean I'm removing those objects, because I use the mask for that below
    lenscat = np.c_['0',lenscat,sep_lens.reshape((1, sep_lens.shape[0]))] # inserting as the last column of the catalogue
    lenscat = np.delete(lenscat,np.where(msk_lens[0].data[lenscat[y_lens].astype(int),lenscat[x_lens].astype(int)] != 0),axis=1) # remove the masked objects; I tested that this is the correct order of x and y. x and y in lenscat are the actual coordinates in the natural reading of the .fits file (not the reading of python, which inverts axes)
    #print msk_lens[0].data[400,485] # testing
    lenscat = np.delete(lenscat,np.where(lenscat[classify] < 0),axis=1) # removes all stars from the catalogue
    #print lenscat[RA_lens],lenscat[DEC_lens],SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5').separation(center_lens).arcsec,lenscat[classify]
    if remove != 'no':
        removecat = np.loadtxt('%s/%s.cat' % (rootlenscat,remove))
        for i in range(len(removecat)):
            coord_lensinit = SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5') # redefine this whenever I remove rows
            coordremove = SkyCoord('%s %s' %(removecat[i][0],removecat[i][1]), frame='fk5', unit=(u.deg, u.deg))
            #print coordremove
            #print np.where(coord_lensinit.separation(coordremove).arcsec < 0.5)
            lenscat = np.delete(lenscat,np.where(coord_lensinit.separation(coordremove).arcsec < 0.5),axis=1) # removes selected objects from the catalogue, within a selection radius of 0.5 arcsec
            #print np.shape(lenscat)
        global suffix
        suffix = '_%s' % remove
    #for i in range(np.shape(lensbpz)[1]): # used for testing
        #print lensbpz[:,i][x_lens],lensbpz[:,i][y_lens]
    #print lenscat[RA_lens],lenscat[DEC_lens],SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5').separation(center_lens).arcsec,lenscat[classify]
    lenscat[x_lens] = lenscat[x_lens] - (pixnr - (2/pixlens.value) * radius)/2   # rescale the pixel coordinates in the original catalogue as appropriate for the relevant aperture radius; since the original is fractional, I am not using ((pixnr - (2/pixlens.value) * radius)/2).astype(int)
    lenscat[y_lens] = lenscat[y_lens] - (pixnr - (2/pixlens.value) * radius)/2
    return lenscat

lensbpz = lensprep(lensbpz)
lenseazy = lensprep(lenseazy)
#print lensbpz[RA_lens],lensbpz[DEC_lens],SkyCoord(ra=lensbpz[RA_lens]*u.degree, dec=lensbpz[DEC_lens]*u.degree, frame='fk5').separation(center_lens).arcsec,lensbpz[classify]

# Open the lens and field masks and set up the cell grid
print "Masking the lens catalogue corresponding to each CFHTLenS field cell..."

pixCFHT = 0.187 * u.arcsec
msk = fits.open('%s' % fmask)
worldfield = WCS('%s' % fmask)

if radius == 45:
    overlap = 2
if radius == 60:
    overlap = 3
if radius == 90:
    overlap = 4
if radius == 120:
    overlap = 5

mskfrac = msk_lens[0].data[int(round((pixnr - (2/pixlens.value) * radius)/2)) : int(round(pixnr - (pixnr - (2/pixlens.value) * radius)/2)), int(round((pixnr - (2/pixlens.value) * radius)/2)) : int(round(pixnr - (pixnr - (2/pixlens.value) * radius)/2))]  # will be used later in computing the mask cover fraction
unmaskedlens = np.where(mskfrac == 0) # 2D coordinates of mskfrac (the central subsection of msk_lens)
cells_on_a_side = int((len(msk[0].data[1]) * pixCFHT.value) / (2 * radius)) - 1
print "Cells on a side: ", cells_on_a_side
centerfieldspix = np.zeros((overlap,overlap,2))
unmaskedcell = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
    
# Declare the weighted counts:

lens_gal_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_gal_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_gal_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_gal_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))

lens_zweight_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zweight_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zoverr_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zoverr_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_massoverr_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_massoverr_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overr_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overr_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overr_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overr_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2rms_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2rms_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3rms_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3rms_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overrms_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overrms_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overrms_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overrms_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_flexion_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_flexion_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_tidal_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_tidal_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergence_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergence_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergencehalo_24bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergencehalo_23bpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zweight_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zweight_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zoverr_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zoverr_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_massoverr_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_massoverr_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overr_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overr_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overr_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overr_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2rms_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2rms_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3rms_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3rms_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overrms_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overrms_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overrms_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overrms_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_flexion_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_flexion_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_tidal_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_tidal_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergence_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergence_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergencehalo_24eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergencehalo_23eazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))

print("Initialization of lens catalogue was completed in %0.1f seconds" % (time.time() - start_time))

######################################################################
# WORKING ON THE FIELD CATALOGUE
######################################################################

start_timefield = time.time()
print "Reading ", fieldID, "..."

field = np.loadtxt('%s' % fieldID,usecols = [1,2,4,14,16,20,22],unpack=True)
# 1->0 RA
# 2->1 DEC
# 4->2 z_bpz
# 14->3 i
# 16->4 y
# 20->5 mass_BEST
# 22->6 mass_MED
RA_field = 0
DEC_field = 1
z_field = 2
i_field = 3
y_field = 4
mass_BEST_field = 5
mass_MED_field = 6
# will define shortly:
Mhalo_field = 7
cell_xpix = -3
cell_ypix = -2
cell_sep = -1

field = np.delete(field,np.where(field[z_field] > z_s),axis=1) # eliminate objects at redshifts higher than the source
#print 'shape', np.shape(field)
field = np.delete(field,np.where((field[z_field] >= zinf) & (field[z_field] <= zsup)),axis=1) # eliminate objects corresponding to the redshift slice
#print 'shape', np.shape(field)
field[i_field][field[i_field] < 0] = field[y_field][field[i_field] < 0] # for objects with y mags, use those
field = np.delete(field,np.where(field[i_field] < brightmag),axis = 1) # eliminate objects brighter than the upper brightness limit
field[mass_BEST_field][field[mass_BEST_field] < 0] = 9.0 # fix the very rare unphysical masses
field[mass_MED_field][field[mass_MED_field] < 0] = field[mass_BEST_field][field[mass_MED_field] < 0] # fix the more common unphysical med masses with best masses

# compute halo masses

field = np.c_['0',field,np.zeros(field.shape[1]).reshape(1,np.zeros(field.shape[1]).shape[0])] # add one column for halo mass
a = 1 / (1 + field[z_field][field[z_field] <= 1])
logM1a = M10_ + M1a_ * (a - 1)
logMs0a = Ms00_ + Ms0a_ * (a-1)
notlogMs0a = 10 ** logMs0a
b = b0_ + ba_ * (a-1)
d = d0_ + da_ * (a-1)
g = g0_ + ga_ * (a-1)
field[Mhalo_field][field[z_field] <= 1] = logM1a + b * (field[mass_MED_field][field[z_field] <= 1] - logMs0a) + ((10 ** field[mass_MED_field][field[z_field] <= 1]/notlogMs0a)**d)/(1+(10 ** field[mass_MED_field][field[z_field] <= 1]/notlogMs0a)**(-g)) - 1/2
a = 1 / (1 + field[z_field][field[z_field] > 1])
del logM1a
del logMs0a
del notlogMs0a
del b
del d
del g
logM1a = M10 + M1a * (a-1)
logMs0a = Ms00 + Ms0a * (a-1)
notlogMs0a = 10 ** logMs0a
b = b0 + ba * (a-1)
d = d0 + da * (a-1)
g = g0 + ga * (a-1)
field[Mhalo_field][field[z_field] > 1] = logM1a + b * (field[mass_MED_field][field[z_field] > 1] - logMs0a) + ((10 ** field[mass_MED_field][field[z_field] > 1]/notlogMs0a)**d)/(1+(10 ** field[mass_MED_field][field[z_field] > 1]/notlogMs0a)**(-g)) - 1/2

#field[Mhalo_field][field[Mhalo_field] - field[mass_MED_field] > 4] = field[mass_MED_field][field[Mhalo_field] - field[mass_MED_field] > 4] + 4 # halo masses can get very large, resulting in negative convergence

fieldpix = np.zeros((3,field.shape[1])) # this will contain the fieldpix[1], fieldpix[0] and the separation in arcsec from the cell center
field = np.c_['0',field,fieldpix.reshape(3,fieldpix.shape[1])] # adding fieldpix as extension to the catalogue
fieldpix = worldfield.wcs_world2pix(field[RA_field], field[DEC_field], 1) # find pixel coordinates of each field galaxy, relative to the fmask .fits file (therefore with CFHTLenS pixel size); first column is the physical image x axis
field[cell_xpix] = fieldpix[1] # the physical image y axis, in agreement with the lens section above
field[cell_ypix] = fieldpix[0]

# declarations

field_gal_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_gal_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zweight_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zweight_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_oneoverr_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_oneoverr_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zoverr_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zoverr_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_massoverr_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_massoverr_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overr_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overr_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overr_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overr_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2rms_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2rms_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3rms_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3rms_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overrms_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overrms_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overrms_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overrms_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_flexion_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_flexion_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_tidal_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_tidal_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergence_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergence_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergencehalo_24 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergencehalo_23 = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

print("Initialization of field catalogue was completed in %0.1f seconds" % (time.time() - start_timefield))

######################################################################
# COMPUTING THE WEIGHTS
######################################################################

start_timeweights = time.time()

for k in range(overlap):
    for l in range(overlap):
        print k,l
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                '''Define the center of each cell as a matrix of pixel centers and choose which cells to discard because of large masked area'''
                
                xlow = (2.0 * radius * (u.arcsec / pixCFHT).value) * i + (2.0 * radius * (u.arcsec / pixCFHT).value) * k / overlap   # coordinates of cell edges in the field mask; here x axis refers to the array axis, so the image y axis
                xhigh = (2.0 * radius * (u.arcsec / pixCFHT).value) * i + (2.0 * radius * (u.arcsec / pixCFHT).value) + (2.0 * radius * (u.arcsec / pixCFHT).value) * k / overlap
                ylow = (2.0 * radius * (u.arcsec / pixCFHT).value) * j + (2.0 * radius * (u.arcsec / pixCFHT).value) * l / overlap
                yhigh = (2.0 * radius * (u.arcsec / pixCFHT).value) * j + (2.0 * radius * (u.arcsec / pixCFHT).value) + (2.0 * radius * (u.arcsec / pixCFHT).value) * l / overlap
                centerfieldspix[k][l][0] = xlow + (xhigh - xlow) / 2.0 # coordinates of cell centers in the field mask, in CFTLENS-size pixels
                centerfieldspix[k][l][1] = ylow + (yhigh - ylow) / 2.0
                #centerfields120[i][j] = SkyCoord(ra=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2.0, ylow + (yhigh-ylow)/2.0, 0)[0]*u.degree, dec=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2.0, ylow + (yhigh-ylow)/2.0,0)[1]*u.degree, frame='fk5')
                unmaskedfieldx = xlow + (1.0 * unmaskedlens[0]/((2/pixlens.value) * radius)) * (xhigh - xlow) # the pixels in the field mask that correspond to each unmasked pixel in the lens mask; this is used only to compute the number of pixels in the lens field which are not masked by either the lens of field masks; both x and unmaskedlens[0] refer to the image y axis
                unmaskedfieldy = ylow + (1.0 * unmaskedlens[1]/((2/pixlens.value) * radius)) * (yhigh - ylow)
                unmaskedfield = msk[0].data[unmaskedfieldx.astype(int),unmaskedfieldy.astype(int)] # even if the pixel scales for field and lens are different, the size of this array is the number of unmasked pixels in the lens mask; int is needed here to read from msk
                unmaskedcell[k][l][i][j] = unmaskedfield[unmaskedfield == 0].size / float(np.pi * ((radius * (u.arcsec/pixlens).value) ** 2)) # specifies the fraction of the surface of the circular aperture inside the cell, that is not covered by any masks
    
                if unmaskedcell[k][l][i][j] >= 0.5: # only work on fields with > 50% of their surface not covered by masks; among those, also work later on > 75%
                    ''' Mask the lens catalogue using the field submask corresponding to each cell'''
        
                    lenscoords = np.copy(lensbpz[x_lens:y_lens + 1]) # copy the pixel coordinates
                    lenscoords_fieldx = float(xlow)+(1.0 * lenscoords[1]/((2/pixlens.value) * radius))*(float(xhigh)-float(xlow)) # project the lens catalogue onto the field mask; this is good, because it matches the formula for unmaskedfieldx
                    lenscoords_fieldy = float(ylow)+(1.0 * lenscoords[0]/((2/pixlens.value) * radius))*(float(yhigh)-float(ylow))
                    lenscoords_field = msk[0].data[lenscoords_fieldx.astype(int),lenscoords_fieldy.astype(int)]
                    lensbpz_masked = np.c_['0',lensbpz,lenscoords_field.reshape(1,lenscoords_field.shape[0])] # check if the lens catalogue objects fall inside field masks; I tested that this expresion gives the correct order of the rows
                    lensbpz_masked = np.delete(lensbpz_masked,np.where(lensbpz_masked[-1] != 0),axis=1) # remove objects inside a field mask
                    lensbpz_masked = lensbpz_masked[:-1] # delete the last column
                    
                    lenscoords = np.copy(lenseazy[x_lens:y_lens + 1]) # copy the pixel coordinates
                    lenscoords_fieldx = float(xlow)+(1.0 * lenscoords[1]/((2/pixlens.value) * radius))*(float(xhigh)-float(xlow))
                    lenscoords_fieldy = float(ylow)+(1.0 * lenscoords[0]/((2/pixlens.value) * radius))*(float(yhigh)-float(ylow))
                    lenscoords_field = msk[0].data[lenscoords_fieldx.astype(int),lenscoords_fieldy.astype(int)]
                    lenseazy_masked = np.c_['0',lenseazy,lenscoords_field.reshape(1,lenscoords_field.shape[0])] # check if the lens catalogue objects fall inside field masks
                    lenseazy_masked = np.delete(lenseazy_masked,np.where(lenseazy_masked[-1] != 0),axis=1) # remove objects inside a field mask
                    lenseazy_masked = lenseazy_masked[:-1] # delete the last column
                
                    for n in range(samples):
                        lensbpz_masked24 = np.copy(lensbpz_masked)
                        lenseazy_masked24 = np.copy(lenseazy_masked)
                        if n == 0:
                            lensbpz_masked24 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[i_lens] > 24) | (lensbpz_masked24[z_lens] > z_s)),axis=1)
                            lensbpz_masked24 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[z_lens] >= zinf) & (lensbpz_masked24[z_lens] <= zsup)),axis=1) # remove the redshift slice
                            lenseazy_masked24 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[i_lens] > 24) | (lenseazy_masked24[z_lens] > z_s)),axis=1)
                            lenseazy_masked24 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[z_lens] >= zinf) & (lenseazy_masked24[z_lens] <= zsup)),axis=1)
                            lensbpz_masked23 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[i_lens] > 23) | (lensbpz_masked24[z_lens] > z_s)),axis=1)
                            lensbpz_masked23 = np.delete(lensbpz_masked23,np.where((lensbpz_masked23[z_lens] >= zinf) & (lensbpz_masked23[z_lens] <= zsup)),axis=1)
                            lenseazy_masked23 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[i_lens] > 23) | (lenseazy_masked24[z_lens] > z_s)),axis=1)
                            lenseazy_masked23 = np.delete(lenseazy_masked23,np.where((lenseazy_masked23[z_lens] >= zinf) & (lenseazy_masked23[z_lens] <= zsup)),axis=1)
                        else:
                            for o in range(lensbpz_masked24.shape[1]):
                                lensbpz_masked24[i_lens][o] = np.random.normal(lensbpz_masked24[i_lens][o], np.max([lensbpz_masked24[i_err_lens][o],0.005]),1)[0]
                            for o in range(lenseazy_masked24.shape[1]):
                                lenseazy_masked24[i_lens][o] = np.random.normal(lenseazy_masked24[i_lens][o], np.max([lenseazy_masked24[i_err_lens][o],0.005]),1)[0]
                            lensbpz_masked24 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[i_lens] > 24) | (lensbpz_masked24[z_lens + n * 3] > z_s)),axis=1)
                            lensbpz_masked24 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[z_lens + n * 3] >= zinf) & (lensbpz_masked24[z_lens + n * 3] <= zsup)),axis=1)  # remove the redshift slice
                            lenseazy_masked24 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[i_lens] > 24) | (lenseazy_masked24[z_lens + n * 3] > z_s)),axis=1)
                            lenseazy_masked24 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[z_lens + n * 3] >= zinf) & (lenseazy_masked24[z_lens + n * 3] <= zsup)),axis=1)
                            lensbpz_masked23 = np.delete(lensbpz_masked24,np.where((lensbpz_masked24[i_lens] > 23) | (lensbpz_masked24[z_lens + n * 3] > z_s)),axis=1)
                            lensbpz_masked23 = np.delete(lensbpz_masked23,np.where((lensbpz_masked23[z_lens + n * 3] >= zinf) & (lensbpz_masked23[z_lens + n * 3] <= zsup)),axis=1)
                            lenseazy_masked23 = np.delete(lenseazy_masked24,np.where((lenseazy_masked24[i_lens] > 23) | (lenseazy_masked24[z_lens + n * 3] > z_s)),axis=1)
                            lenseazy_masked23 = np.delete(lenseazy_masked23,np.where((lenseazy_masked23[z_lens + n * 3] >= zinf) & (lenseazy_masked23[z_lens + n * 3] <= zsup)),axis=1)
                        lens_gal_24bpz,lens_zweight_24bpz,lens_mass_24bpz,lens_mass2_24bpz,lens_mass3_24bpz,lens_oneoverr_24bpz,lens_zoverr_24bpz,lens_massoverr_24bpz,lens_mass2overr_24bpz,lens_mass3overr_24bpz,lens_flexion_24bpz,lens_tidal_24bpz,lens_convergence_24bpz,lens_convergencehalo_24bpz,lens_mass2rms_24bpz,lens_mass3rms_24bpz,lens_mass2overrms_24bpz,lens_mass3overrms_24bpz = lensinit(lensbpz_masked24,lens_gal_24bpz,lens_zweight_24bpz,lens_mass_24bpz,lens_mass2_24bpz,lens_mass3_24bpz,lens_oneoverr_24bpz,lens_zoverr_24bpz,lens_massoverr_24bpz,lens_mass2overr_24bpz,lens_mass3overr_24bpz,lens_flexion_24bpz,lens_tidal_24bpz,lens_convergence_24bpz,lens_convergencehalo_24bpz,lens_mass2rms_24bpz,lens_mass3rms_24bpz,lens_mass2overrms_24bpz,lens_mass3overrms_24bpz)
                        #print "d ",lens_gal_24bpz[k][l][i][j][n]# test to check if the function actually returns the result globally
                        lens_gal_23bpz,lens_zweight_23bpz,lens_mass_23bpz,lens_mass2_23bpz,lens_mass3_23bpz,lens_oneoverr_23bpz,lens_zoverr_23bpz,lens_massoverr_23bpz,lens_mass2overr_23bpz,lens_mass3overr_23bpz,lens_flexion_23bpz,lens_tidal_23bpz,lens_convergence_23bpz,lens_convergencehalo_23bpz,lens_mass2rms_23bpz,lens_mass3rms_23bpz,lens_mass2overrms_23bpz,lens_mass3overrms_23bpz = lensinit(lensbpz_masked23,lens_gal_23bpz,lens_zweight_23bpz,lens_mass_23bpz,lens_mass2_23bpz,lens_mass3_23bpz,lens_oneoverr_23bpz,lens_zoverr_23bpz,lens_massoverr_23bpz,lens_mass2overr_23bpz,lens_mass3overr_23bpz,lens_flexion_23bpz,lens_tidal_23bpz,lens_convergence_23bpz,lens_convergencehalo_23bpz,lens_mass2rms_23bpz,lens_mass3rms_23bpz,lens_mass2overrms_23bpz,lens_mass3overrms_23bpz)
                        #print "d ",lens_gal_23bpz[k][l][i][j][n]# test to check if the function actually returns the result globally
                        lens_gal_24eazy,lens_zweight_24eazy,lens_mass_24eazy,lens_mass2_24eazy,lens_mass3_24eazy,lens_oneoverr_24eazy,lens_zoverr_24eazy,lens_massoverr_24eazy,lens_mass2overr_24eazy,lens_mass3overr_24eazy,lens_flexion_24eazy,lens_tidal_24eazy,lens_convergence_24eazy,lens_convergencehalo_24eazy,lens_mass2rms_24eazy,lens_mass3rms_24eazy,lens_mass2overrms_24eazy,lens_mass3overrms_24eazy = lensinit(lenseazy_masked24,lens_gal_24eazy,lens_zweight_24eazy,lens_mass_24eazy,lens_mass2_24eazy,lens_mass3_24eazy,lens_oneoverr_24eazy,lens_zoverr_24eazy,lens_massoverr_24eazy,lens_mass2overr_24eazy,lens_mass3overr_24eazy,lens_flexion_24eazy,lens_tidal_24eazy,lens_convergence_24eazy,lens_convergencehalo_24eazy,lens_mass2rms_24eazy,lens_mass3rms_24eazy,lens_mass2overrms_24eazy,lens_mass3overrms_24eazy)
                        #print "d ",lens_gal_24eazy[k][l][i][j][n]# test to check if the function actually returns the result globally
                        lens_gal_23eazy,lens_zweight_23eazy,lens_mass_23eazy,lens_mass2_23eazy,lens_mass3_23eazy,lens_oneoverr_23eazy,lens_zoverr_23eazy,lens_massoverr_23eazy,lens_mass2overr_23eazy,lens_mass3overr_23eazy,lens_flexion_23eazy,lens_tidal_23eazy,lens_convergence_23eazy,lens_convergencehalo_23eazy,lens_mass2rms_23eazy,lens_mass3rms_23eazy,lens_mass2overrms_23eazy,lens_mass3overrms_23eazy = lensinit(lenseazy_masked23,lens_gal_23eazy,lens_zweight_23eazy,lens_mass_23eazy,lens_mass2_23eazy,lens_mass3_23eazy,lens_oneoverr_23eazy,lens_zoverr_23eazy,lens_massoverr_23eazy,lens_mass2overr_23eazy,lens_mass3overr_23eazy,lens_flexion_23eazy,lens_tidal_23eazy,lens_convergence_23eazy,lens_convergencehalo_23eazy,lens_mass2rms_23eazy,lens_mass3rms_23eazy,lens_mass2overrms_23eazy,lens_mass3overrms_23eazy)
                        #print "d ",lens_gal_23eazy[k][l][i][j][n]# test to check if the function actually returns the result globally

                    '''Compute weights for the field catalogue'''
                    field_masked = np.copy(field)
                    field_masked = np.delete(field_masked,np.where((field_masked[cell_xpix] < xlow) | (field_masked[cell_xpix] > xhigh) | (field_masked[cell_ypix] < ylow) | (field_masked[cell_ypix] > yhigh)),axis=1) # isolate the field objects matching in position the field coordinates on which the lens catalogue is projected
                    #field_masked = np.delete(field_masked,np.where((field_masked[cell_xpix] > xlow))
                    field_masked[cell_xpix] = (field_masked[cell_xpix] - (centerfieldspix[k][l][0] - radius * (u.arcsec / pixCFHT).value)) * (pixCFHT / pixlens).value # the last multiplication brings the pixels to the lens pixel size; work on the physical image y axis
                    field_masked[cell_ypix] = (field_masked[cell_ypix] - (centerfieldspix[k][l][1] - radius * (u.arcsec / pixCFHT).value)) * (pixCFHT / pixlens).value
                    field_masked[cell_sep] = np.sqrt((field_masked[cell_xpix] - (radius * 2/pixlens.value)/2) ** 2 + (field_masked[cell_ypix] - (radius * 2/pixlens.value)/2) ** 2) * (pixlens / u.arcsec).value # separation in arcsec from the center of the cell
                    field_masked[cell_sep][field_masked[cell_sep] < 10] = 10
                    field_masked_24 = np.delete(field_masked,np.where(msk_lens[0].data[np.array([field_masked[cell_xpix] + (pixnr - (2/pixlens.value) * radius)/2]).astype(int),np.array([field_masked[cell_ypix] + (pixnr - (2/pixlens.value) * radius)/2]).astype(int)] != 0),axis=1) # remove field objects falling inside the lens mask; needed to account for the fact that the mask is always pixnr pixels on a side
                    #if np.max(field_masked_24[cell_sep]) > radius: print "a", np.max(field_masked_24[cell_sep]) # testing; due to the pixelated mask, there are objects very close to the radius limit but just very slightly (subarcsec) away; I will not do anything about this
                    field_masked_23 = np.delete(field_masked_24,np.where(field_masked_24[i_field] > 23),axis=1)
                    w_gal_24 = np.shape(field_masked_24)[1]
                    w_gal_23 = np.shape(field_masked_23)[1]
                    #mmm = np.copy(field_convergencehalo_23)
                    field_gal_24,field_zweight_24,field_mass_24,field_mass2_24,field_mass3_24,field_oneoverr_24,field_zoverr_24,field_massoverr_24,field_mass2overr_24,field_mass3overr_24,field_mass2rms_24,field_mass3rms_24,field_mass2overrms_24,field_mass3overrms_24,field_flexion_24,field_tidal_24,field_convergence_24,field_convergencehalo_24 = fieldinit(field_masked_24,w_gal_24,field_gal_24,field_zweight_24,field_mass_24,field_mass2_24,field_mass3_24,field_oneoverr_24,field_zoverr_24,field_massoverr_24,field_mass2overr_24,field_mass3overr_24,field_mass2rms_24,field_mass3rms_24,field_mass2overrms_24,field_mass3overrms_24,field_flexion_24,field_tidal_24,field_convergence_24,field_convergencehalo_24)
                    field_gal_23,field_zweight_23,field_mass_23,field_mass2_23,field_mass3_23,field_oneoverr_23,field_zoverr_23,field_massoverr_23,field_mass2overr_23,field_mass3overr_23,field_mass2rms_23,field_mass3rms_23,field_mass2overrms_23,field_mass3overrms_23,field_flexion_23,field_tidal_23,field_convergence_23,field_convergencehalo_23 = fieldinit(field_masked_23,w_gal_23,field_gal_23,field_zweight_23,field_mass_23,field_mass2_23,field_mass3_23,field_oneoverr_23,field_zoverr_23,field_massoverr_23,field_mass2overr_23,field_mass3overr_23,field_mass2rms_23,field_mass3rms_23,field_mass2overrms_23,field_mass3overrms_23,field_flexion_23,field_tidal_23,field_convergence_23,field_convergencehalo_23)
                    #print np.min(mmm - field_convergencehalo_23)

                    ''' TEST THAT BOTH THE LENS AND FIELD CATALOGUE ARE PROPERLY MASKED AND MATCHED AGAINST EACHOTHER (BROUGHT THE LENS CATALOGUE TO FIELD PIXEL SCALE)'''
                    #if (k == 0) and (l == 0) and (i == 18) and (j == 15): # tested for radius 45, not appropriate for 120
                        #print "TESTING %d %d %d %d" %(k,l,i,j)
                        #xxx = lenscoords_fieldx - xlow # corresponds to the physical image y axis
                        #yyy = lenscoords_fieldy - ylow
                        #xy_wcs = SkyCoord(ra=worldfield.wcs_pix2world(lenscoords_fieldy, lenscoords_fieldx, 1)[0]*u.degree, dec=worldfield.wcs_pix2world(lenscoords_fieldy, lenscoords_fieldx, 1)[1]*u.degree, frame='fk5') # I checked that the order of the axes is correct, if lenscoords_fieldy represents x axis in the physical image
                        #ra_wcs = np.copy(xxx)
                        #dec_wcs = np.copy(yyy)
                        #for ii in range(len(xxx)):
                            #ra_wcs[ii] = xy_wcs[ii].ra.value
                            #dec_wcs[ii] = xy_wcs[ii].dec.value
                        #np.savetxt("/Users/perseus/Desktop/CFHTLenSmasks/test_lens_all.cat",np.c_[yyy,xxx]) # so that it corresponds to the natural image axes
                        #np.savetxt("/Users/perseus/Desktop/CFHTLenSmasks/test_lens_unmasked.cat",np.c_[yyy[lenscoords_field == 0],xxx[lenscoords_field == 0]])
                        #np.savetxt("/Users/perseus/Desktop/CFHTLenSmasks/test_lens_unmasked_wcs.reg",np.c_[ra_wcs[lenscoords_field == 0],dec_wcs[lenscoords_field == 0]],fmt='j2000; circle %s %s 2\"')
                        #print ylow,yhigh,xlow,xhigh,centerfieldspix[k][l]
                        #data = msk[0].data[xlow:xhigh,ylow:yhigh]
                        #newf = fits.PrimaryHDU()
                        #newf.data = msk[0].data[xlow:xhigh,ylow:yhigh]
                        #newf.header = msk[0].header
                        #newf.header.update(worldfield[xlow:xhigh,ylow:yhigh].to_header())
                        #fits.writeto("/Users/perseus/Desktop/CFHTLenSmasks/test_lens_fieldmask.fits", newf.data, newf.header, clobber=True)
                        #data, header = fits.getdata('%s/msk%s_asecrad%s_no%sarcsec.fits' % (rootlenscat,lensID,radius,inner), header=True)
                        ##data, header = fits.getdata('%s/testmask45.fits' % rootlenscat, header=True)
                        #data = data[(pixnr - (2/pixlens.value) * radius)/2 : pixnr - (pixnr - (2/pixlens.value) * radius)/2,(pixnr - (2/pixlens.value) * radius)/2 : pixnr - (pixnr - (2/pixlens.value) * radius)/2]
                        #data = scipy.ndimage.interpolation.zoom(data,0.2626/0.187)
                        #fits.writeto("/Users/perseus/Desktop/CFHTLenSmasks/test_lens_lensmaskrescaled.fits", data, header, clobber=True)
                        #field_masked_24[RA_field],field_masked_24[DEC_field]
                        #np.savetxt("/Users/perseus/Desktop/CFHTLenSmasks/test_field_wcs.reg",np.c_[field_masked_24[RA_field],field_masked_24[DEC_field]],fmt='j2000; circle %s %s 2\"')

msk_lens.close()
msk.close()

print("Finishing lens and field catalogues in %0.1f seconds" % (time.time() - start_timeweights))

######################################################################
# WRITING OUTPUT
######################################################################

start_write = time.time()

print "Writing output..."

count = ">75 percent and >50 percent: %d %d, %d %d" % (unmaskedcell[unmaskedcell>=0.75].shape[0], unmaskedcell.shape[0] * unmaskedcell.shape[1] * unmaskedcell.shape[2] * unmaskedcell.shape[3], unmaskedcell[unmaskedcell>=0.5].shape[0], unmaskedcell.shape[0] * unmaskedcell.shape[1] * unmaskedcell.shape[2] * unmaskedcell.shape[3])
print count
mskname = 'msk%sarcsecrad%sarcsecgap.fits'[0:-5] % (radius,inner)
fcount = open('%s/%s_%s_%s_%s_%s_%s_zgap%s_%s%s_count_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'w') # [-25:-4] corresponds to strings of the form W1m0m0_24galphotmstar
fcount.write(count)
fcount.close()

os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_0_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix)) # '-f' ignores non-existent files
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_0_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_1_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_1_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_2_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_2_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_3_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_3_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_4_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_4_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_5_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_5_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_6_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_6_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_7_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_7_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_8_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_8_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_9_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))
os.system('rm -f %s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_9_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix))

file50_0 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_0_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_0 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_0_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_1 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_1_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_1 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_1_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_2 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_2_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_2 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_2_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_3 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_3_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_3 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_3_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_4 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_4_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_4 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_4_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_5 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_5_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_5 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_5_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_6 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_6_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_6 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_6_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_7 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_7_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_7 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_7_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_8 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_8_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_8 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_8_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file50_9 = open('%s/%s_50_%s_%s_%s_%s_%s_zgap%s_%s%s_9_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')
file75_9 = open('%s/%s_75_%s_%s_%s_%s_%s_zgap%s_%s%s_9_testothersamples.lst' % (output,fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),'a')

header = "#      1_overlap_x  2_overlap_y  3_cell_x  4_cell_y  5_lens_gal_24bpz  6_lens_gal_23bpz  7_lens_zweight_24bpz  8_lens_zweight_23bpz 9_lens_mass_24bpz 10_lens_mass_23bpz 11_lens_mass2_24bpz 12_lens_mass2_23bpz 13_lens_mass3_24bpz 14_lens_mass3_23bpz 15_lens_oneoverr_24bpz 16_lens_oneoverr_23bpz 17_lens_zoverr_24bpz 18_lens_zoverr_23bpz 19_lens_massoverr_24bpz 20_lens_massoverr_23bpz 21_lens_mass2overr_24bpz 22_lens_mass2overr_23bpz 23_lens_mass3overr_24bpz 24_lens_mass3overr_23bpz 25_lens_mass2rms_24bpz 26_lens_mass2rms_23bpz 27_lens_mass3rms_24bpz 28_lens_mass3rms_23bpz 29_lens_mass2overrms_24bpz 30_lens_mass2overrms_23bpz 31_lens_mass3overrms_24bpz 32_lens_mass3overrms_23bpz 33_lens_flexion_24bpz 34_lens_flexion_23bpz 35_lens_tidal_24bpz 36_lens_tidal_23bpz 37_lens_convergence_24bpz 38_lens_convergence_23bpz 39_lens_convergencehalo_24bpz 40_lens_convergencehalo_23bpz 41_lens_gal_24eazy 42_lens_gal_23eazy 43_lens_zweight_24eazy 44_lens_zweight_23eazy 45_lens_mass_24eazy 46_lens_mass_23eazy 47_lens_mass2_24eazy 48_lens_mass2_23eazy 49_lens_mass3_24eazy 50_lens_mass3_23eazy 51_lens_oneoverr_24eazy 52_lens_oneoverr_23eazy 53_lens_zoverr_24eazy 54_lens_zoverr_23eazy  55_lens_massoverr_24eazy  56_lens_massoverr_23eazy  57_lens_mass2overr_24eazy  58_lens_mass2overr_23eazy  59_lens_mass3overr_24eazy  60_lens_mass3overr_23eazy  61_lens_mass2rms_24eazy  62_lens_mass2rms_23eazy  63_lens_mass3rms_24eazy  64_lens_mass3rms_23eazy  65_lens_mass2overrms_24eazy  66_lens_mass2overrms_23eazy  67_lens_mass3overrms_24eazy  68_lens_mass3overrms_23eazy  69_lens_flexion_24eazy  70_lens_flexion_23eazy  71_lens_tidal_24eazy  72_lens_tidal_23eazy  73_lens_convergence_24eazy  74_lens_convergence_23eazy  75_lens_convergencehalo_24eazy  76_lens_convergencehalo_23eazy \n"

file50_0.write(header)
file75_0.write(header)
file50_1.write(header)
file75_1.write(header)
file50_2.write(header)
file75_2.write(header)
file50_3.write(header)
file75_3.write(header)
file50_4.write(header)
file75_4.write(header)
file50_5.write(header)
file75_5.write(header)
file50_6.write(header)
file75_6.write(header)
file50_7.write(header)
file75_7.write(header)
file50_8.write(header)
file75_8.write(header)
file50_9.write(header)
file75_9.write(header)

for k in range(overlap):
    for l in range(overlap):
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                for n in range(10):
                    #print "d",i,j,k,l,field_gal_24[k][l][i][j] # test to check if the function actually returns the result globally
                    #print "d",i,j,k,l,field_gal_23[k][l][i][j] # test to check if the function actually returns the result globally
                    if (unmaskedcell[k][l][i][j] >= 0.50) & (np.min(lens_gal_24bpz[k][l][i][j]) != 0) & (np.min(lens_gal_23bpz[k][l][i][j]) != 0) & (np.min(lens_gal_24eazy[k][l][i][j]) != 0) & (np.min(lens_gal_23eazy[k][l][i][j]) != 0) & (field_gal_24[k][l][i][j] != 0) & (field_gal_23[k][l][i][j] != 0):
                        string='%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n'%(k,l,i,j,1.0*lens_gal_24bpz[k][l][i][j][n]/field_gal_24[k][l][i][j],1.0*lens_gal_23bpz[k][l][i][j][n]/field_gal_23[k][l][i][j],1.0*lens_zweight_24bpz[k][l][i][j][n]/field_zweight_24[k][l][i][j],1.0*lens_zweight_23bpz[k][l][i][j][n]/field_zweight_23[k][l][i][j],1.0*lens_mass_24bpz[k][l][i][j][n]/field_mass_24[k][l][i][j],1.0*lens_mass_23bpz[k][l][i][j][n]/field_mass_23[k][l][i][j],1.0*lens_mass2_24bpz[k][l][i][j][n]/field_mass2_24[k][l][i][j],1.0*lens_mass2_23bpz[k][l][i][j][n]/field_mass2_23[k][l][i][j],1.0*lens_mass3_24bpz[k][l][i][j][n]/field_mass3_24[k][l][i][j],1.0*lens_mass3_23bpz[k][l][i][j][n]/field_mass3_23[k][l][i][j],1.0*lens_oneoverr_24bpz[k][l][i][j][n]/field_oneoverr_24[k][l][i][j],1.0*lens_oneoverr_23bpz[k][l][i][j][n]/field_oneoverr_23[k][l][i][j],1.0*lens_zoverr_24bpz[k][l][i][j][n]/field_zoverr_24[k][l][i][j],1.0*lens_zoverr_23bpz[k][l][i][j][n]/field_zoverr_23[k][l][i][j],1.0*lens_massoverr_24bpz[k][l][i][j][n]/field_massoverr_24[k][l][i][j],1.0*lens_massoverr_23bpz[k][l][i][j][n]/field_massoverr_23[k][l][i][j],1.0*lens_mass2overr_24bpz[k][l][i][j][n]/field_mass2overr_24[k][l][i][j],1.0*lens_mass2overr_23bpz[k][l][i][j][n]/field_mass2overr_23[k][l][i][j],1.0*lens_mass3overr_24bpz[k][l][i][j][n]/field_mass3overr_24[k][l][i][j],1.0*lens_mass3overr_23bpz[k][l][i][j][n]/field_mass3overr_23[k][l][i][j],1.0*lens_mass2rms_24bpz[k][l][i][j][n]/field_mass2rms_24[k][l][i][j],1.0*lens_mass2rms_23bpz[k][l][i][j][n]/field_mass2rms_23[k][l][i][j],1.0*lens_mass3rms_24bpz[k][l][i][j][n]/field_mass3rms_24[k][l][i][j],1.0*lens_mass3rms_23bpz[k][l][i][j][n]/field_mass3rms_23[k][l][i][j],1.0*lens_mass2overrms_24bpz[k][l][i][j][n]/field_mass2overrms_24[k][l][i][j],1.0*lens_mass2overrms_23bpz[k][l][i][j][n]/field_mass2overrms_23[k][l][i][j],1.0*lens_mass3overrms_24bpz[k][l][i][j][n]/field_mass3overrms_24[k][l][i][j],1.0*lens_mass3overrms_23bpz[k][l][i][j][n]/field_mass3overrms_23[k][l][i][j],1.0*lens_flexion_24bpz[k][l][i][j][n]/field_flexion_24[k][l][i][j],1.0*lens_flexion_23bpz[k][l][i][j][n]/field_flexion_23[k][l][i][j],1.0*lens_tidal_24bpz[k][l][i][j][n]/field_tidal_24[k][l][i][j],1.0*lens_tidal_23bpz[k][l][i][j][n]/field_tidal_23[k][l][i][j],1.0*lens_convergence_24bpz[k][l][i][j][n]/field_convergence_24[k][l][i][j],1.0*lens_convergence_23bpz[k][l][i][j][n]/field_convergence_23[k][l][i][j],1.0*lens_convergencehalo_24bpz[k][l][i][j][n]/field_convergencehalo_24[k][l][i][j],1.0*lens_convergencehalo_23bpz[k][l][i][j][n]/field_convergencehalo_23[k][l][i][j],1.0*lens_gal_24eazy[k][l][i][j][n]/field_gal_24[k][l][i][j],1.0*lens_gal_23eazy[k][l][i][j][n]/field_gal_23[k][l][i][j],1.0*lens_zweight_24eazy[k][l][i][j][n]/field_zweight_24[k][l][i][j],1.0*lens_zweight_23eazy[k][l][i][j][n]/field_zweight_23[k][l][i][j],1.0*lens_mass_24eazy[k][l][i][j][n]/field_mass_24[k][l][i][j],1.0*lens_mass_23eazy[k][l][i][j][n]/field_mass_23[k][l][i][j],1.0*lens_mass2_24eazy[k][l][i][j][n]/field_mass2_24[k][l][i][j],1.0*lens_mass2_23eazy[k][l][i][j][n]/field_mass2_23[k][l][i][j],1.0*lens_mass3_24eazy[k][l][i][j][n]/field_mass3_24[k][l][i][j],1.0*lens_mass3_23eazy[k][l][i][j][n]/field_mass3_23[k][l][i][j],1.0*lens_oneoverr_24eazy[k][l][i][j][n]/field_oneoverr_24[k][l][i][j],1.0*lens_oneoverr_23eazy[k][l][i][j][n]/field_oneoverr_23[k][l][i][j],1.0*lens_zoverr_24eazy[k][l][i][j][n]/field_zoverr_24[k][l][i][j],1.0*lens_zoverr_23eazy[k][l][i][j][n]/field_zoverr_23[k][l][i][j],1.0*lens_massoverr_24eazy[k][l][i][j][n]/field_massoverr_24[k][l][i][j],1.0*lens_massoverr_23eazy[k][l][i][j][n]/field_massoverr_23[k][l][i][j],1.0*lens_mass2overr_24eazy[k][l][i][j][n]/field_mass2overr_24[k][l][i][j],1.0*lens_mass2overr_23eazy[k][l][i][j][n]/field_mass2overr_23[k][l][i][j],1.0*lens_mass3overr_24eazy[k][l][i][j][n]/field_mass3overr_24[k][l][i][j],1.0*lens_mass3overr_23eazy[k][l][i][j][n]/field_mass3overr_23[k][l][i][j],1.0*lens_mass2rms_24eazy[k][l][i][j][n]/field_mass2rms_24[k][l][i][j],1.0*lens_mass2rms_23eazy[k][l][i][j][n]/field_mass2rms_23[k][l][i][j],1.0*lens_mass3rms_24eazy[k][l][i][j][n]/field_mass3rms_24[k][l][i][j],1.0*lens_mass3rms_23eazy[k][l][i][j][n]/field_mass3rms_23[k][l][i][j],1.0*lens_mass2overrms_24eazy[k][l][i][j][n]/field_mass2overrms_24[k][l][i][j],1.0*lens_mass2overrms_23eazy[k][l][i][j][n]/field_mass2overrms_23[k][l][i][j],1.0*lens_mass3overrms_24eazy[k][l][i][j][n]/field_mass3overrms_24[k][l][i][j],1.0*lens_mass3overrms_23eazy[k][l][i][j][n]/field_mass3overrms_23[k][l][i][j],1.0*lens_flexion_24eazy[k][l][i][j][n]/field_flexion_24[k][l][i][j],1.0*lens_flexion_23eazy[k][l][i][j][n]/field_flexion_23[k][l][i][j],1.0*lens_tidal_24eazy[k][l][i][j][n]/field_tidal_24[k][l][i][j],1.0*lens_tidal_23eazy[k][l][i][j][n]/field_tidal_23[k][l][i][j],1.0*lens_convergence_24eazy[k][l][i][j][n]/field_convergence_24[k][l][i][j],1.0*lens_convergence_23eazy[k][l][i][j][n]/field_convergence_23[k][l][i][j],1.0*lens_convergencehalo_24eazy[k][l][i][j][n]/field_convergencehalo_24[k][l][i][j],1.0*lens_convergencehalo_23eazy[k][l][i][j][n]/field_convergencehalo_23[k][l][i][j])
                        if n == 0:
                            file50_0.write(string)
                        if n == 1:
                            file50_1.write(string)
                        if n == 2:
                            file50_2.write(string)
                        if n == 3:
                            file50_3.write(string)
                        if n == 4:
                            file50_4.write(string)
                        if n == 5:
                            file50_5.write(string)
                        if n == 6:
                            file50_6.write(string)
                        if n == 7:
                            file50_7.write(string)
                        if n == 8:
                            file50_8.write(string)
                        if n == 9:
                            file50_9.write(string)
for k in range(overlap):
    for l in range(overlap):
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                for n in range(10):
                    if (unmaskedcell[k][l][i][j] >= 0.75) & (np.min(lens_gal_24bpz[k][l][i][j]) != 0) & (np.min(lens_gal_23bpz[k][l][i][j]) != 0) & (np.min(lens_gal_24eazy[k][l][i][j]) != 0) & (np.min(lens_gal_23eazy[k][l][i][j]) != 0) & (field_gal_24[k][l][i][j] != 0) & (field_gal_23[k][l][i][j] != 0):
                        string='%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n'%(k,l,i,j,1.0*lens_gal_24bpz[k][l][i][j][n]/field_gal_24[k][l][i][j],1.0*lens_gal_23bpz[k][l][i][j][n]/field_gal_23[k][l][i][j],1.0*lens_zweight_24bpz[k][l][i][j][n]/field_zweight_24[k][l][i][j],1.0*lens_zweight_23bpz[k][l][i][j][n]/field_zweight_23[k][l][i][j],1.0*lens_mass_24bpz[k][l][i][j][n]/field_mass_24[k][l][i][j],1.0*lens_mass_23bpz[k][l][i][j][n]/field_mass_23[k][l][i][j],1.0*lens_mass2_24bpz[k][l][i][j][n]/field_mass2_24[k][l][i][j],1.0*lens_mass2_23bpz[k][l][i][j][n]/field_mass2_23[k][l][i][j],1.0*lens_mass3_24bpz[k][l][i][j][n]/field_mass3_24[k][l][i][j],1.0*lens_mass3_23bpz[k][l][i][j][n]/field_mass3_23[k][l][i][j],1.0*lens_oneoverr_24bpz[k][l][i][j][n]/field_oneoverr_24[k][l][i][j],1.0*lens_oneoverr_23bpz[k][l][i][j][n]/field_oneoverr_23[k][l][i][j],1.0*lens_zoverr_24bpz[k][l][i][j][n]/field_zoverr_24[k][l][i][j],1.0*lens_zoverr_23bpz[k][l][i][j][n]/field_zoverr_23[k][l][i][j],1.0*lens_massoverr_24bpz[k][l][i][j][n]/field_massoverr_24[k][l][i][j],1.0*lens_massoverr_23bpz[k][l][i][j][n]/field_massoverr_23[k][l][i][j],1.0*lens_mass2overr_24bpz[k][l][i][j][n]/field_mass2overr_24[k][l][i][j],1.0*lens_mass2overr_23bpz[k][l][i][j][n]/field_mass2overr_23[k][l][i][j],1.0*lens_mass3overr_24bpz[k][l][i][j][n]/field_mass3overr_24[k][l][i][j],1.0*lens_mass3overr_23bpz[k][l][i][j][n]/field_mass3overr_23[k][l][i][j],1.0*lens_mass2rms_24bpz[k][l][i][j][n]/field_mass2rms_24[k][l][i][j],1.0*lens_mass2rms_23bpz[k][l][i][j][n]/field_mass2rms_23[k][l][i][j],1.0*lens_mass3rms_24bpz[k][l][i][j][n]/field_mass3rms_24[k][l][i][j],1.0*lens_mass3rms_23bpz[k][l][i][j][n]/field_mass3rms_23[k][l][i][j],1.0*lens_mass2overrms_24bpz[k][l][i][j][n]/field_mass2overrms_24[k][l][i][j],1.0*lens_mass2overrms_23bpz[k][l][i][j][n]/field_mass2overrms_23[k][l][i][j],1.0*lens_mass3overrms_24bpz[k][l][i][j][n]/field_mass3overrms_24[k][l][i][j],1.0*lens_mass3overrms_23bpz[k][l][i][j][n]/field_mass3overrms_23[k][l][i][j],1.0*lens_flexion_24bpz[k][l][i][j][n]/field_flexion_24[k][l][i][j],1.0*lens_flexion_23bpz[k][l][i][j][n]/field_flexion_23[k][l][i][j],1.0*lens_tidal_24bpz[k][l][i][j][n]/field_tidal_24[k][l][i][j],1.0*lens_tidal_23bpz[k][l][i][j][n]/field_tidal_23[k][l][i][j],1.0*lens_convergence_24bpz[k][l][i][j][n]/field_convergence_24[k][l][i][j],1.0*lens_convergence_23bpz[k][l][i][j][n]/field_convergence_23[k][l][i][j],1.0*lens_convergencehalo_24bpz[k][l][i][j][n]/field_convergencehalo_24[k][l][i][j],1.0*lens_convergencehalo_23bpz[k][l][i][j][n]/field_convergencehalo_23[k][l][i][j],1.0*lens_gal_24eazy[k][l][i][j][n]/field_gal_24[k][l][i][j],1.0*lens_gal_23eazy[k][l][i][j][n]/field_gal_23[k][l][i][j],1.0*lens_zweight_24eazy[k][l][i][j][n]/field_zweight_24[k][l][i][j],1.0*lens_zweight_23eazy[k][l][i][j][n]/field_zweight_23[k][l][i][j],1.0*lens_mass_24eazy[k][l][i][j][n]/field_mass_24[k][l][i][j],1.0*lens_mass_23eazy[k][l][i][j][n]/field_mass_23[k][l][i][j],1.0*lens_mass2_24eazy[k][l][i][j][n]/field_mass2_24[k][l][i][j],1.0*lens_mass2_23eazy[k][l][i][j][n]/field_mass2_23[k][l][i][j],1.0*lens_mass3_24eazy[k][l][i][j][n]/field_mass3_24[k][l][i][j],1.0*lens_mass3_23eazy[k][l][i][j][n]/field_mass3_23[k][l][i][j],1.0*lens_oneoverr_24eazy[k][l][i][j][n]/field_oneoverr_24[k][l][i][j],1.0*lens_oneoverr_23eazy[k][l][i][j][n]/field_oneoverr_23[k][l][i][j],1.0*lens_zoverr_24eazy[k][l][i][j][n]/field_zoverr_24[k][l][i][j],1.0*lens_zoverr_23eazy[k][l][i][j][n]/field_zoverr_23[k][l][i][j],1.0*lens_massoverr_24eazy[k][l][i][j][n]/field_massoverr_24[k][l][i][j],1.0*lens_massoverr_23eazy[k][l][i][j][n]/field_massoverr_23[k][l][i][j],1.0*lens_mass2overr_24eazy[k][l][i][j][n]/field_mass2overr_24[k][l][i][j],1.0*lens_mass2overr_23eazy[k][l][i][j][n]/field_mass2overr_23[k][l][i][j],1.0*lens_mass3overr_24eazy[k][l][i][j][n]/field_mass3overr_24[k][l][i][j],1.0*lens_mass3overr_23eazy[k][l][i][j][n]/field_mass3overr_23[k][l][i][j],1.0*lens_mass2rms_24eazy[k][l][i][j][n]/field_mass2rms_24[k][l][i][j],1.0*lens_mass2rms_23eazy[k][l][i][j][n]/field_mass2rms_23[k][l][i][j],1.0*lens_mass3rms_24eazy[k][l][i][j][n]/field_mass3rms_24[k][l][i][j],1.0*lens_mass3rms_23eazy[k][l][i][j][n]/field_mass3rms_23[k][l][i][j],1.0*lens_mass2overrms_24eazy[k][l][i][j][n]/field_mass2overrms_24[k][l][i][j],1.0*lens_mass2overrms_23eazy[k][l][i][j][n]/field_mass2overrms_23[k][l][i][j],1.0*lens_mass3overrms_24eazy[k][l][i][j][n]/field_mass3overrms_24[k][l][i][j],1.0*lens_mass3overrms_23eazy[k][l][i][j][n]/field_mass3overrms_23[k][l][i][j],1.0*lens_flexion_24eazy[k][l][i][j][n]/field_flexion_24[k][l][i][j],1.0*lens_flexion_23eazy[k][l][i][j][n]/field_flexion_23[k][l][i][j],1.0*lens_tidal_24eazy[k][l][i][j][n]/field_tidal_24[k][l][i][j],1.0*lens_tidal_23eazy[k][l][i][j][n]/field_tidal_23[k][l][i][j],1.0*lens_convergence_24eazy[k][l][i][j][n]/field_convergence_24[k][l][i][j],1.0*lens_convergence_23eazy[k][l][i][j][n]/field_convergence_23[k][l][i][j],1.0*lens_convergencehalo_24eazy[k][l][i][j][n]/field_convergencehalo_24[k][l][i][j],1.0*lens_convergencehalo_23eazy[k][l][i][j][n]/field_convergencehalo_23[k][l][i][j])
                        if n == 0:
                            file75_0.write(string)
                        if n == 1:
                            file75_1.write(string)
                        if n == 2:
                            file75_2.write(string)
                        if n == 3:
                            file75_3.write(string)
                        if n == 4:
                            file75_4.write(string)
                        if n == 5:
                            file75_5.write(string)
                        if n == 6:
                            file75_6.write(string)
                        if n == 7:
                            file75_7.write(string)
                        if n == 8:
                            file75_8.write(string)
                        if n == 9:
                            file75_9.write(string)
file50_0.close()
file75_0.close()
file50_1.close()
file75_1.close()
file50_2.close()
file75_2.close()
file50_3.close()
file75_3.close()
file50_4.close()
file75_4.close()
file50_5.close()
file75_5.close()
file50_6.close()
file75_6.close()
file50_7.close()
file75_7.close()
file50_8.close()
file75_8.close()
file50_9.close()
file75_9.close()

print("Writing output completed in %0.1f seconds" % ((time.time() - start_write)))

print("Total time for %s: --- %0.2f seconds ---" % ('%s_%s_%s_%s_%s_%s_zgap%s_%s%s' % (fieldID[-25:-4],mskname,lensID,det,irac,type,zinf,zsup,suffix),(time.time() - start_time)))

print 'Done!'
