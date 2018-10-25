# CE Rusu, Feb 14 2018
# This code makes use of the CFHTLens *galphotmstar.cat files and the lens photometric+Mstar+photoz catalogue; it computes weighted ratios (lens/field) with proper masking, for various radii, limiting mag, number of samples and classification scheme.
# run as: python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_overlap_sampling_nobeta.py WFI2033 /Volumes/LaCieDavis/CFHTcatalogues/W1m0m0_limmaggalphotmstar.fits /Volumes/LaCieDavis/CFHTLenSmasks/W1m0m0_izrgu_finalmask_mosaic.fits /Users/cerusu/Dropbox/Davis_work/code /Volumes/LaCieSubaru/weightedcounts/WFI2033 45 5 IRAC deti meds removegrouphandpicked -1 -1
# the code produces output for two different limiting mags: limbright and 24, or for just limbright (typically 23) in case limmag < 24
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
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
import astropy.table as table
from astropy.table import Table

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
    limmag = 24
    limbright = 23
    photoz = 'bpzeazy'
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "HE0435":
    z_s = 1.69
    z_l = 0.455
    limmag = 24
    limbright = 23
    brightmag = 17.48
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "HE1104":
    z_s = 2.32
    z_l = 0.73
    limmag = 24
    limbright = 23
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "RX1131":
    z_s = 0.66
    z_l = 0.295
    limmag = 24
    limbright = 23
    pixnr = 1200
    pixlens = 0.200 * u.arcsec
if lensID == "WFI2033":
    z_s = 1.66
    z_l = 0.66
    limmag = 22.5
    limbright = 22.5
    brightmag = 16.90
    photoz = 'bpzeazy'
    pixnr = 915
    pixlens = 0.2625 * u.arcsec
if lensID == "J1206":
    z_s = 1.80
    z_l = 0.75
    limmag = 23
    limbright = 23
    brightmag = 18.05
    photoz = 'bpz'
    pixnr = 1283
    pixlens = 0.187 * u.arcsec

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
if lensID == "J1206":
    center_lens = SkyCoord('12:06:29.650 +43:32:19.90', frame='fk5', unit=(u.hourangle, u.deg))

lensbpz = np.loadtxt('%s/%s%sbpz_nobeta_%s.cat' % (rootlenscat,lensID,irac,det), unpack=True)
if photoz == 'bpzeazy': lenseazy = np.loadtxt('%s/%s%seazy_nobeta_%s.cat' % (rootlenscat,lensID,irac,det), unpack=True)
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
    #print np.shape(lenscat)
    lenscat = np.c_['0',lenscat,sep_lens.reshape((1, sep_lens.shape[0]))] # inserting as the last column of the catalogue
    #print np.shape(lenscat)
    lenscat = np.delete(lenscat,np.where(msk_lens[0].data[lenscat[y_lens].astype(int),lenscat[x_lens].astype(int)] != 0),axis=1) # remove the masked objects; I tested that this is the correct order of x and y. x and y in lenscat are the actual coordinates in the natural reading of the .fits file (not the reading of python, which inverts axes)
    #print msk_lens[0].data[400,485] # testing
    #print np.shape(lenscat)
    lenscat = np.delete(lenscat,np.where(lenscat[classify] < 0),axis=1) # removes all stars from the catalogue
    #print np.shape(lenscat)
    #print lenscat[RA_lens],lenscat[DEC_lens],SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5').separation(center_lens).arcsec,lenscat[classify]
    if remove != 'no':
        removecat = np.loadtxt('%s/%s.cat' % (rootlenscat,remove))
        for i in range(len(removecat)):
            coord_lensinit = SkyCoord(ra=lenscat[RA_lens]*u.degree, dec=lenscat[DEC_lens]*u.degree, frame='fk5') # redefine this whenever I remove rows
            coordremove = SkyCoord('%s %s' %(removecat[i][0],removecat[i][1]), frame='fk5', unit=(u.deg, u.deg))
            #print coordremove
            #print np.shape(lenscat)
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
if photoz == 'bpzeazy': lenseazy = lensprep(lenseazy)
#print lensbpz[RA_lens],lensbpz[DEC_lens],SkyCoord(ra=lensbpz[RA_lens]*u.degree, dec=lensbpz[DEC_lens]*u.degree, frame='fk5').separation(center_lens).arcsec,lensbpz[classify]

# Open the lens and field masks and set up the cell grid
print "Masking the lens catalogue corresponding to each CFHTLenS field cell..."

pixCFHT = 0.187 * u.arcsec
msk = fits.open('%s' % fmask)
worldfield = WCS('%s' % fmask)

if radius == 45: overlap = 2
if radius == 60: overlap = 3
if radius == 90: overlap = 4
if radius == 120: overlap = 5

mskfrac = msk_lens[0].data[int(round((pixnr - (2/pixlens.value) * radius)/2)) : int(round(pixnr - (pixnr - (2/pixlens.value) * radius)/2)), int(round((pixnr - (2/pixlens.value) * radius)/2)) : int(round(pixnr - (pixnr - (2/pixlens.value) * radius)/2))]  # will be used later in computing the mask cover fraction
unmaskedlens = np.where(mskfrac == 0) # 2D coordinates of mskfrac (the central subsection of msk_lens)
cells_on_a_side = int((len(msk[0].data[1]) * pixCFHT.value) / (2 * radius)) - 1
print "Cells on a side: ", cells_on_a_side
centerfieldspix = np.zeros((overlap,overlap,2))
unmaskedcell = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

# Declare the weighted counts:

if limmag != limbright: lens_gal_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_gal_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_zweight_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zweight_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass2_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass3_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_oneoverr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_zoverr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_zoverr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_massoverr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_massoverr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass2overr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass3overr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass2rms_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2rms_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass3rms_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3rms_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass2overrms_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass2overrms_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_mass3overrms_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_mass3overrms_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_flexion_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_flexion_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_tidal_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_tidal_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_convergence_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergence_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
if limmag != limbright: lens_convergencehalo_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_convergencehalo_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))

if photoz == 'bpzeazy':
    if limmag != limbright: lens_gal_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_gal_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_zweight_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_zweight_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass2_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass2_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass3_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass3_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_oneoverr_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_oneoverr_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_zoverr_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_zoverr_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_massoverr_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_massoverr_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass2overr_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass2overr_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass3overr_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass3overr_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass2rms_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass2rms_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass3rms_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass3rms_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass2overrms_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass2overrms_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_mass3overrms_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_mass3overrms_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_flexion_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_flexion_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_tidal_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_tidal_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_convergence_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_convergence_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    if limmag != limbright: lens_convergencehalo_limmageazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
    lens_convergencehalo_limbrighteazy = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))

print("Initialization of lens catalogue was completed in %0.1f seconds" % (time.time() - start_time))

######################################################################
# WORKING ON THE FIELD CATALOGUE
######################################################################

start_timefield = time.time()
print "Reading ", fieldID, "..."

field_ = Table.read('%s' % fieldID)
field = np.c_[field_['RA'],field_['DEC'],field_['photoz'],field_['i'],field_['y'],field_['mass_BEST'],field_['mass_MED']].T
#print np.shape(field)
del field_
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
if limmag != limbright: field_gal_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_gal_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_zweight_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zweight_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass2_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass3_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_oneoverr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_oneoverr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_zoverr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_zoverr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_massoverr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_massoverr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass2overr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass3overr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass2rms_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2rms_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass3rms_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3rms_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass2overrms_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass2overrms_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_mass3overrms_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_mass3overrms_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_flexion_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_flexion_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_tidal_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_tidal_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_convergence_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergence_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_convergencehalo_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_convergencehalo_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

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

                    if photoz == 'bpzeazy':
                        lenscoords = np.copy(lenseazy[x_lens:y_lens + 1]) # copy the pixel coordinates
                        lenscoords_fieldx = float(xlow)+(1.0 * lenscoords[1]/((2/pixlens.value) * radius))*(float(xhigh)-float(xlow))
                        lenscoords_fieldy = float(ylow)+(1.0 * lenscoords[0]/((2/pixlens.value) * radius))*(float(yhigh)-float(ylow))
                        lenscoords_field = msk[0].data[lenscoords_fieldx.astype(int),lenscoords_fieldy.astype(int)]
                        lenseazy_masked = np.c_['0',lenseazy,lenscoords_field.reshape(1,lenscoords_field.shape[0])] # check if the lens catalogue objects fall inside field masks
                        lenseazy_masked = np.delete(lenseazy_masked,np.where(lenseazy_masked[-1] != 0),axis=1) # remove objects inside a field mask
                        lenseazy_masked = lenseazy_masked[:-1] # delete the last column

                    for n in range(samples):
                        lensbpz_maskedlimmag = np.copy(lensbpz_masked)
                        if photoz == 'bpzeazy': lenseazy_maskedlimmag = np.copy(lenseazy_masked)
                        if n == 0:
                            if limmag != limbright:
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[i_lens] > limmag) | (lensbpz_maskedlimmag[z_lens] > z_s)),axis=1)
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[z_lens] >= zinf) & (lensbpz_maskedlimmag[z_lens] <= zsup)),axis=1) # remove the redshift slice
                                if photoz == 'bpzeazy': lenseazy_maskedlimmag = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[i_lens] > limmag) | (lenseazy_maskedlimmag[z_lens] > z_s)),axis=1)
                                if photoz == 'bpzeazy': lenseazy_maskedlimmag = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[z_lens] >= zinf) & (lenseazy_maskedlimmag[z_lens] <= zsup)),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[i_lens] > limbright) | (lensbpz_maskedlimmag[z_lens] > z_s)),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimbright,np.where((lensbpz_maskedlimbright[z_lens] >= zinf) & (lensbpz_maskedlimbright[z_lens] <= zsup)),axis=1)
                            if photoz == 'bpzeazy': lenseazy_maskedlimbright = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[i_lens] > limbright) | (lenseazy_maskedlimmag[z_lens] > z_s)),axis=1)
                            if photoz == 'bpzeazy': lenseazy_maskedlimbright = np.delete(lenseazy_maskedlimbright,np.where((lenseazy_maskedlimbright[z_lens] >= zinf) & (lenseazy_maskedlimbright[z_lens] <= zsup)),axis=1)
                        else:
                            for o in range(lensbpz_maskedlimmag.shape[1]):
                                lensbpz_maskedlimmag[i_lens][o] = np.random.normal(lensbpz_maskedlimmag[i_lens][o], np.max([lensbpz_maskedlimmag[i_err_lens][o],0.005]),1)[0]
                            if photoz == 'bpzeazy':
                                for o in range(lenseazy_maskedlimmag.shape[1]):
                                    lenseazy_maskedlimmag[i_lens][o] = np.random.normal(lenseazy_maskedlimmag[i_lens][o], np.max([lenseazy_maskedlimmag[i_err_lens][o],0.005]),1)[0]
                            if limmag != limbright:
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[i_lens] > limmag) | (lensbpz_maskedlimmag[z_lens + n * 3] > z_s)),axis=1)
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[z_lens + n * 3] >= zinf) & (lensbpz_maskedlimmag[z_lens + n * 3] <= zsup)),axis=1)  # remove the redshift slice
                                if photoz == 'bpzeazy': lenseazy_maskedlimmag = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[i_lens] > limmag) | (lenseazy_maskedlimmag[z_lens + n * 3] > z_s)),axis=1)
                                if photoz == 'bpzeazy': lenseazy_maskedlimmag = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[z_lens + n * 3] >= zinf) & (lenseazy_maskedlimmag[z_lens + n * 3] <= zsup)),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimmag,np.where((lensbpz_maskedlimmag[i_lens] > limbright) | (lensbpz_maskedlimmag[z_lens + n * 3] > z_s)),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimbright,np.where((lensbpz_maskedlimbright[z_lens + n * 3] >= zinf) & (lensbpz_maskedlimbright[z_lens + n * 3] <= zsup)),axis=1)
                            if photoz == 'bpzeazy': lenseazy_maskedlimbright = np.delete(lenseazy_maskedlimmag,np.where((lenseazy_maskedlimmag[i_lens] > limbright) | (lenseazy_maskedlimmag[z_lens + n * 3] > z_s)),axis=1)
                            if photoz == 'bpzeazy': lenseazy_maskedlimbright = np.delete(lenseazy_maskedlimbright,np.where((lenseazy_maskedlimbright[z_lens + n * 3] >= zinf) & (lenseazy_maskedlimbright[z_lens + n * 3] <= zsup)),axis=1)
                        if limmag != limbright:
                            lens_gal_limmagbpz,lens_zweight_limmagbpz,lens_mass_limmagbpz,lens_mass2_limmagbpz,lens_mass3_limmagbpz,lens_oneoverr_limmagbpz,lens_zoverr_limmagbpz,lens_massoverr_limmagbpz,lens_mass2overr_limmagbpz,lens_mass3overr_limmagbpz,lens_flexion_limmagbpz,lens_tidal_limmagbpz,lens_convergence_limmagbpz,lens_convergencehalo_limmagbpz,lens_mass2rms_limmagbpz,lens_mass3rms_limmagbpz,lens_mass2overrms_limmagbpz,lens_mass3overrms_limmagbpz = lensinit(lensbpz_maskedlimmag,lens_gal_limmagbpz,lens_zweight_limmagbpz,lens_mass_limmagbpz,lens_mass2_limmagbpz,lens_mass3_limmagbpz,lens_oneoverr_limmagbpz,lens_zoverr_limmagbpz,lens_massoverr_limmagbpz,lens_mass2overr_limmagbpz,lens_mass3overr_limmagbpz,lens_flexion_limmagbpz,lens_tidal_limmagbpz,lens_convergence_limmagbpz,lens_convergencehalo_limmagbpz,lens_mass2rms_limmagbpz,lens_mass3rms_limmagbpz,lens_mass2overrms_limmagbpz,lens_mass3overrms_limmagbpz)
                        #print "d ",lens_gal_limmagbpz[k][l][i][j][n]# test to check if the function actually returns the result globally
                        lens_gal_limbrightbpz,lens_zweight_limbrightbpz,lens_mass_limbrightbpz,lens_mass2_limbrightbpz,lens_mass3_limbrightbpz,lens_oneoverr_limbrightbpz,lens_zoverr_limbrightbpz,lens_massoverr_limbrightbpz,lens_mass2overr_limbrightbpz,lens_mass3overr_limbrightbpz,lens_flexion_limbrightbpz,lens_tidal_limbrightbpz,lens_convergence_limbrightbpz,lens_convergencehalo_limbrightbpz,lens_mass2rms_limbrightbpz,lens_mass3rms_limbrightbpz,lens_mass2overrms_limbrightbpz,lens_mass3overrms_limbrightbpz = lensinit(lensbpz_maskedlimbright,lens_gal_limbrightbpz,lens_zweight_limbrightbpz,lens_mass_limbrightbpz,lens_mass2_limbrightbpz,lens_mass3_limbrightbpz,lens_oneoverr_limbrightbpz,lens_zoverr_limbrightbpz,lens_massoverr_limbrightbpz,lens_mass2overr_limbrightbpz,lens_mass3overr_limbrightbpz,lens_flexion_limbrightbpz,lens_tidal_limbrightbpz,lens_convergence_limbrightbpz,lens_convergencehalo_limbrightbpz,lens_mass2rms_limbrightbpz,lens_mass3rms_limbrightbpz,lens_mass2overrms_limbrightbpz,lens_mass3overrms_limbrightbpz)
                        #print "d ",lens_gal_limbrightbpz[k][l][i][j][n]# test to check if the function actually returns the result globally
                        if photoz == 'bpzeazy':
                            if limmag != limbright:
                                lens_gal_limmageazy,lens_zweight_limmageazy,lens_mass_limmageazy,lens_mass2_limmageazy,lens_mass3_limmageazy,lens_oneoverr_limmageazy,lens_zoverr_limmageazy,lens_massoverr_limmageazy,lens_mass2overr_limmageazy,lens_mass3overr_limmageazy,lens_flexion_limmageazy,lens_tidal_limmageazy,lens_convergence_limmageazy,lens_convergencehalo_limmageazy,lens_mass2rms_limmageazy,lens_mass3rms_limmageazy,lens_mass2overrms_limmageazy,lens_mass3overrms_limmageazy = lensinit(lenseazy_maskedlimmag,lens_gal_limmageazy,lens_zweight_limmageazy,lens_mass_limmageazy,lens_mass2_limmageazy,lens_mass3_limmageazy,lens_oneoverr_limmageazy,lens_zoverr_limmageazy,lens_massoverr_limmageazy,lens_mass2overr_limmageazy,lens_mass3overr_limmageazy,lens_flexion_limmageazy,lens_tidal_limmageazy,lens_convergence_limmageazy,lens_convergencehalo_limmageazy,lens_mass2rms_limmageazy,lens_mass3rms_limmageazy,lens_mass2overrms_limmageazy,lens_mass3overrms_limmageazy)
                            #print "d ",lens_gal_limmageazy[k][l][i][j][n]# test to check if the function actually returns the result globally
                            lens_gal_limbrighteazy,lens_zweight_limbrighteazy,lens_mass_limbrighteazy,lens_mass2_limbrighteazy,lens_mass3_limbrighteazy,lens_oneoverr_limbrighteazy,lens_zoverr_limbrighteazy,lens_massoverr_limbrighteazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limbrighteazy,lens_flexion_limbrighteazy,lens_tidal_limbrighteazy,lens_convergence_limbrighteazy,lens_convergencehalo_limbrighteazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limbrighteazy = lensinit(lenseazy_maskedlimbright,lens_gal_limbrighteazy,lens_zweight_limbrighteazy,lens_mass_limbrighteazy,lens_mass2_limbrighteazy,lens_mass3_limbrighteazy,lens_oneoverr_limbrighteazy,lens_zoverr_limbrighteazy,lens_massoverr_limbrighteazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limbrighteazy,lens_flexion_limbrighteazy,lens_tidal_limbrighteazy,lens_convergence_limbrighteazy,lens_convergencehalo_limbrighteazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limbrighteazy)
                            #print "d ",lens_gal_limbrighteazy[k][l][i][j][n]# test to check if the function actually returns the result globally

                    '''Compute weights for the field catalogue'''
                    field_masked = np.copy(field)
                    field_masked = np.delete(field_masked,np.where((field_masked[cell_xpix] < xlow) | (field_masked[cell_xpix] > xhigh) | (field_masked[cell_ypix] < ylow) | (field_masked[cell_ypix] > yhigh)),axis=1) # isolate the field objects matching in position the field coordinates on which the lens catalogue is projected
                    #field_masked = np.delete(field_masked,np.where((field_masked[cell_xpix] > xlow))
                    field_masked[cell_xpix] = (field_masked[cell_xpix] - (centerfieldspix[k][l][0] - radius * (u.arcsec / pixCFHT).value)) * (pixCFHT / pixlens).value # the last multiplication brings the pixels to the lens pixel size; work on the physical image y axis
                    field_masked[cell_ypix] = (field_masked[cell_ypix] - (centerfieldspix[k][l][1] - radius * (u.arcsec / pixCFHT).value)) * (pixCFHT / pixlens).value
                    field_masked[cell_sep] = np.sqrt((field_masked[cell_xpix] - (radius * 2/pixlens.value)/2) ** 2 + (field_masked[cell_ypix] - (radius * 2/pixlens.value)/2) ** 2) * (pixlens / u.arcsec).value # separation in arcsec from the center of the cell
                    field_masked[cell_sep][field_masked[cell_sep] < 10] = 10
                    temp1 = np.array([field_masked[cell_xpix] + (pixnr - (2/pixlens.value) * radius)/2]).astype(int)
                    temp2 = np.array([field_masked[cell_ypix] + (pixnr - (2/pixlens.value) * radius)/2]).astype(int)
                    temp1[temp1==pixnr] = pixnr - 1
                    temp2[temp2==pixnr] = pixnr - 1
                    field_masked_limmag = np.delete(field_masked,np.where(msk_lens[0].data[temp1,temp2] != 0),axis=1) # remove field objects falling inside the lens mask; needed to account for the fact that the mask is always pixnr pixels on a side
                    #if np.max(field_masked_limmag[cell_sep]) > radius: print "a", np.max(field_masked_limmag[cell_sep]) # testing; due to the pixelated mask, there are objects very close to the radius limit but just very slightly (subarcsec) away; I will not do anything about this
                    field_masked_limmag = np.delete(field_masked_limmag,np.where(field_masked_limmag[i_field] > limmag),axis=1)
                    field_masked_limbright = np.delete(field_masked_limmag,np.where(field_masked_limmag[i_field] > limbright),axis=1)
                    w_gal_limmag = np.shape(field_masked_limmag)[1]
                    w_gal_limbright = np.shape(field_masked_limbright)[1]
                    #mmm = np.copy(field_convergencehalo_limbright)
                    if limmag != limbright:
                       field_gal_limmag,field_zweight_limmag,field_mass_limmag,field_mass2_limmag,field_mass3_limmag,field_oneoverr_limmag,field_zoverr_limmag,field_massoverr_limmag,field_mass2overr_limmag,field_mass3overr_limmag,field_mass2rms_limmag,field_mass3rms_limmag,field_mass2overrms_limmag,field_mass3overrms_limmag,field_flexion_limmag,field_tidal_limmag,field_convergence_limmag,field_convergencehalo_limmag = fieldinit(field_masked_limmag,w_gal_limmag,field_gal_limmag,field_zweight_limmag,field_mass_limmag,field_mass2_limmag,field_mass3_limmag,field_oneoverr_limmag,field_zoverr_limmag,field_massoverr_limmag,field_mass2overr_limmag,field_mass3overr_limmag,field_mass2rms_limmag,field_mass3rms_limmag,field_mass2overrms_limmag,field_mass3overrms_limmag,field_flexion_limmag,field_tidal_limmag,field_convergence_limmag,field_convergencehalo_limmag)
                    field_gal_limbright,field_zweight_limbright,field_mass_limbright,field_mass2_limbright,field_mass3_limbright,field_oneoverr_limbright,field_zoverr_limbright,field_massoverr_limbright,field_mass2overr_limbright,field_mass3overr_limbright,field_mass2rms_limbright,field_mass3rms_limbright,field_mass2overrms_limbright,field_mass3overrms_limbright,field_flexion_limbright,field_tidal_limbright,field_convergence_limbright,field_convergencehalo_limbright = fieldinit(field_masked_limbright,w_gal_limbright,field_gal_limbright,field_zweight_limbright,field_mass_limbright,field_mass2_limbright,field_mass3_limbright,field_oneoverr_limbright,field_zoverr_limbright,field_massoverr_limbright,field_mass2overr_limbright,field_mass3overr_limbright,field_mass2rms_limbright,field_mass3rms_limbright,field_mass2overrms_limbright,field_mass3overrms_limbright,field_flexion_limbright,field_tidal_limbright,field_convergence_limbright,field_convergencehalo_limbright)
                    #print np.min(mmm - field_convergencehalo_limbright)

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
                        #field_masked_limmag[RA_field],field_masked_limmag[DEC_field]
                        #np.savetxt("/Users/perseus/Desktop/CFHTLenSmasks/test_field_wcs.reg",np.c_[field_masked_limmag[RA_field],field_masked_limmag[DEC_field]],fmt='j2000; circle %s %s 2\"')

msk_lens.close()
msk.close()

print("Finishing lens and field catalogues in %0.1f seconds" % (time.time() - start_timeweights))

######################################################################
# WRITING OUTPUT
######################################################################

start_write = time.time()

print "Writing output..."

count = "Fields above .75 and .50 limits: %d %d, %d %d" % (unmaskedcell[unmaskedcell>=0.75].shape[0], unmaskedcell.shape[0] * unmaskedcell.shape[1] * unmaskedcell.shape[2] * unmaskedcell.shape[3], unmaskedcell[unmaskedcell>=0.5].shape[0], unmaskedcell.shape[0] * unmaskedcell.shape[1] * unmaskedcell.shape[2] * unmaskedcell.shape[3])
print count
mskname = 'msk%sarcsecrad%sarcsecgap.fits'[0:-5] % (radius,inner)
fcount = open('%s/%s_wghtratios_%s_%s_%s_%s_%s_zgap%s_%s%s_count.cat' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix),'w') # [-25:-4] corresponds to strings of the form W1m0m0_limmaggalphotmstar
fcount.write(count)
fcount.close()

fout50_bpzlimbright_0 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_0 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_1 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_1 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_2 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_2 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_3 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_3 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_4 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_4 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_5 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_5 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_6 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_6 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_7 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_7 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_8 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_8 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout50_bpzlimbright_9 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
fout75_bpzlimbright_9 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
#os.system('rm -f %s' % fout50_0) # '-f' ignores non-existent files
if photoz == 'bpzeazy':
    fout50_eazylimbright_0 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_0 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_1 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_1 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_2 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_2 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_3 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_3 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_4 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_4 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_5 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_5 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_6 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_6 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_7 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_7 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_8 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_8 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout50_eazylimbright_9 = '%s/%s_wghtratios_50_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
    fout75_eazylimbright_9 = '%s/%s_wghtratios_75_%s_%s_%s_eazy_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,det,irac,type,zinf,zsup,suffix)
if limmag != brightmag:
    fout50_bpzlimmag_0 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_0 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_1 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_1 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_2 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_2 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_3 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_3 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_4 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_4 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_5 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_5 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_6 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_6 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_7 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_7 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_8 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_8 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout50_bpzlimmag_9 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    fout75_bpzlimmag_9 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
    if photoz == 'bpzeazy':
        fout50_eazylimmag_0 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_0 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_1 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_1 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_2 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_2 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_3 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_3 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_4 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_4 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_5 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_5 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_6 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_6 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_7 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_7 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_8 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_8 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout50_eazylimmag_9 = '%s/%s_wghtratios_50_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)
        fout75_eazylimmag_9 = '%s/%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix)

def outputfunc(*argv):
    start0 = True; start1 = True; start2 = True; start3 = True; start4 = True; start5 = True; start6 = True; start7 = True; start8 = True; start9 = True
    if len(argv) == 109:
        '''limbright limmag bpz eazy'''
        # frac,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limmageazy,lens_gal_limbrighteazy,lens_zweight_limmageazy,lens_zweight_limbrighteazy,lens_mass_limmageazy,lens_mass_limbrighteazy,lens_mass2_limmageazy,lens_mass2_limbrighteazy,lens_mass3_limmageazy,lens_mass3_limbrighteazy,lens_oneoverr_limmageazy,lens_oneoverr_limbrighteazy,lens_zoverr_limmageazy,lens_zoverr_limbrighteazy,lens_massoverr_limmageazy,lens_massoverr_limbrighteazy,lens_mass2overr_limmageazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limmageazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limmageazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limmageazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limmageazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limmageazy,lens_mass3overrms_limbrighteazy,lens_flexion_limmageazy,lens_flexion_limbrighteazy,lens_tidal_limmageazy,lens_tidal_limbrighteazy,lens_convergence_limmageazy,lens_convergence_limbrighteazy,lens_convergencehalo_limmageazy,lens_convergencehalo_limbrighteazy
        frac_ = argv[0]; lens_gal_limmagbpz_ = argv[1]; field_gal_limmag_ = argv[2]; lens_gal_limbrightbpz_ = argv[3]; field_gal_limbright_ = argv[4]; lens_zweight_limmagbpz_ = argv[5]; field_zweight_limmag_ = argv[6]; lens_zweight_limbrightbpz_ = argv[7]; field_zweight_limbright_ = argv[8]; lens_mass_limmagbpz_ = argv[9]; field_mass_limmag_ = argv[10]; lens_mass_limbrightbpz_ = argv[11]; field_mass_limbright_ = argv[12]; lens_mass2_limmagbpz_ = argv[13]; field_mass2_limmag_ = argv[14]; lens_mass2_limbrightbpz_ = argv[15]; field_mass2_limbright_ = argv[16]; lens_mass3_limmagbpz_ = argv[17]; field_mass3_limmag_ = argv[18]; lens_mass3_limbrightbpz_ = argv[19]; field_mass3_limbright_ = argv[20]; lens_oneoverr_limmagbpz_ = argv[21]; field_oneoverr_limmag_ = argv[22]; lens_oneoverr_limbrightbpz_ = argv[23]; field_oneoverr_limbright_ = argv[24]; lens_zoverr_limmagbpz_ = argv[25]; field_zoverr_limmag_ = argv[26]; lens_zoverr_limbrightbpz_ = argv[27]; field_zoverr_limbright_ = argv[28]; lens_massoverr_limmagbpz_ = argv[29]; field_massoverr_limmag_ = argv[30]; lens_massoverr_limbrightbpz_ = argv[31]; field_massoverr_limbright_ = argv[32]; lens_mass2overr_limmagbpz_ = argv[33]; field_mass2overr_limmag_ = argv[34]; lens_mass2overr_limbrightbpz_ = argv[35]; field_mass2overr_limbright_ = argv[36]; lens_mass3overr_limmagbpz_ = argv[37]; field_mass3overr_limmag_ = argv[38]; lens_mass3overr_limbrightbpz_ = argv[39]; field_mass3overr_limbright_ = argv[40]; lens_mass2rms_limmagbpz_ = argv[41]; field_mass2rms_limmag_ = argv[42]; lens_mass2rms_limbrightbpz_ = argv[43]; field_mass2rms_limbright_ = argv[44]; lens_mass3rms_limmagbpz_ = argv[45]; field_mass3rms_limmag_ = argv[46]; lens_mass3rms_limbrightbpz_ = argv[47]; field_mass3rms_limbright_ = argv[48]; lens_mass2overrms_limmagbpz_ = argv[49]; field_mass2overrms_limmag_ = argv[50]; lens_mass2overrms_limbrightbpz_ = argv[51]; field_mass2overrms_limbright_ = argv[52]; lens_mass3overrms_limmagbpz_ = argv[53]; field_mass3overrms_limmag_ = argv[54]; lens_mass3overrms_limbrightbpz_ = argv[55]; field_mass3overrms_limbright_ = argv[56]; lens_flexion_limmagbpz_ = argv[57]; field_flexion_limmag_ = argv[58]; lens_flexion_limbrightbpz_ = argv[59]; field_flexion_limbright_ = argv[60]; lens_tidal_limmagbpz_ = argv[61]; field_tidal_limmag_ = argv[62]; lens_tidal_limbrightbpz_ = argv[63]; field_tidal_limbright_ = argv[64]; lens_convergence_limmagbpz_ = argv[65]; field_convergence_limmag_ = argv[66]; lens_convergence_limbrightbpz_ = argv[67]; field_convergence_limbright_ = argv[68]; lens_convergencehalo_limmagbpz_ = argv[69]; field_convergencehalo_limmag_ = argv[70]; lens_convergencehalo_limbrightbpz_ = argv[71]; field_convergencehalo_limbright_ = argv[72]; lens_gal_limmageazy_ = argv[73]; lens_gal_limbrighteazy_ = argv[74]; lens_zweight_limmageazy_ = argv[75]; lens_zweight_limbrighteazy_ = argv[76]; lens_mass_limmageazy_ = argv[77]; lens_mass_limbrighteazy_ = argv[78]; lens_mass2_limmageazy_ = argv[79]; lens_mass2_limbrighteazy_ = argv[80]; lens_mass3_limmageazy_ = argv[81]; lens_mass3_limbrighteazy_ = argv[82]; lens_oneoverr_limmageazy_ = argv[83]; lens_oneoverr_limbrighteazy_ = argv[84]; lens_zoverr_limmageazy_ = argv[85]; lens_zoverr_limbrighteazy_ = argv[86]; lens_massoverr_limmageazy_ = argv[87]; lens_massoverr_limbrighteazy_ = argv[88]; lens_mass2overr_limmageazy_ = argv[89]; lens_mass2overr_limbrighteazy_ = argv[90]; lens_mass3overr_limmageazy_ = argv[91]; lens_mass3overr_limbrighteazy_ = argv[92]; lens_mass2rms_limmageazy_ = argv[93]; lens_mass2rms_limbrighteazy_ = argv[94]; lens_mass3rms_limmageazy_ = argv[95]; lens_mass3rms_limbrighteazy_ = argv[96]; lens_mass2overrms_limmageazy_ = argv[97]; lens_mass2overrms_limbrighteazy_ = argv[98]; lens_mass3overrms_limmageazy_ = argv[99]; lens_mass3overrms_limbrighteazy_ = argv[100]; lens_flexion_limmageazy_ = argv[101]; lens_flexion_limbrighteazy_ = argv[102]; lens_tidal_limmageazy_ = argv[103]; lens_tidal_limbrighteazy_ = argv[104]; lens_convergence_limmageazy_ = argv[105]; lens_convergence_limbrighteazy_ = argv[106]; lens_convergencehalo_limmageazy_ = argv[107]; lens_convergencehalo_limbrighteazy_ = argv[108]
        outbpzlimbright = np.zeros(22)
        outbpzlimmag = np.zeros(22)
        outeazylimmag = np.zeros(22)
        outeazylimbright = np.zeros(22)
    if len(argv) == 73:
        '''limbright limmag bpz'''
        # frac,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright
        frac_ = argv[0]; lens_gal_limmagbpz_ = argv[1]; field_gal_limmag_ = argv[2]; lens_gal_limbrightbpz_ = argv[3]; field_gal_limbright_ = argv[4]; lens_zweight_limmagbpz_ = argv[5]; field_zweight_limmag_ = argv[6]; lens_zweight_limbrightbpz_ = argv[7]; field_zweight_limbright_ = argv[8]; lens_mass_limmagbpz_ = argv[9]; field_mass_limmag_ = argv[10]; lens_mass_limbrightbpz_ = argv[11]; field_mass_limbright_ = argv[12]; lens_mass2_limmagbpz_ = argv[13]; field_mass2_limmag_ = argv[14]; lens_mass2_limbrightbpz_ = argv[15]; field_mass2_limbright_ = argv[16]; lens_mass3_limmagbpz_ = argv[17]; field_mass3_limmag_ = argv[18]; lens_mass3_limbrightbpz_ = argv[19]; field_mass3_limbright_ = argv[20]; lens_oneoverr_limmagbpz_ = argv[21]; field_oneoverr_limmag_ = argv[22]; lens_oneoverr_limbrightbpz_ = argv[23]; field_oneoverr_limbright_ = argv[24]; lens_zoverr_limmagbpz_ = argv[25]; field_zoverr_limmag_ = argv[26]; lens_zoverr_limbrightbpz_ = argv[27]; field_zoverr_limbright_ = argv[28]; lens_massoverr_limmagbpz_ = argv[29]; field_massoverr_limmag_ = argv[30]; lens_massoverr_limbrightbpz_ = argv[31]; field_massoverr_limbright_ = argv[32]; lens_mass2overr_limmagbpz_ = argv[33]; field_mass2overr_limmag_ = argv[34]; lens_mass2overr_limbrightbpz_ = argv[35]; field_mass2overr_limbright_ = argv[36]; lens_mass3overr_limmagbpz_ = argv[37]; field_mass3overr_limmag_ = argv[38]; lens_mass3overr_limbrightbpz_ = argv[39]; field_mass3overr_limbright_ = argv[40]; lens_mass2rms_limmagbpz_ = argv[41]; field_mass2rms_limmag_ = argv[42]; lens_mass2rms_limbrightbpz_ = argv[43]; field_mass2rms_limbright_ = argv[44]; lens_mass3rms_limmagbpz_ = argv[45]; field_mass3rms_limmag_ = argv[46]; lens_mass3rms_limbrightbpz_ = argv[47]; field_mass3rms_limbright_ = argv[48]; lens_mass2overrms_limmagbpz_ = argv[49]; field_mass2overrms_limmag_ = argv[50]; lens_mass2overrms_limbrightbpz_ = argv[51]; field_mass2overrms_limbright_ = argv[52]; lens_mass3overrms_limmagbpz_ = argv[53]; field_mass3overrms_limmag_ = argv[54]; lens_mass3overrms_limbrightbpz_ = argv[55]; field_mass3overrms_limbright_ = argv[56]; lens_flexion_limmagbpz_ = argv[57]; field_flexion_limmag_ = argv[58]; lens_flexion_limbrightbpz_ = argv[59]; field_flexion_limbright_ = argv[60]; lens_tidal_limmagbpz_ = argv[61]; field_tidal_limmag_ = argv[62]; lens_tidal_limbrightbpz_ = argv[63]; field_tidal_limbright_ = argv[64]; lens_convergence_limmagbpz_ = argv[65]; field_convergence_limmag_ = argv[66]; lens_convergence_limbrightbpz_ = argv[67]; field_convergence_limbright_ = argv[68]; lens_convergencehalo_limmagbpz_ = argv[69]; field_convergencehalo_limmag_ = argv[70]; lens_convergencehalo_limbrightbpz_ = argv[71]; field_convergencehalo_limbright_ = argv[72]
        outbpzlimbright = np.zeros(22)
        outbpzlimmag = np.zeros(22)
    if len(argv) == 55:
        '''limbright bpz eazy'''
        # frac,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limbrighteazy,lens_zweight_limbrighteazy,lens_mass_limbrighteazy,lens_mass2_limbrighteazy,lens_mass3_limbrighteazy,lens_oneoverr_limbrighteazy,lens_zoverr_limbrighteazy,lens_massoverr_limbrighteazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limbrighteazy,lens_flexion_limbrighteazy,lens_tidal_limbrighteazy,lens_convergence_limbrighteazy,lens_convergencehalo_limbrighteazy
        frac_ = argv[0]; lens_gal_limbrightbpz_ = argv[1]; field_gal_limbright_ = argv[2]; lens_zweight_limbrightbpz_ = argv[3]; field_zweight_limbright_ = argv[4]; lens_mass_limbrightbpz_ = argv[5]; field_mass_limbright_ = argv[6]; lens_mass2_limbrightbpz_ = argv[7]; field_mass2_limbright_ = argv[8]; lens_mass3_limbrightbpz_ = argv[9]; field_mass3_limbright_ = argv[10]; lens_oneoverr_limbrightbpz_ = argv[11]; field_oneoverr_limbright_ = argv[12]; lens_zoverr_limbrightbpz_ = argv[13]; field_zoverr_limbright_ = argv[14]; lens_massoverr_limbrightbpz_ = argv[15]; field_massoverr_limbright_ = argv[16]; lens_mass2overr_limbrightbpz_ = argv[17]; field_mass2overr_limbright_ = argv[18]; lens_mass3overr_limbrightbpz_ = argv[19]; field_mass3overr_limbright_ = argv[20]; lens_mass2rms_limbrightbpz_ = argv[21]; field_mass2rms_limbright_ = argv[22]; lens_mass3rms_limbrightbpz_ = argv[23]; field_mass3rms_limbright_ = argv[24]; lens_mass2overrms_limbrightbpz_ = argv[25]; field_mass2overrms_limbright_ = argv[26]; lens_mass3overrms_limbrightbpz_ = argv[27]; field_mass3overrms_limbright_ = argv[28]; lens_flexion_limbrightbpz_ = argv[29]; field_flexion_limbright_ = argv[30]; lens_tidal_limbrightbpz_ = argv[31]; field_tidal_limbright_ = argv[32]; lens_convergence_limbrightbpz_ = argv[33]; field_convergence_limbright_ = argv[34]; lens_convergencehalo_limbrightbpz_ = argv[35]; field_convergencehalo_limbright_ = argv[36]; lens_gal_limbrighteazy_ = argv[37]; lens_zweight_limbrighteazy_ = argv[38]; lens_mass_limbrighteazy_ = argv[39]; lens_mass2_limbrighteazy_ = argv[40]; lens_mass3_limbrighteazy_ = argv[41]; lens_oneoverr_limbrighteazy_ = argv[42]; lens_zoverr_limbrighteazy_ = argv[43]; lens_massoverr_limbrighteazy_ = argv[44]; lens_mass2overr_limbrighteazy_ = argv[45]; lens_mass3overr_limbrighteazy_ = argv[46]; lens_mass2rms_limbrighteazy_ = argv[47]; lens_mass3rms_limbrighteazy_ = argv[48]; lens_mass2overrms_limbrighteazy_ = argv[49]; lens_mass3overrms_limbrighteazy_ = argv[50]; lens_flexion_limbrighteazy_ = argv[51]; lens_tidal_limbrighteazy_ = argv[52]; lens_convergence_limbrighteazy_ = argv[53]; lens_convergencehalo_limbrighteazy_ = argv[54]
        outbpzlimbright = np.zeros(22)
        outeazylimbright = np.zeros(22)
    if len(argv) == 37:
        '''limbright bpz'''
        # frac,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright
        frac_ = argv[0]; lens_gal_limbrightbpz_ = argv[1]; field_gal_limbright_ = argv[2]; lens_zweight_limbrightbpz_ = argv[3]; field_zweight_limbright_ = argv[4]; lens_mass_limbrightbpz_ = argv[5]; field_mass_limbright_ = argv[6]; lens_mass2_limbrightbpz_ = argv[7]; field_mass2_limbright_ = argv[8]; lens_mass3_limbrightbpz_ = argv[9]; field_mass3_limbright_ = argv[10]; lens_oneoverr_limbrightbpz_ = argv[11]; field_oneoverr_limbright_ = argv[12]; lens_zoverr_limbrightbpz_ = argv[13]; field_zoverr_limbright_ = argv[14]; lens_massoverr_limbrightbpz_ = argv[15]; field_massoverr_limbright_ = argv[16]; lens_mass2overr_limbrightbpz_ = argv[17]; field_mass2overr_limbright_ = argv[18]; lens_mass3overr_limbrightbpz_ = argv[19]; field_mass3overr_limbright_ = argv[20]; lens_mass2rms_limbrightbpz_ = argv[21]; field_mass2rms_limbright_ = argv[22]; lens_mass3rms_limbrightbpz_ = argv[23]; field_mass3rms_limbright_ = argv[24]; lens_mass2overrms_limbrightbpz_ = argv[25]; field_mass2overrms_limbright_ = argv[26]; lens_mass3overrms_limbrightbpz_ = argv[27]; field_mass3overrms_limbright_ = argv[28]; lens_flexion_limbrightbpz_ = argv[29]; field_flexion_limbright_ = argv[30]; lens_tidal_limbrightbpz_ = argv[31]; field_tidal_limbright_ = argv[32]; lens_convergence_limbrightbpz_ = argv[33]; field_convergence_limbright_ = argv[34]; lens_convergencehalo_limbrightbpz_ = argv[35]; field_convergencehalo_limbright_ = argv[36]
        outbpzlimbright = np.zeros(22)

    for k in range(overlap):
        for l in range(overlap):
            for i in range(cells_on_a_side):
                for j in range(cells_on_a_side):
                    for n in range(10):
                        condition = False
                        if len(argv) == 109:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limmagbpz_[k][l][i][j]) != 0) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (np.min(lens_gal_limmageazy_[k][l][i][j]) != 0) & (np.min(lens_gal_limbrighteazy_[k][l][i][j]) != 0) & (field_gal_limmag_[k][l][i][j] != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if len(argv) == 73:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limmagbpz_[k][l][i][j]) != 0) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (field_gal_limmag_[k][l][i][j] != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if len(argv) == 55:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (np.min(lens_gal_limbrighteazy_[k][l][i][j]) != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if len(argv) == 37:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if condition == True:
                            if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109:
                                outbpzlimbright[0] = np.int16(k)
                                outbpzlimbright[1] = np.int16(l)
                                outbpzlimbright[2] = np.int16(i)
                                outbpzlimbright[3] = np.int16(j)
                                outbpzlimbright[4] = np.float32(1.0*lens_gal_limbrightbpz_[k][l][i][j][n]/field_gal_limbright_[k][l][i][j])
                                outbpzlimbright[5] = np.float32(1.0*lens_zweight_limbrightbpz_[k][l][i][j][n]/field_zweight_limbright_[k][l][i][j])
                                outbpzlimbright[6] = np.float32(1.0*lens_mass_limbrightbpz_[k][l][i][j][n]/field_mass_limbright_[k][l][i][j])
                                outbpzlimbright[7] = np.float32(1.0*lens_mass2_limbrightbpz_[k][l][i][j][n]/field_mass2_limbright_[k][l][i][j])
                                outbpzlimbright[8] = np.float32(1.0*lens_mass3_limbrightbpz_[k][l][i][j][n]/field_mass3_limbright_[k][l][i][j])
                                outbpzlimbright[9] = np.float32(1.0*lens_oneoverr_limbrightbpz_[k][l][i][j][n]/field_oneoverr_limbright_[k][l][i][j])
                                outbpzlimbright[10] = np.float32(1.0*lens_zoverr_limbrightbpz_[k][l][i][j][n]/field_zoverr_limbright_[k][l][i][j])
                                outbpzlimbright[11] = np.float32(1.0*lens_massoverr_limbrightbpz_[k][l][i][j][n]/field_massoverr_limbright_[k][l][i][j])
                                outbpzlimbright[12] = np.float32(1.0*lens_mass2overr_limbrightbpz_[k][l][i][j][n]/field_mass2overr_limbright_[k][l][i][j])
                                outbpzlimbright[13] = np.float32(1.0*lens_mass3overr_limbrightbpz_[k][l][i][j][n]/field_mass3overr_limbright_[k][l][i][j])
                                outbpzlimbright[14] = np.float32(1.0*lens_mass2rms_limbrightbpz_[k][l][i][j][n]/field_mass2rms_limbright_[k][l][i][j])
                                outbpzlimbright[15] = np.float32(1.0*lens_mass3rms_limbrightbpz_[k][l][i][j][n]/field_mass3rms_limbright_[k][l][i][j])
                                outbpzlimbright[16] = np.float32(1.0*lens_mass2overrms_limbrightbpz_[k][l][i][j][n]/field_mass2overrms_limbright_[k][l][i][j])
                                outbpzlimbright[17] = np.float32(1.0*lens_mass3overrms_limbrightbpz_[k][l][i][j][n]/field_mass3overrms_limbright_[k][l][i][j])
                                outbpzlimbright[18] = np.float32(1.0*lens_flexion_limbrightbpz_[k][l][i][j][n]/field_flexion_limbright_[k][l][i][j])
                                outbpzlimbright[19] = np.float32(1.0*lens_tidal_limbrightbpz_[k][l][i][j][n]/field_tidal_limbright_[k][l][i][j])
                                outbpzlimbright[20] = np.float32(1.0*lens_convergence_limbrightbpz_[k][l][i][j][n]/field_convergence_limbright_[k][l][i][j])
                                outbpzlimbright[21] = np.float32(1.0*lens_convergencehalo_limbrightbpz_[k][l][i][j][n]/field_convergencehalo_limbright_[k][l][i][j])
                            if len(argv) == 55 or len(argv) == 109:
                                outeazylimbright[0] = np.int16(k)
                                outeazylimbright[1] = np.int16(l)
                                outeazylimbright[2] = np.int16(i)
                                outeazylimbright[3] = np.int16(j)
                                outeazylimbright[4] = np.float32(1.0*lens_gal_limbrighteazy_[k][l][i][j][n]/field_gal_limbright_[k][l][i][j])
                                outeazylimbright[5] = np.float32(1.0*lens_zweight_limbrighteazy_[k][l][i][j][n]/field_zweight_limbright_[k][l][i][j])
                                outeazylimbright[6] = np.float32(1.0*lens_mass_limbrighteazy_[k][l][i][j][n]/field_mass_limbright_[k][l][i][j])
                                outeazylimbright[7] = np.float32(1.0*lens_mass2_limbrighteazy_[k][l][i][j][n]/field_mass2_limbright_[k][l][i][j])
                                outeazylimbright[8] = np.float32(1.0*lens_mass3_limbrighteazy_[k][l][i][j][n]/field_mass3_limbright_[k][l][i][j])
                                outeazylimbright[9] = np.float32(1.0*lens_oneoverr_limbrighteazy_[k][l][i][j][n]/field_oneoverr_limbright_[k][l][i][j])
                                outeazylimbright[10] = np.float32(1.0*lens_zoverr_limbrighteazy_[k][l][i][j][n]/field_zoverr_limbright_[k][l][i][j])
                                outeazylimbright[11] = np.float32(1.0*lens_massoverr_limbrighteazy_[k][l][i][j][n]/field_massoverr_limbright_[k][l][i][j])
                                outeazylimbright[12] = np.float32(1.0*lens_mass2overr_limbrighteazy_[k][l][i][j][n]/field_mass2overr_limbright_[k][l][i][j])
                                outeazylimbright[13] = np.float32(1.0*lens_mass3overr_limbrighteazy_[k][l][i][j][n]/field_mass3overr_limbright_[k][l][i][j])
                                outeazylimbright[14] = np.float32(1.0*lens_mass2rms_limbrighteazy_[k][l][i][j][n]/field_mass2rms_limbright_[k][l][i][j])
                                outeazylimbright[15] = np.float32(1.0*lens_mass3rms_limbrighteazy_[k][l][i][j][n]/field_mass3rms_limbright_[k][l][i][j])
                                outeazylimbright[16] = np.float32(1.0*lens_mass2overrms_limbrighteazy_[k][l][i][j][n]/field_mass2overrms_limbright_[k][l][i][j])
                                outeazylimbright[17] = np.float32(1.0*lens_mass3overrms_limbrighteazy_[k][l][i][j][n]/field_mass3overrms_limbright_[k][l][i][j])
                                outeazylimbright[18] = np.float32(1.0*lens_flexion_limbrighteazy_[k][l][i][j][n]/field_flexion_limbright_[k][l][i][j])
                                outeazylimbright[19] = np.float32(1.0*lens_tidal_limbrighteazy_[k][l][i][j][n]/field_tidal_limbright_[k][l][i][j])
                                outeazylimbright[20] = np.float32(1.0*lens_convergence_limbrighteazy_[k][l][i][j][n]/field_convergence_limbright_[k][l][i][j])
                                outeazylimbright[21] = np.float32(1.0*lens_convergencehalo_limbrighteazy_[k][l][i][j][n]/field_convergencehalo_limbright_[k][l][i][j])
                            if len(argv) == 73 or len(argv) == 109:
                                outbpzlimmag[0] = np.int16(k)
                                outbpzlimmag[1] = np.int16(l)
                                outbpzlimmag[2] = np.int16(i)
                                outbpzlimmag[3] = np.int16(j)
                                outbpzlimmag[4] = np.float32(1.0*lens_gal_limmagbpz_[k][l][i][j][n]/field_gal_limmag_[k][l][i][j])
                                outbpzlimmag[5] = np.float32(1.0*lens_zweight_limmagbpz_[k][l][i][j][n]/field_zweight_limmag_[k][l][i][j])
                                outbpzlimmag[6] = np.float32(1.0*lens_mass_limmagbpz_[k][l][i][j][n]/field_mass_limmag_[k][l][i][j])
                                outbpzlimmag[7] = np.float32(1.0*lens_mass2_limmagbpz_[k][l][i][j][n]/field_mass2_limmag_[k][l][i][j])
                                outbpzlimmag[8] = np.float32(1.0*lens_mass3_limmagbpz_[k][l][i][j][n]/field_mass3_limmag_[k][l][i][j])
                                outbpzlimmag[9] = np.float32(1.0*lens_oneoverr_limmagbpz_[k][l][i][j][n]/field_oneoverr_limmag_[k][l][i][j])
                                outbpzlimmag[10] = np.float32(1.0*lens_zoverr_limmagbpz_[k][l][i][j][n]/field_zoverr_limmag_[k][l][i][j])
                                outbpzlimmag[11] = np.float32(1.0*lens_massoverr_limmagbpz_[k][l][i][j][n]/field_massoverr_limmag_[k][l][i][j])
                                outbpzlimmag[12] = np.float32(1.0*lens_mass2overr_limmagbpz_[k][l][i][j][n]/field_mass2overr_limmag_[k][l][i][j])
                                outbpzlimmag[13] = np.float32(1.0*lens_mass3overr_limmagbpz_[k][l][i][j][n]/field_mass3overr_limmag_[k][l][i][j])
                                outbpzlimmag[14] = np.float32(1.0*lens_mass2rms_limmagbpz_[k][l][i][j][n]/field_mass2rms_limmag_[k][l][i][j])
                                outbpzlimmag[15] = np.float32(1.0*lens_mass3rms_limmagbpz_[k][l][i][j][n]/field_mass3rms_limmag_[k][l][i][j])
                                outbpzlimmag[16] = np.float32(1.0*lens_mass2overrms_limmagbpz_[k][l][i][j][n]/field_mass2overrms_limmag_[k][l][i][j])
                                outbpzlimmag[17] = np.float32(1.0*lens_mass3overrms_limmagbpz_[k][l][i][j][n]/field_mass3overrms_limmag_[k][l][i][j])
                                outbpzlimmag[18] = np.float32(1.0*lens_flexion_limmagbpz_[k][l][i][j][n]/field_flexion_limmag_[k][l][i][j])
                                outbpzlimmag[19] = np.float32(1.0*lens_tidal_limmagbpz_[k][l][i][j][n]/field_tidal_limmag_[k][l][i][j])
                                outbpzlimmag[20] = np.float32(1.0*lens_convergence_limmagbpz_[k][l][i][j][n]/field_convergence_limmag_[k][l][i][j])
                                outbpzlimmag[21] = np.float32(1.0*lens_convergencehalo_limmagbpz_[k][l][i][j][n]/field_convergencehalo_limmag_[k][l][i][j])
                            if len(argv) == 109:
                                outeazylimmag[0] = np.int16(k)
                                outeazylimmag[1] = np.int16(l)
                                outeazylimmag[2] = np.int16(i)
                                outeazylimmag[3] = np.int16(j)
                                outeazylimmag[4] = np.float32(1.0*lens_gal_limmageazy_[k][l][i][j][n]/field_gal_limmag_[k][l][i][j])
                                outeazylimmag[5] = np.float32(1.0*lens_zweight_limmageazy_[k][l][i][j][n]/field_zweight_limmag_[k][l][i][j])
                                outeazylimmag[6] = np.float32(1.0*lens_mass_limmageazy_[k][l][i][j][n]/field_mass_limmag_[k][l][i][j])
                                outeazylimmag[7] = np.float32(1.0*lens_mass2_limmageazy_[k][l][i][j][n]/field_mass2_limmag_[k][l][i][j])
                                outeazylimmag[8] = np.float32(1.0*lens_mass3_limmageazy_[k][l][i][j][n]/field_mass3_limmag_[k][l][i][j])
                                outeazylimmag[9] = np.float32(1.0*lens_oneoverr_limmageazy_[k][l][i][j][n]/field_oneoverr_limmag_[k][l][i][j])
                                outeazylimmag[10] = np.float32(1.0*lens_zoverr_limmageazy_[k][l][i][j][n]/field_zoverr_limmag_[k][l][i][j])
                                outeazylimmag[11] = np.float32(1.0*lens_massoverr_limmageazy_[k][l][i][j][n]/field_massoverr_limmag_[k][l][i][j])
                                outeazylimmag[12] = np.float32(1.0*lens_mass2overr_limmageazy_[k][l][i][j][n]/field_mass2overr_limmag_[k][l][i][j])
                                outeazylimmag[13] = np.float32(1.0*lens_mass3overr_limmageazy_[k][l][i][j][n]/field_mass3overr_limmag_[k][l][i][j])
                                outeazylimmag[14] = np.float32(1.0*lens_mass2rms_limmageazy_[k][l][i][j][n]/field_mass2rms_limmag_[k][l][i][j])
                                outeazylimmag[15] = np.float32(1.0*lens_mass3rms_limmageazy_[k][l][i][j][n]/field_mass3rms_limmag_[k][l][i][j])
                                outeazylimmag[16] = np.float32(1.0*lens_mass2overrms_limmageazy_[k][l][i][j][n]/field_mass2overrms_limmag_[k][l][i][j])
                                outeazylimmag[17] = np.float32(1.0*lens_mass3overrms_limmageazy_[k][l][i][j][n]/field_mass3overrms_limmag_[k][l][i][j])
                                outeazylimmag[18] = np.float32(1.0*lens_flexion_limmageazy_[k][l][i][j][n]/field_flexion_limmag_[k][l][i][j])
                                outeazylimmag[19] = np.float32(1.0*lens_tidal_limmageazy_[k][l][i][j][n]/field_tidal_limmag_[k][l][i][j])
                                outeazylimmag[20] = np.float32(1.0*lens_convergence_limmageazy_[k][l][i][j][n]/field_convergence_limmag_[k][l][i][j])
                                outeazylimmag[21] = np.float32(1.0*lens_convergencehalo_limmageazy_[k][l][i][j][n]/field_convergencehalo_limmag_[k][l][i][j])
                            if n == 0:
                                if start0 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_0 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_0 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_0 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_0 = outeazylimmag
                                    start0 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_0 = np.c_[outlistbpzlimbright_0,outbpzlimbright] # normally this line should start with else, but for some reason that makes it lose the first row and duplicate the second; this way it will duplicate the first line, and at the end I will mask its first instance
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_0 = np.c_[outlisteazylimbright_0,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_0 = np.c_[outlistbpzlimmag_0,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_0 = np.c_[outlisteazylimmag_0,outeazylimmag]
                            if n == 1:
                                if start1 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_1 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_1 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_1 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_1 = outeazylimmag
                                    start1 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_1 = np.c_[outlistbpzlimbright_1,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_1 = np.c_[outlisteazylimbright_1,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_1 = np.c_[outlistbpzlimmag_1,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_1 = np.c_[outlisteazylimmag_1,outeazylimmag]
                            if n == 2:
                                if start2 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_2 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_2 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_2 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_2 = outeazylimmag
                                    start2 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_2 = np.c_[outlistbpzlimbright_2,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_2 = np.c_[outlisteazylimbright_2,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_2 = np.c_[outlistbpzlimmag_2,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_2 = np.c_[outlisteazylimmag_2,outeazylimmag]
                            if n == 3:
                                if start3 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_3 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_3 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_3 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_3 = outeazylimmag
                                    start3 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_3 = np.c_[outlistbpzlimbright_3,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_3 = np.c_[outlisteazylimbright_3,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_3 = np.c_[outlistbpzlimmag_3,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_3 = np.c_[outlisteazylimmag_3,outeazylimmag]
                            if n == 4:
                                if start4 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_4 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_4 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_4 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_4 = outeazylimmag
                                    start4 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_4 = np.c_[outlistbpzlimbright_4,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_4 = np.c_[outlisteazylimbright_4,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_4 = np.c_[outlistbpzlimmag_4,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_4 = np.c_[outlisteazylimmag_4,outeazylimmag]
                            if n == 5:
                                if start5 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_5 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_5 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_5 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_5 = outeazylimmag
                                    start5 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_5 = np.c_[outlistbpzlimbright_5,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_5 = np.c_[outlisteazylimbright_5,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_5 = np.c_[outlistbpzlimmag_5,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_5 = np.c_[outlisteazylimmag_5,outeazylimmag]
                            if n == 6:
                                if start6 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_6 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_6 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_6 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_6 = outeazylimmag
                                    start6 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_6 = np.c_[outlistbpzlimbright_6,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_6 = np.c_[outlisteazylimbright_6,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_6 = np.c_[outlistbpzlimmag_6,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_6 = np.c_[outlisteazylimmag_6,outeazylimmag]
                            if n == 7:
                                if start7 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_7 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_7 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_7 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_7 = outeazylimmag
                                    start7 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_7 = np.c_[outlistbpzlimbright_7,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_7 = np.c_[outlisteazylimbright_7,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_7 = np.c_[outlistbpzlimmag_7,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_7 = np.c_[outlisteazylimmag_7,outeazylimmag]
                            if n == 8:
                                if start8 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_8 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_8 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_8 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_8 = outeazylimmag
                                    start8 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_8 = np.c_[outlistbpzlimbright_8,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_8 = np.c_[outlisteazylimbright_8,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_8 = np.c_[outlistbpzlimmag_8,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_8 = np.c_[outlisteazylimmag_8,outeazylimmag]
                            if n == 9:
                                if start9 == True:
                                    if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_9 = outbpzlimbright
                                    if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_9 = outeazylimbright
                                    if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_9 = outbpzlimmag
                                    if len(argv) == 109: outlisteazylimmag_9 = outeazylimmag
                                    start9 = False
                                if len(argv) == 37 or len(argv) == 55 or len(argv) == 73 or len(argv) == 109: outlistbpzlimbright_9 = np.c_[outlistbpzlimbright_9,outbpzlimbright]
                                if len(argv) == 55 or len(argv) == 109: outlisteazylimbright_9 = np.c_[outlisteazylimbright_9,outeazylimbright]
                                if len(argv) == 73 or len(argv) == 109: outlistbpzlimmag_9 = np.c_[outlistbpzlimmag_9,outbpzlimmag]
                                if len(argv) == 109: outlisteazylimmag_9 = np.c_[outlisteazylimmag_9,outeazylimmag]

    if len(argv) == 109: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:],outlisteazylimbright_0[:,1:],outlisteazylimbright_1[:,1:],outlisteazylimbright_2[:,1:],outlisteazylimbright_3[:,1:],outlisteazylimbright_4[:,1:],outlisteazylimbright_5[:,1:],outlisteazylimbright_6[:,1:],outlisteazylimbright_7[:,1:],outlisteazylimbright_8[:,1:],outlisteazylimbright_9[:,1:],outlistbpzlimmag_0[:,1:],outlistbpzlimmag_1[:,1:],outlistbpzlimmag_2[:,1:],outlistbpzlimmag_3[:,1:],outlistbpzlimmag_4[:,1:],outlistbpzlimmag_5[:,1:],outlistbpzlimmag_6[:,1:],outlistbpzlimmag_7[:,1:],outlistbpzlimmag_8[:,1:],outlistbpzlimmag_9[:,1:],outlisteazylimmag_0[:,1:],outlisteazylimmag_1[:,1:],outlisteazylimmag_2[:,1:],outlisteazylimmag_3[:,1:],outlisteazylimmag_4[:,1:],outlisteazylimmag_5[:,1:],outlisteazylimmag_6[:,1:],outlisteazylimmag_7[:,1:],outlisteazylimmag_8[:,1:],outlisteazylimmag_9[:,1:]
    if len(argv) == 73: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:],outlistbpzlimmag_0[:,1:],outlistbpzlimmag_1[:,1:],outlistbpzlimmag_2[:,1:],outlistbpzlimmag_3[:,1:],outlistbpzlimmag_4[:,1:],outlistbpzlimmag_5[:,1:],outlistbpzlimmag_6[:,1:],outlistbpzlimmag_7[:,1:],outlistbpzlimmag_8[:,1:],outlistbpzlimmag_9[:,1:]
    if len(argv) == 55: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:],outlisteazylimbright_0[:,1:],outlisteazylimbright_1[:,1:],outlisteazylimbright_2[:,1:],outlisteazylimbright_3[:,1:],outlisteazylimbright_4[:,1:],outlisteazylimbright_5[:,1:],outlisteazylimbright_6[:,1:],outlisteazylimbright_7[:,1:],outlisteazylimbright_8[:,1:],outlisteazylimbright_9[:,1:]
    if len(argv) == 37: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:]

names = ('1_overlap_x','2_overlap_y','3_cell_x','4_cell_y','5_lens_gal','6_lens_zweight','7_lens_mass','8_lens_mass2','9_lens_mass3','10_lens_oneoverr','11_lens_zoverr','12_lens_massoverr','13_lens_mass2overr','14_lens_mass3overr','15_lens_mass2rms','16_lens_mass3rms','17_lens_mass2overrms','18_lens_mass3overrms','19_lens_flexion','20_lens_tidal','21_lens_convergence','22_lens_convergencehalo')
dtype=(np.int16,np.int16,np.int16,np.int16,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32)

if limmag != limbright and photoz == 'bpzeazy':
        outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9,outlist50eazylimbright_0,outlist50eazylimbright_1,outlist50eazylimbright_2,outlist50eazylimbright_3,outlist50eazylimbright_4,outlist50eazylimbright_5,outlist50eazylimbright_6,outlist50eazylimbright_7,outlist50eazylimbright_8,outlist50eazylimbright_9,outlist50bpzlimmag_0,outlist50bpzlimmag_1,outlist50bpzlimmag_2,outlist50bpzlimmag_3,outlist50bpzlimmag_4,outlist50bpzlimmag_5,outlist50bpzlimmag_6,outlist50bpzlimmag_7,outlist50bpzlimmag_8,outlist50bpzlimmag_9,outlist50eazylimmag_0,outlist50eazylimmag_1,outlist50eazylimmag_2,outlist50eazylimmag_3,outlist50eazylimmag_4,outlist50eazylimmag_5,outlist50eazylimmag_6,outlist50eazylimmag_7,outlist50eazylimmag_8,outlist50eazylimmag_9 = outputfunc(0.5,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limmageazy,lens_gal_limbrighteazy,lens_zweight_limmageazy,lens_zweight_limbrighteazy,lens_mass_limmageazy,lens_mass_limbrighteazy,lens_mass2_limmageazy,lens_mass2_limbrighteazy,lens_mass3_limmageazy,lens_mass3_limbrighteazy,lens_oneoverr_limmageazy,lens_oneoverr_limbrighteazy,lens_zoverr_limmageazy,lens_zoverr_limbrighteazy,lens_massoverr_limmageazy,lens_massoverr_limbrighteazy,lens_mass2overr_limmageazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limmageazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limmageazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limmageazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limmageazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limmageazy,lens_mass3overrms_limbrighteazy,lens_flexion_limmageazy,lens_flexion_limbrighteazy,lens_tidal_limmageazy,lens_tidal_limbrighteazy,lens_convergence_limmageazy,lens_convergence_limbrighteazy,lens_convergencehalo_limmageazy,lens_convergencehalo_limbrighteazy)
        outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9,outlist75eazylimbright_0,outlist75eazylimbright_1,outlist75eazylimbright_2,outlist75eazylimbright_3,outlist75eazylimbright_4,outlist75eazylimbright_5,outlist75eazylimbright_6,outlist75eazylimbright_7,outlist75eazylimbright_8,outlist75eazylimbright_9,outlist75bpzlimmag_0,outlist75bpzlimmag_1,outlist75bpzlimmag_2,outlist75bpzlimmag_3,outlist75bpzlimmag_4,outlist75bpzlimmag_5,outlist75bpzlimmag_6,outlist75bpzlimmag_7,outlist75bpzlimmag_8,outlist75bpzlimmag_9,outlist75eazylimmag_0,outlist75eazylimmag_1,outlist75eazylimmag_2,outlist75eazylimmag_3,outlist75eazylimmag_4,outlist75eazylimmag_5,outlist75eazylimmag_6,outlist75eazylimmag_7,outlist75eazylimmag_8,outlist75eazylimmag_9 = outputfunc(0.75,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limmageazy,lens_gal_limbrighteazy,lens_zweight_limmageazy,lens_zweight_limbrighteazy,lens_mass_limmageazy,lens_mass_limbrighteazy,lens_mass2_limmageazy,lens_mass2_limbrighteazy,lens_mass3_limmageazy,lens_mass3_limbrighteazy,lens_oneoverr_limmageazy,lens_oneoverr_limbrighteazy,lens_zoverr_limmageazy,lens_zoverr_limbrighteazy,lens_massoverr_limmageazy,lens_massoverr_limbrighteazy,lens_mass2overr_limmageazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limmageazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limmageazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limmageazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limmageazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limmageazy,lens_mass3overrms_limbrighteazy,lens_flexion_limmageazy,lens_flexion_limbrighteazy,lens_tidal_limmageazy,lens_tidal_limbrighteazy,lens_convergence_limmageazy,lens_convergence_limbrighteazy,lens_convergencehalo_limmageazy,lens_convergencehalo_limbrighteazy)

if limmag != limbright and photoz == 'bpz':
        outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9,outlist50bpzlimmag_0,outlist50bpzlimmag_1,outlist50bpzlimmag_2,outlist50bpzlimmag_3,outlist50bpzlimmag_4,outlist50bpzlimmag_5,outlist50bpzlimmag_6,outlist50bpzlimmag_7,outlist50bpzlimmag_8,outlist50bpzlimmag_9 = outputfunc(0.5,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)
        outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9,outlist75bpzlimmag_0,outlist75bpzlimmag_1,outlist75bpzlimmag_2,outlist75bpzlimmag_3,outlist75bpzlimmag_4,outlist75bpzlimmag_5,outlist75bpzlimmag_6,outlist75bpzlimmag_7,outlist75bpzlimmag_8,outlist75bpzlimmag_9 = outputfunc(0.75,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)

if limmag == limbright and photoz == 'bpzeazy':
        outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9,outlist50eazylimbright_0,outlist50eazylimbright_1,outlist50eazylimbright_2,outlist50eazylimbright_3,outlist50eazylimbright_4,outlist50eazylimbright_5,outlist50eazylimbright_6,outlist50eazylimbright_7,outlist50eazylimbright_8,outlist50eazylimbright_9 = outputfunc(0.50,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limbrighteazy,lens_zweight_limbrighteazy,lens_mass_limbrighteazy,lens_mass2_limbrighteazy,lens_mass3_limbrighteazy,lens_oneoverr_limbrighteazy,lens_zoverr_limbrighteazy,lens_massoverr_limbrighteazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limbrighteazy,lens_flexion_limbrighteazy,lens_tidal_limbrighteazy,lens_convergence_limbrighteazy,lens_convergencehalo_limbrighteazy)
        outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9,outlist75eazylimbright_0,outlist75eazylimbright_1,outlist75eazylimbright_2,outlist75eazylimbright_3,outlist75eazylimbright_4,outlist75eazylimbright_5,outlist75eazylimbright_6,outlist75eazylimbright_7,outlist75eazylimbright_8,outlist75eazylimbright_9 = outputfunc(0.75,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright,lens_gal_limbrighteazy,lens_zweight_limbrighteazy,lens_mass_limbrighteazy,lens_mass2_limbrighteazy,lens_mass3_limbrighteazy,lens_oneoverr_limbrighteazy,lens_zoverr_limbrighteazy,lens_massoverr_limbrighteazy,lens_mass2overr_limbrighteazy,lens_mass3overr_limbrighteazy,lens_mass2rms_limbrighteazy,lens_mass3rms_limbrighteazy,lens_mass2overrms_limbrighteazy,lens_mass3overrms_limbrighteazy,lens_flexion_limbrighteazy,lens_tidal_limbrighteazy,lens_convergence_limbrighteazy,lens_convergencehalo_limbrighteazy)

if limmag == limbright and photoz == 'bpz':
        outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9 = outputfunc(0.50,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)
        outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9 = outputfunc(0.75,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)

if photoz == 'bpzeazy' or photoz == 'bpz':
        t = table.Table(outlist50bpzlimbright_0.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_0,overwrite=True)
        t = table.Table(outlist75bpzlimbright_0.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_0,overwrite=True)
        t = table.Table(outlist50bpzlimbright_1.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_1,overwrite=True)
        t = table.Table(outlist75bpzlimbright_1.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_1,overwrite=True)
        t = table.Table(outlist50bpzlimbright_2.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_2,overwrite=True)
        t = table.Table(outlist75bpzlimbright_2.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_2,overwrite=True)
        t = table.Table(outlist50bpzlimbright_3.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_3,overwrite=True)
        t = table.Table(outlist75bpzlimbright_3.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_3,overwrite=True)
        t = table.Table(outlist50bpzlimbright_4.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_4,overwrite=True)
        t = table.Table(outlist75bpzlimbright_4.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_4,overwrite=True)
        t = table.Table(outlist50bpzlimbright_5.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_5,overwrite=True)
        t = table.Table(outlist75bpzlimbright_5.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_5,overwrite=True)
        t = table.Table(outlist50bpzlimbright_6.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_6,overwrite=True)
        t = table.Table(outlist75bpzlimbright_6.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_6,overwrite=True)
        t = table.Table(outlist50bpzlimbright_7.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_7,overwrite=True)
        t = table.Table(outlist75bpzlimbright_7.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_7,overwrite=True)
        t = table.Table(outlist50bpzlimbright_8.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_8,overwrite=True)
        t = table.Table(outlist75bpzlimbright_8.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_8,overwrite=True)
        t = table.Table(outlist50bpzlimbright_9.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_9,overwrite=True)
        t = table.Table(outlist75bpzlimbright_9.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_9,overwrite=True)
if photoz == 'bpzeazy':
        t = table.Table(outlist50eazylimbright_0.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_0,overwrite=True)
        t = table.Table(outlist75eazylimbright_0.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_0,overwrite=True)
        t = table.Table(outlist50eazylimbright_1.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_1,overwrite=True)
        t = table.Table(outlist75eazylimbright_1.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_1,overwrite=True)
        t = table.Table(outlist50eazylimbright_2.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_2,overwrite=True)
        t = table.Table(outlist75eazylimbright_2.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_2,overwrite=True)
        t = table.Table(outlist50eazylimbright_3.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_3,overwrite=True)
        t = table.Table(outlist75eazylimbright_3.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_3,overwrite=True)
        t = table.Table(outlist50eazylimbright_4.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_4,overwrite=True)
        t = table.Table(outlist75eazylimbright_4.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_4,overwrite=True)
        t = table.Table(outlist50eazylimbright_5.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_5,overwrite=True)
        t = table.Table(outlist75eazylimbright_5.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_5,overwrite=True)
        t = table.Table(outlist50eazylimbright_6.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_6,overwrite=True)
        t = table.Table(outlist75eazylimbright_6.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_6,overwrite=True)
        t = table.Table(outlist50eazylimbright_7.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_7,overwrite=True)
        t = table.Table(outlist75eazylimbright_7.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_7,overwrite=True)
        t = table.Table(outlist50eazylimbright_8.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_8,overwrite=True)
        t = table.Table(outlist75eazylimbright_8.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_8,overwrite=True)
        t = table.Table(outlist50eazylimbright_9.T, names=names, dtype=dtype); t.write(fout50_eazylimbright_9,overwrite=True)
        t = table.Table(outlist75eazylimbright_9.T, names=names, dtype=dtype); t.write(fout75_eazylimbright_9,overwrite=True)
if (limmag != limbright and photoz == 'bpz') or (limmag != limbright and photoz == 'bpzeazy'):
        t = table.Table(outlist50bpzlimmag_0.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_0,overwrite=True)
        t = table.Table(outlist75bpzlimmag_0.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_0,overwrite=True)
        t = table.Table(outlist50bpzlimmag_1.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_1,overwrite=True)
        t = table.Table(outlist75bpzlimmag_1.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_1,overwrite=True)
        t = table.Table(outlist50bpzlimmag_2.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_2,overwrite=True)
        t = table.Table(outlist75bpzlimmag_2.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_2,overwrite=True)
        t = table.Table(outlist50bpzlimmag_3.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_3,overwrite=True)
        t = table.Table(outlist75bpzlimmag_3.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_3,overwrite=True)
        t = table.Table(outlist50bpzlimmag_4.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_4,overwrite=True)
        t = table.Table(outlist75bpzlimmag_4.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_4,overwrite=True)
        t = table.Table(outlist50bpzlimmag_5.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_5,overwrite=True)
        t = table.Table(outlist75bpzlimmag_5.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_5,overwrite=True)
        t = table.Table(outlist50bpzlimmag_6.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_6,overwrite=True)
        t = table.Table(outlist75bpzlimmag_6.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_6,overwrite=True)
        t = table.Table(outlist50bpzlimmag_7.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_7,overwrite=True)
        t = table.Table(outlist75bpzlimmag_7.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_7,overwrite=True)
        t = table.Table(outlist50bpzlimmag_8.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_8,overwrite=True)
        t = table.Table(outlist75bpzlimmag_8.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_8,overwrite=True)
        t = table.Table(outlist50bpzlimmag_9.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_9,overwrite=True)
        t = table.Table(outlist75bpzlimmag_9.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_9,overwrite=True)
if limmag != limbright and photoz == 'bpzeazy':
        t = table.Table(outlist50eazylimmag_0.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_0,overwrite=True)
        t = table.Table(outlist75eazylimmag_0.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_0,overwrite=True)
        t = table.Table(outlist50eazylimmag_1.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_1,overwrite=True)
        t = table.Table(outlist75eazylimmag_1.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_1,overwrite=True)
        t = table.Table(outlist50eazylimmag_2.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_2,overwrite=True)
        t = table.Table(outlist75eazylimmag_2.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_2,overwrite=True)
        t = table.Table(outlist50eazylimmag_3.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_3,overwrite=True)
        t = table.Table(outlist75eazylimmag_3.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_3,overwrite=True)
        t = table.Table(outlist50eazylimmag_4.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_4,overwrite=True)
        t = table.Table(outlist75eazylimmag_4.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_4,overwrite=True)
        t = table.Table(outlist50eazylimmag_5.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_5,overwrite=True)
        t = table.Table(outlist75eazylimmag_5.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_5,overwrite=True)
        t = table.Table(outlist50eazylimmag_6.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_6,overwrite=True)
        t = table.Table(outlist75eazylimmag_6.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_6,overwrite=True)
        t = table.Table(outlist50eazylimmag_7.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_7,overwrite=True)
        t = table.Table(outlist75eazylimmag_7.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_7,overwrite=True)
        t = table.Table(outlist50eazylimmag_8.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_8,overwrite=True)
        t = table.Table(outlist75eazylimmag_8.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_8,overwrite=True)
        t = table.Table(outlist50eazylimmag_9.T, names=names, dtype=dtype); t.write(fout50_eazylimmag_9,overwrite=True)
        t = table.Table(outlist75eazylimmag_9.T, names=names, dtype=dtype); t.write(fout75_eazylimmag_9,overwrite=True)

print("Writing output completed in %0.1f seconds" % ((time.time() - start_write)))
print("Total time for %s: --- %0.2f seconds ---" % ('%s_wghtratios_75_%s_%s_limmag_eazy_%s_%s_%s_zgap%s_%s%s' % (fieldID[-26:-20],mskname,lensID,det,irac,type,zinf,zsup,suffix), (time.time() - start_time)))

print 'Done!'
