# CE Rusu, Feb 14 2018
# This code makes use of the CFHTLens *galphotmstar.cat files and the lens photometric+Mstar+photoz catalogue; it computes weighted ratios (lens/field) with proper masking, for various radii, limiting mag, number of samples and classification scheme.
# run as: python /Users/cerusu/GITHUB/zMstarPDF/python/catalogue_utilities/weightinguniversal_overlap_sampling_nobetasingleband.py PG1115 /Volumes/LaCieDavis/CFHTcatalogues/W1m0m0_r24galphot.fits /Volumes/LaCieDavis/CFHTLenSmasks/W1m0m0_izrgu_finalmask_mosaic.fits /Users/cerusu/Dropbox/Davis_work/code /Volumes/LaCieSubaru/weightedcounts/PG1115 45 5 meds removehandpicked
# the scripts to run en masse are in /Users/cerusu/GITHUB/zMstarPDF/python/scripts/DESKTOP/
# the code is optimized for speed, but may be memory intensive because it stores the input catalogue in memory
# definitions:
#field refers to one of the 171 CHFT fields
#lens refers to the name of the H0licow lens
#field mask refers to the mask .fits file of the corresponding CHFT field
#cell refers to a rectangular surface of length on a side 2*radius
#radius is the aperture radius around the actual lens object around which the environmant is explored
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
type = str(sys.argv[8])
remove = str(sys.argv[9])

suffix = '' # added to the name of the output catalogues
samples = 20 + 1
#samples = 4 + 1

#print("Arguments: \n Lens: %s \n Field file: %s \n Field mask file: %s \n Limiting radius: %s \n" % (lensID, fieldID, fmask, radius))

def lensinit(lensbpz_masked_,lens_gal_,lens_oneoverr_):
    if (lensbpz_masked_.shape[1] != 0):
        lens_gal_[k][l][i][j][n] = lensbpz_masked_.shape[1]
        #print lens_gal_[k][l][i][j][n] # test to check if the function actually returns the result globally
        lens_oneoverr = 1.0/lensbpz_masked_[sep]
        if type == "meds":
            lens_oneoverr_[k][l][i][j][n] = np.median(lens_oneoverr) * lens_gal_[k][l][i][j][n]
        if type == "sum":
            lens_oneoverr_[k][l][i][j][n] = np.sum(lens_oneoverr)
    return lens_gal_,lens_oneoverr_

def fieldinit(field_masked_,w_gal_,field_gal_,field_oneoverr_):
    if (w_gal_ != 0):
        field_gal_[k][l][i][j] = w_gal_
        field_oneoverr = 1.0/field_masked_[cell_sep]
        if type == "meds":
            field_oneoverr_[k][l][i][j] = np.median(field_oneoverr) * w_gal_
        if type == "sum":
            field_oneoverr_[k][l][i][j] = np.sum(field_oneoverr)
    return field_gal_,field_oneoverr_

start_time = time.time()

do50 = False # whether or not to consider fields with less than 50% of the surface covered by masks, or just less than 75%
if lensID == "PG1115":
    limmag = 23
    limbright = 23
    brightmag = 15.3
    pixnr = 1011
    pixlens = 0.237388724 * u.arcsec
    boosterr = 0.10

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

######################################################################
# WORKING ON THE LENS CATALOGUE
######################################################################

print "Initializing the lens catalogue..."

if lensID == "PG1115":
    center_lens = SkyCoord('11:18:16.900 +07:45:59.00', frame='fk5', unit=(u.hourangle, u.deg))

lensbpz = np.loadtxt('%s.cat' % rootlenscat, unpack=True)
if np.shape(lensbpz)[0] != 8 + 3 * samples:
    print("Number of declared samples and the ones inside the input file do not match!!!")
    sys.exit()

msk_lens = fits.open('%s/msk%sarcsecrad%sarcsecgap.fits' % (rootlenscat,radius,inner))
# defining the columns:
x_lens = 0
y_lens = 1
RA_lens = 2
DEC_lens = 3
flux_radius = 8
r_lens = 11
r_err_lens = 12
# soon to be defined:
sep = 16 # distance in arcsec between the lens center and the galaxies

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
    #lenscat = np.delete(lenscat,np.where(lenscat[classify] < 0),axis=1) # removes all stars from the catalogue
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
if limmag != limbright: lens_oneoverr_limmagbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))
lens_oneoverr_limbrightbpz = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side,samples))

print("Initialization of lens catalogue was completed in %0.1f seconds" % (time.time() - start_time))

######################################################################
# WORKING ON THE FIELD CATALOGUE
######################################################################

start_timefield = time.time()
print "Reading ", fieldID, "..."
['MAG_r','MAGERR_r','id','ALPHA_J2000','DELTA_J2000','Flag']
field_ = Table.read('%s' % fieldID)
field = np.c_[field_['ALPHA_J2000'],field_['DELTA_J2000'],field_['MAG_r']].T
#print np.shape(field)
del field_
RA_field = 0
DEC_field = 1
r_field = 2
# will define shortly:
cell_xpix = -3
cell_ypix = -2
cell_sep = -1

#field = np.delete(field,np.where(field[z_field] > z_s),axis=1) # eliminate objects at redshifts higher than the source
#print 'shape', np.shape(field)
#field = np.delete(field,np.where((field[z_field] >= zinf) & (field[z_field] <= zsup)),axis=1) # eliminate objects corresponding to the redshift slice
#print 'shape', np.shape(field)
#field[r_field][field[r_field] < 0] = field[y_field][field[r_field] < 0] # for objects with y mags, use those
field = np.delete(field,np.where(field[r_field] < brightmag),axis = 1) # eliminate objects brighter than the upper brightness limit
#field[mass_BEST_field][field[mass_BEST_field] < 0] = 9.0 # fix the very rare unphysical masses
#field[mass_MED_field][field[mass_MED_field] < 0] = field[mass_BEST_field][field[mass_MED_field] < 0] # fix the more common unphysical med masses with best masses

fieldpix = np.zeros((3,field.shape[1])) # this will contain the fieldpix[1], fieldpix[0] and the separation in arcsec from the cell center
field = np.c_['0',field,fieldpix.reshape(3,fieldpix.shape[1])] # adding fieldpix as extension to the catalogue
fieldpix = worldfield.wcs_world2pix(field[RA_field], field[DEC_field], 1) # find pixel coordinates of each field galaxy, relative to the fmask .fits file (therefore with CFHTLenS pixel size); first column is the physical image x axis
field[cell_xpix] = fieldpix[1] # the physical image y axis, in agreement with the lens section above
field[cell_ypix] = fieldpix[0]

# declarations
if limmag != limbright: field_gal_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_gal_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
if limmag != limbright: field_oneoverr_limmag = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))
field_oneoverr_limbright = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

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

                    for n in range(samples):
                        lensbpz_maskedlimmag = np.copy(lensbpz_masked)

                        if n == 0:
                            if limmag != limbright:
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where(lensbpz_maskedlimmag[r_lens] > limmag),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimmag,np.where(lensbpz_maskedlimmag[r_lens] > limbright),axis=1)
                        else:
                            for o in range(lensbpz_maskedlimmag.shape[1]):
                                lensbpz_maskedlimmag[r_lens][o] = np.random.normal(lensbpz_maskedlimmag[r_lens][o], np.max([lensbpz_maskedlimmag[i_err_lens][o],0.005]),1)[0]
                            if limmag != limbright:
                                lensbpz_maskedlimmag = np.delete(lensbpz_maskedlimmag,np.where(lensbpz_maskedlimmag[r_lens] > limmag),axis=1)
                            lensbpz_maskedlimbright = np.delete(lensbpz_maskedlimmag,np.where(lensbpz_maskedlimmag[r_lens] > limbright),axis=1)
                        if limmag != limbright:
                            lens_gal_limmagbpz,lens_oneoverr_limmagbpz = lensinit(lensbpz_maskedlimmag,lens_gal_limmagbpz,lens_oneoverr_limmagbpz)
                        #print "d ",lens_gal_limmagbpz[k][l][i][j][n]# test to check if the function actually returns the result globally
                        lens_gal_limbrightbpz,lens_oneoverr_limbrightbpz = lensinit(lensbpz_maskedlimbright,lens_gal_limbrightbpz)
                        #print "d ",lens_gal_limbrightbpz[k][l][i][j][n]# test to check if the function actually returns the result globally

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
                    field_masked_limmag = np.delete(field_masked_limmag,np.where(field_masked_limmag[r_field] > limmag),axis=1)
                    field_masked_limbright = np.delete(field_masked_limmag,np.where(field_masked_limmag[r_field] > limbright),axis=1)
                    w_gal_limmag = np.shape(field_masked_limmag)[1]
                    w_gal_limbright = np.shape(field_masked_limbright)[1]
                    #mmm = np.copy(field_convergencehalo_limbright)
                    if limmag != limbright:
                       field_gal_limmag,field_oneoverr_limmag = fieldinit(field_masked_limmag,w_gal_limmag,field_gal_limmag,field_oneoverr_limmag)
                    field_gal_limbright,field_oneoverr_limbright = fieldinit(field_masked_limbright,w_gal_limbright,field_gal_limbright,field_oneoverr_limbright)
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
fcount = open('%s/%s_wghtratios_%s_%s_%s_zgap%s_%s_count.cat' % (output,fieldID[-26:-20],mskname,lensID,det,irac,type,suffix),'w') # [-25:-4] corresponds to strings of the form W1m0m0_limmaggalphotmstar
fcount.write(count)
fcount.close()

if do50 == True: fout50_bpzlimbright_0 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
fout75_bpzlimbright_0 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
if do50 == True: fout50_bpzlimbright_1 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
fout75_bpzlimbright_1 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
if do50 == True: fout50_bpzlimbright_2 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
fout75_bpzlimbright_2 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
if do50 == True: fout50_bpzlimbright_3 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
fout75_bpzlimbright_3 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
if do50 == True: fout50_bpzlimbright_4 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
fout75_bpzlimbright_4 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
if samples > 5:
    if do50 == True: fout50_bpzlimbright_5 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_5 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_6 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_6 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_7 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_7 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_8 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_8 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_9 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_9 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_10 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_10.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_10 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_10.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_11 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_11.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_11 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_11.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_12 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_12.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_12 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_12.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_13 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_13.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_13 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_13.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_14 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_14.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_14 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_14.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_15 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_15.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_15 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_15.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_16 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_16.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_16 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_16.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_17 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_17.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_17 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_17.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_18 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_18.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_18 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_18.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_19 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_19.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_19 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_19.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimbright_20 = '%s/%s_wghtratios_50_%s_%s_%s_bpz_%s_zgap%s_20.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimbright_20 = '%s/%s_wghtratios_75_%s_%s_%s_bpz_%s_zgap%s_20.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
#os.system('rm -f %s' % fout50_0) # '-f' ignores non-existent files
if limmag != brightmag:
    if do50 == True: fout50_bpzlimmag_0 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimmag_0 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_0.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimmag_1 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimmag_1 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_1.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimmag_2 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimmag_2 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_2.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimmag_3 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimmag_3 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_3.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if do50 == True: fout50_bpzlimmag_4 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    fout75_bpzlimmag_4 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_4.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
    if samples > 5:
        if do50 == True: fout50_bpzlimmag_5 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_5 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_5.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_6 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_6 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_6.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_7 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_7 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_7.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_8 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_8 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_8.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_9 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_9 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_9.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_10 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_10.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_10 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_10.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_11 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_11.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_11 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_11.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_12 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_12.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_12 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_12.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_13 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_13.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_13 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_13.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_14 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_14.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_14 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_14.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_15 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_15.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_15 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_15.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_16 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_16.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_16 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_16.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_17 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_17.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_17 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_17.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_18 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_18.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_18 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_18.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_19 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_19.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_19 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_19.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        if do50 == True: fout50_bpzlimmag_20 = '%s/%s_wghtratios_50_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_20.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)
        fout75_bpzlimmag_20 = '%s/%s_wghtratios_75_%s_%s_limmag_bpz_%s_%s_%s_zgap%s_%s%s_20.fits' % (output,fieldID[-26:-20],mskname,lensID,limbright,type,suffix)

def outputfunc(*argv):
    start0 = True; start1 = True; start2 = True; start3 = True; start4 = True; start5 = True; start6 = True; start7 = True; start8 = True; start9 = True; start10 = True; start11 = True; start12 = True; start13 = True; start14 = True; start15 = True; start16 = True; start17 = True; start18 = True; start19 = True; start20 = True
    if len(argv) == 9:
        '''limbright limmag bpz'''
        # frac,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright
        frac_ = argv[0]; lens_gal_limmagbpz_ = argv[1]; field_gal_limmag_ = argv[2]; lens_gal_limbrightbpz_ = argv[3]; field_gal_limbright_ = argv[4]; lens_oneoverr_limmagbpz_ = argv[5]; field_oneoverr_limmag_ = argv[6]; lens_oneoverr_limbrightbpz_ = argv[7]; field_oneoverr_limbright_ = argv[8]
        outbpzlimbright = np.zeros(22)
        outbpzlimmag = np.zeros(22)
    if len(argv) == 5:
        '''limbright bpz'''
        # frac,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright
        frac_ = argv[0]; lens_gal_limbrightbpz_ = argv[1]; field_gal_limbright_ = argv[2]; lens_oneoverr_limbrightbpz_ = argv[3]; field_oneoverr_limbright_ = argv[4]
        outbpzlimbright = np.zeros(22)

    for k in range(overlap):
        for l in range(overlap):
            for i in range(cells_on_a_side):
                for j in range(cells_on_a_side):
                    for n in range(samples):
                        condition = False
                        if len(argv) == 9:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limmagbpz_[k][l][i][j]) != 0) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (field_gal_limmag_[k][l][i][j] != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if len(argv) == 5:
                            if (unmaskedcell[k][l][i][j] >= frac_) & (np.min(lens_gal_limbrightbpz_[k][l][i][j]) != 0) & (field_gal_limbright_[k][l][i][j] != 0): condition = True
                        if condition == True:
                            if len(argv) == 9 or len(argv) == 5:
                                outbpzlimbright[0] = np.int16(k)
                                outbpzlimbright[1] = np.int16(l)
                                outbpzlimbright[2] = np.int16(i)
                                outbpzlimbright[3] = np.int16(j)
                                outbpzlimbright[4] = np.float32(1.0*lens_gal_limbrightbpz_[k][l][i][j][n]/field_gal_limbright_[k][l][i][j])
                                outbpzlimbright[5] = np.float32(1.0*lens_oneoverr_limbrightbpz_[k][l][i][j][n]/field_oneoverr_limbright_[k][l][i][j])
                            if len(argv) == 5:
                                outbpzlimmag[0] = np.int16(k)
                                outbpzlimmag[1] = np.int16(l)
                                outbpzlimmag[2] = np.int16(i)
                                outbpzlimmag[3] = np.int16(j)
                                outbpzlimmag[4] = np.float32(1.0*lens_gal_limmagbpz_[k][l][i][j][n]/field_gal_limmag_[k][l][i][j])
                                outbpzlimmag[5] = np.float32(1.0*lens_oneoverr_limmagbpz_[k][l][i][j][n]/field_oneoverr_limmag_[k][l][i][j])
                            if n == 0:
                                if start0 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_0 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_0 = outbpzlimmag
                                    start0 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_0 = np.c_[outlistbpzlimbright_0,outbpzlimbright] # normally this line should start with else, but for some reason that makes it lose the first row and duplicate the second; this way it will duplicate the first line, and at the end I will mask its first instance
                                if len(argv) == 5: outlistbpzlimmag_0 = np.c_[outlistbpzlimmag_0,outbpzlimmag]
                            if n == 1:
                                if start1 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_1 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_1 = outbpzlimmag
                                    start1 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_1 = np.c_[outlistbpzlimbright_1,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_1 = np.c_[outlistbpzlimmag_1,outbpzlimmag]
                            if n == 2:
                                if start2 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_2 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_2 = outbpzlimmag
                                    start2 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_2 = np.c_[outlistbpzlimbright_2,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_2 = np.c_[outlistbpzlimmag_2,outbpzlimmag]
                            if n == 3:
                                if start3 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_3 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_3 = outbpzlimmag
                                    start3 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_3 = np.c_[outlistbpzlimbright_3,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_3 = np.c_[outlistbpzlimmag_3,outbpzlimmag]
                            if n == 4:
                                if start4 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_4 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_4 = outbpzlimmag
                                    start4 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_4 = np.c_[outlistbpzlimbright_4,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_4 = np.c_[outlistbpzlimmag_4,outbpzlimmag]
                            if n == 5:
                                if start5 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_5 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_5 = outbpzlimmag
                                    start5 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_5 = np.c_[outlistbpzlimbright_5,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_5 = np.c_[outlistbpzlimmag_5,outbpzlimmag]
                            if n == 6:
                                if start6 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_6 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_6 = outbpzlimmag
                                    start6 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_6 = np.c_[outlistbpzlimbright_6,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_6 = np.c_[outlistbpzlimmag_6,outbpzlimmag]
                            if n == 7:
                                if start7 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_7 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_7 = outbpzlimmag
                                    start7 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_7 = np.c_[outlistbpzlimbright_7,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_7 = np.c_[outlistbpzlimmag_7,outbpzlimmag]
                            if n == 8:
                                if start8 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_8 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_8 = outbpzlimmag
                                    start8 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_8 = np.c_[outlistbpzlimbright_8,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_8 = np.c_[outlistbpzlimmag_8,outbpzlimmag]
                            if n == 9:
                                if start9 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_9 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_9 = outbpzlimmag
                                    start9 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_9 = np.c_[outlistbpzlimbright_9,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_9 = np.c_[outlistbpzlimmag_9,outbpzlimmag]
                            if n == 10:
                                if start10 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_10 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_10 = outbpzlimmag
                                    start10 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_10 = np.c_[outlistbpzlimbright_10,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_10 = np.c_[outlistbpzlimmag_10,outbpzlimmag]
                            if n == 11:
                                if start11 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_11 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_11 = outbpzlimmag
                                    start11 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_11 = np.c_[outlistbpzlimbright_11,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_11 = np.c_[outlistbpzlimmag_11,outbpzlimmag]
                            if n == 12:
                                if start12 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_12 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_12 = outbpzlimmag
                                    start12 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_12 = np.c_[outlistbpzlimbright_12,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_12 = np.c_[outlistbpzlimmag_12,outbpzlimmag]
                            if n == 13:
                                if start13 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_13 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_13 = outbpzlimmag
                                    start13 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_13 = np.c_[outlistbpzlimbright_13,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_13 = np.c_[outlistbpzlimmag_13,outbpzlimmag]
                            if n == 14:
                                if start14 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_14 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_14 = outbpzlimmag
                                    start14 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_14 = np.c_[outlistbpzlimbright_14,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_14 = np.c_[outlistbpzlimmag_14,outbpzlimmag]
                            if n == 15:
                                if start15 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_15 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_15 = outbpzlimmag
                                    start15 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_15 = np.c_[outlistbpzlimbright_15,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_15 = np.c_[outlistbpzlimmag_15,outbpzlimmag]
                            if n == 16:
                                if start16 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_16 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_16 = outbpzlimmag
                                    start16 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_16 = np.c_[outlistbpzlimbright_16,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_16 = np.c_[outlistbpzlimmag_16,outbpzlimmag]
                            if n == 17:
                                if start17 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_17 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_17 = outbpzlimmag
                                    start17 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_17 = np.c_[outlistbpzlimbright_17,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_17 = np.c_[outlistbpzlimmag_17,outbpzlimmag]
                            if n == 18:
                                if start18 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_18 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_18 = outbpzlimmag
                                    start18 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_18 = np.c_[outlistbpzlimbright_18,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_18 = np.c_[outlistbpzlimmag_18,outbpzlimmag]
                            if n == 19:
                                if start19 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_19 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag_19 = outbpzlimmag
                                    start19 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_19 = np.c_[outlistbpzlimbright_19,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_19 = np.c_[outlistbpzlimmag_19,outbpzlimmag]
                            if n == 20:
                                if start20 == True:
                                    if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_20 = outbpzlimbright
                                    if len(argv) == 5: outlistbpzlimmag20 = outbpzlimmag
                                    start20 = False
                                if len(argv) == 9 or len(argv) == 5: outlistbpzlimbright_20 = np.c_[outlistbpzlimbright_20,outbpzlimbright]
                                if len(argv) == 5: outlistbpzlimmag_20 = np.c_[outlistbpzlimmag_20,outbpzlimmag]

    if len(argv) == 9:
        if samples > 5: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:],outlistbpzlimbright_10[:,1:],outlistbpzlimbright_11[:,1:],outlistbpzlimbright_12[:,1:],outlistbpzlimbright_13[:,1:],outlistbpzlimbright_14[:,1:],outlistbpzlimbright_15[:,1:],outlistbpzlimbright_16[:,1:],outlistbpzlimbright_17[:,1:],outlistbpzlimbright_18[:,1:],outlistbpzlimbright_19[:,1:],outlistbpzlimbright_20[:,1:],outlistbpzlimmag_0[:,1:],outlistbpzlimmag_1[:,1:],outlistbpzlimmag_2[:,1:],outlistbpzlimmag_3[:,1:],outlistbpzlimmag_4[:,1:],outlistbpzlimmag_5[:,1:],outlistbpzlimmag_6[:,1:],outlistbpzlimmag_7[:,1:],outlistbpzlimmag_8[:,1:],outlistbpzlimmag_9[:,1:],outlistbpzlimmag_10[:,1:],outlistbpzlimmag_11[:,1:],outlistbpzlimmag_12[:,1:],outlistbpzlimmag_13[:,1:],outlistbpzlimmag_14[:,1:],outlistbpzlimmag_15[:,1:],outlistbpzlimmag_16[:,1:],outlistbpzlimmag_17[:,1:],outlistbpzlimmag_18[:,1:],outlistbpzlimmag_19[:,1:],outlistbpzlimmag_20[:,1:]
        else: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimmag_0[:,1:],outlistbpzlimmag_1[:,1:],outlistbpzlimmag_2[:,1:],outlistbpzlimmag_3[:,1:],outlistbpzlimmag_4[:,1:]
    if len(argv) == 5:
        if samples > 5: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:],outlistbpzlimbright_5[:,1:],outlistbpzlimbright_6[:,1:],outlistbpzlimbright_7[:,1:],outlistbpzlimbright_8[:,1:],outlistbpzlimbright_9[:,1:],outlistbpzlimbright_10[:,1:],outlistbpzlimbright_11[:,1:],outlistbpzlimbright_12[:,1:],outlistbpzlimbright_13[:,1:],outlistbpzlimbright_14[:,1:],outlistbpzlimbright_15[:,1:],outlistbpzlimbright_16[:,1:],outlistbpzlimbright_17[:,1:],outlistbpzlimbright_18[:,1:],outlistbpzlimbright_19[:,1:],outlistbpzlimbright_20[:,1:]
        else: return outlistbpzlimbright_0[:,1:],outlistbpzlimbright_1[:,1:],outlistbpzlimbright_2[:,1:],outlistbpzlimbright_3[:,1:],outlistbpzlimbright_4[:,1:]

names = ('1_overlap_x','2_overlap_y','3_cell_x','4_cell_y','5_lens_gal','6_lens_oneoverr')
dtype=(np.int16,np.int16,np.int16,np.int16,np.float32,np.float32)

if limmag != limbright:
        if samples > 5:
            if do50 == True: outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9,outlist50bpzlimbright_10,outlist50bpzlimbright_11,outlist50bpzlimbright_12,outlist50bpzlimbright_13,outlist50bpzlimbright_14,outlist50bpzlimbright_15,outlist50bpzlimbright_16,outlist50bpzlimbright_17,outlist50bpzlimbright_18,outlist50bpzlimbright_19,outlist50bpzlimbright_20,outlist50bpzlimmag_0,outlist50bpzlimmag_1,outlist50bpzlimmag_2,outlist50bpzlimmag_3,outlist50bpzlimmag_4,outlist50bpzlimmag_5,outlist50bpzlimmag_6,outlist50bpzlimmag_7,outlist50bpzlimmag_8,outlist50bpzlimmag_9,outlist50bpzlimmag_10,outlist50bpzlimmag_11,outlist50bpzlimmag_12,outlist50bpzlimmag_13,outlist50bpzlimmag_14,outlist50bpzlimmag_15,outlist50bpzlimmag_16,outlist50bpzlimmag_17,outlist50bpzlimmag_18,outlist50bpzlimmag_19,outlist50bpzlimmag_20 = outputfunc(0.5,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)
            outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9,outlist75bpzlimbright_10,outlist75bpzlimbright_11,outlist75bpzlimbright_12,outlist75bpzlimbright_13,outlist75bpzlimbright_14,outlist75bpzlimbright_15,outlist75bpzlimbright_16,outlist75bpzlimbright_17,outlist75bpzlimbright_18,outlist75bpzlimbright_19,outlist75bpzlimbright_20,outlist75bpzlimmag_0,outlist75bpzlimmag_1,outlist75bpzlimmag_2,outlist75bpzlimmag_3,outlist75bpzlimmag_4,outlist75bpzlimmag_5,outlist75bpzlimmag_6,outlist75bpzlimmag_7,outlist75bpzlimmag_8,outlist75bpzlimmag_9,outlist75bpzlimmag_10,outlist75bpzlimmag_11,outlist75bpzlimmag_12,outlist75bpzlimmag_13,outlist75bpzlimmag_14,outlist75bpzlimmag_15,outlist75bpzlimmag_16,outlist75bpzlimmag_17,outlist75bpzlimmag_18,outlist75bpzlimmag_19,outlist75bpzlimmag_20 = outputfunc(0.75,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)
        else:
            if do50 == True: outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimmag_0,outlist50bpzlimmag_1,outlist50bpzlimmag_2,outlist50bpzlimmag_3,outlist50bpzlimmag_4 = outputfunc(0.5,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)
            outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimmag_0,outlist75bpzlimmag_1,outlist75bpzlimmag_2,outlist75bpzlimmag_3,outlist75bpzlimmag_4 = outputfunc(0.75,lens_gal_limmagbpz,field_gal_limmag,lens_gal_limbrightbpz,field_gal_limbright,lens_zweight_limmagbpz,field_zweight_limmag,lens_zweight_limbrightbpz,field_zweight_limbright,lens_mass_limmagbpz,field_mass_limmag,lens_mass_limbrightbpz,field_mass_limbright,lens_mass2_limmagbpz,field_mass2_limmag,lens_mass2_limbrightbpz,field_mass2_limbright,lens_mass3_limmagbpz,field_mass3_limmag,lens_mass3_limbrightbpz,field_mass3_limbright,lens_oneoverr_limmagbpz,field_oneoverr_limmag,lens_oneoverr_limbrightbpz,field_oneoverr_limbright,lens_zoverr_limmagbpz,field_zoverr_limmag,lens_zoverr_limbrightbpz,field_zoverr_limbright,lens_massoverr_limmagbpz,field_massoverr_limmag,lens_massoverr_limbrightbpz,field_massoverr_limbright,lens_mass2overr_limmagbpz,field_mass2overr_limmag,lens_mass2overr_limbrightbpz,field_mass2overr_limbright,lens_mass3overr_limmagbpz,field_mass3overr_limmag,lens_mass3overr_limbrightbpz,field_mass3overr_limbright,lens_mass2rms_limmagbpz,field_mass2rms_limmag,lens_mass2rms_limbrightbpz,field_mass2rms_limbright,lens_mass3rms_limmagbpz,field_mass3rms_limmag,lens_mass3rms_limbrightbpz,field_mass3rms_limbright,lens_mass2overrms_limmagbpz,field_mass2overrms_limmag,lens_mass2overrms_limbrightbpz,field_mass2overrms_limbright,lens_mass3overrms_limmagbpz,field_mass3overrms_limmag,lens_mass3overrms_limbrightbpz,field_mass3overrms_limbright,lens_flexion_limmagbpz,field_flexion_limmag,lens_flexion_limbrightbpz,field_flexion_limbright,lens_tidal_limmagbpz,field_tidal_limmag,lens_tidal_limbrightbpz,field_tidal_limbright,lens_convergence_limmagbpz,field_convergence_limmag,lens_convergence_limbrightbpz,field_convergence_limbright,lens_convergencehalo_limmagbpz,field_convergencehalo_limmag,lens_convergencehalo_limbrightbpz,field_convergencehalo_limbright)

if limmag == limbright:
        if samples > 5:
            if do50 == True: outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4,outlist50bpzlimbright_5,outlist50bpzlimbright_6,outlist50bpzlimbright_7,outlist50bpzlimbright_8,outlist50bpzlimbright_9,outlist50bpzlimbright_10,outlist50bpzlimbright_11,outlist50bpzlimbright_12,outlist50bpzlimbright_13,outlist50bpzlimbright_14,outlist50bpzlimbright_15,outlist50bpzlimbright_16,outlist50bpzlimbright_17,outlist50bpzlimbright_18,outlist50bpzlimbright_19,outlist50bpzlimbright_20 = outputfunc(0.50,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)
            outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4,outlist75bpzlimbright_5,outlist75bpzlimbright_6,outlist75bpzlimbright_7,outlist75bpzlimbright_8,outlist75bpzlimbright_9,outlist75bpzlimbright_10,outlist75bpzlimbright_11,outlist75bpzlimbright_12,outlist75bpzlimbright_13,outlist75bpzlimbright_14,outlist75bpzlimbright_15,outlist75bpzlimbright_16,outlist75bpzlimbright_17,outlist75bpzlimbright_18,outlist75bpzlimbright_19,outlist75bpzlimbright_20 = outputfunc(0.75,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)
        else:
            if do50 == True:
                outlist50bpzlimbright_0,outlist50bpzlimbright_1,outlist50bpzlimbright_2,outlist50bpzlimbright_3,outlist50bpzlimbright_4 = outputfunc(0.50,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)
            outlist75bpzlimbright_0,outlist75bpzlimbright_1,outlist75bpzlimbright_2,outlist75bpzlimbright_3,outlist75bpzlimbright_4 = outputfunc(0.75,lens_gal_limbrightbpz,field_gal_limbright,lens_oneoverr_limbrightbpz,field_oneoverr_limbright)

if do50 == True: t = table.Table(outlist50bpzlimbright_0.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_0,overwrite=True)
t = table.Table(outlist75bpzlimbright_0.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_0,overwrite=True)
if do50 == True: t = table.Table(outlist50bpzlimbright_1.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_1,overwrite=True)
t = table.Table(outlist75bpzlimbright_1.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_1,overwrite=True)
if do50 == True: t = table.Table(outlist50bpzlimbright_2.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_2,overwrite=True)
t = table.Table(outlist75bpzlimbright_2.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_2,overwrite=True)
if do50 == True: t = table.Table(outlist50bpzlimbright_3.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_3,overwrite=True)
t = table.Table(outlist75bpzlimbright_3.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_3,overwrite=True)
if do50 == True: t = table.Table(outlist50bpzlimbright_4.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_4,overwrite=True)
t = table.Table(outlist75bpzlimbright_4.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_4,overwrite=True)
if samples > 5:
    if do50 == True: t = table.Table(outlist50bpzlimbright_5.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_5,overwrite=True)
    t = table.Table(outlist75bpzlimbright_5.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_5,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_6.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_6,overwrite=True)
    t = table.Table(outlist75bpzlimbright_6.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_6,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_7.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_7,overwrite=True)
    t = table.Table(outlist75bpzlimbright_7.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_7,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_8.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_8,overwrite=True)
    t = table.Table(outlist75bpzlimbright_8.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_8,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_9.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_9,overwrite=True)
    t = table.Table(outlist75bpzlimbright_9.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_9,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_10.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_10,overwrite=True)
    t = table.Table(outlist75bpzlimbright_10.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_10,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_11.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_11,overwrite=True)
    t = table.Table(outlist75bpzlimbright_11.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_11,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_12.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_12,overwrite=True)
    t = table.Table(outlist75bpzlimbright_12.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_12,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_13.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_13,overwrite=True)
    t = table.Table(outlist75bpzlimbright_13.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_13,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_14.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_14,overwrite=True)
    t = table.Table(outlist75bpzlimbright_14.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_14,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_15.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_15,overwrite=True)
    t = table.Table(outlist75bpzlimbright_15.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_15,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_16.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_16,overwrite=True)
    t = table.Table(outlist75bpzlimbright_16.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_16,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_17.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_17,overwrite=True)
    t = table.Table(outlist75bpzlimbright_17.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_17,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_18.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_18,overwrite=True)
    t = table.Table(outlist75bpzlimbright_18.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_18,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_19.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_19,overwrite=True)
    t = table.Table(outlist75bpzlimbright_19.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_19,overwrite=True)
    if do50 == True: t = table.Table(outlist50bpzlimbright_20.T, names=names, dtype=dtype); t.write(fout50_bpzlimbright_20,overwrite=True)
    t = table.Table(outlist75bpzlimbright_20.T, names=names, dtype=dtype); t.write(fout75_bpzlimbright_20,overwrite=True)

if limmag != limbright:
        if do50 == True: t = table.Table(outlist50bpzlimmag_0.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_0,overwrite=True)
        t = table.Table(outlist75bpzlimmag_0.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_0,overwrite=True)
        if do50 == True: t = table.Table(outlist50bpzlimmag_1.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_1,overwrite=True)
        t = table.Table(outlist75bpzlimmag_1.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_1,overwrite=True)
        if do50 == True: t = table.Table(outlist50bpzlimmag_2.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_2,overwrite=True)
        t = table.Table(outlist75bpzlimmag_2.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_2,overwrite=True)
        if do50 == True: t = table.Table(outlist50bpzlimmag_3.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_3,overwrite=True)
        t = table.Table(outlist75bpzlimmag_3.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_3,overwrite=True)
        if do50 == True: t = table.Table(outlist50bpzlimmag_4.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_4,overwrite=True)
        t = table.Table(outlist75bpzlimmag_4.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_4,overwrite=True)
        if samples > 5:
            if do50 == True: t = table.Table(outlist50bpzlimmag_5.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_5,overwrite=True)
            t = table.Table(outlist75bpzlimmag_5.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_5,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_6.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_6,overwrite=True)
            t = table.Table(outlist75bpzlimmag_6.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_6,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_7.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_7,overwrite=True)
            t = table.Table(outlist75bpzlimmag_7.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_7,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_8.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_8,overwrite=True)
            t = table.Table(outlist75bpzlimmag_8.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_8,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_9.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_9,overwrite=True)
            t = table.Table(outlist75bpzlimmag_9.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_9,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_10.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_10,overwrite=True)
            t = table.Table(outlist75bpzlimmag_10.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_10,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_11.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_11,overwrite=True)
            t = table.Table(outlist75bpzlimmag_11.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_11,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_12.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_12,overwrite=True)
            t = table.Table(outlist75bpzlimmag_12.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_12,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_13.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_13,overwrite=True)
            t = table.Table(outlist75bpzlimmag_13.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_13,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_14.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_14,overwrite=True)
            t = table.Table(outlist75bpzlimmag_14.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_14,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_15.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_15,overwrite=True)
            t = table.Table(outlist75bpzlimmag_15.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_15,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_16.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_16,overwrite=True)
            t = table.Table(outlist75bpzlimmag_16.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_16,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_17.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_17,overwrite=True)
            t = table.Table(outlist75bpzlimmag_17.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_17,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_18.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_18,overwrite=True)
            t = table.Table(outlist75bpzlimmag_18.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_18,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_19.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_19,overwrite=True)
            t = table.Table(outlist75bpzlimmag_19.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_19,overwrite=True)
            if do50 == True: t = table.Table(outlist50bpzlimmag_20.T, names=names, dtype=dtype); t.write(fout50_bpzlimmag_20,overwrite=True)
            t = table.Table(outlist75bpzlimmag_20.T, names=names, dtype=dtype); t.write(fout75_bpzlimmag_20,overwrite=True)

print("Writing output completed in %0.1f seconds" % ((time.time() - start_write)))
print("Total time for %s: --- %0.2f seconds ---" % ('%s_wghtratios_75_%s_%s_%s_%s_zgap%s' % (fieldID[-26:-20],mskname,lensID,limmag,det,irac,type,zinf,zsup,suffix), (time.time() - start_time)))

print 'Done!'
