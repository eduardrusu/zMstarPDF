# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid
# run as: python fields#.lst masks#.lst samplesize#
# where samplesize number is 100 or 1000

import numpy as np
import scipy
import sys
from scipy import special
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time


# Open the list of fields and associated masks; open the final HE0435 mask, with correct header; open the catalogue, define coordinates and redshift

start_timefield = time.time()

with open(sys.argv[1]) as f:  # fields#.lst
    listfields = f.readlines()

with open(sys.argv[2]) as f:   # masks#.lst
    listmasks = f.readlines()

print str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3])

msk0435 = fits.open('msk0435.fits')
cat0435 = Table.read('HE0435rugiJK_i24nostarXYlepharephotozstellarmasserrbar.cat', names=('X', 'Y', 'ra','dec','z','z_inf','z_sup','mass','mass_inf','mass_sup'), format='ascii')
dist0435 = Column(np.arange(len(cat0435)), name='dist', dtype=('<f8'))  # this is done in order to implement the 1/r weight convention for < 10 arcsec
cat0435.add_column(dist0435)
center0435 = SkyCoord('04:38:14.871 -12:17:14.96', frame='fk5', unit=(u.hourangle, u.deg))
coord0435=SkyCoord(ra=cat0435['ra']*u.degree, dec=cat0435['dec']*u.degree, frame='fk5')
cat0435['dist'] = coord0435.separation(center0435).arcsec
for i in range(len(cat0435)):
    if cat0435['dist'][i] < 10:
        cat0435['dist'][i] = 10

z_s0435 = 1.39
pixlens = 0.200 * u.arcsec


# Open the subfield mask and catalogue, and set up the cell grid

pixCFHT = 0.187 * u.arcsec

for count in range(len(listfields)):
    #if "W1m" in [x[0:len(listmasks[0])] for x in listmasks][count]:
    
        start_timesubfield = time.time()
        msk = fits.open('%s' % [x[0:len(listmasks[0])-1] for x in listmasks][count])
        cells_on_a_side = int((len(msk[0].data[1]) * pixCFHT) / (1200 * pixlens))
        print ("Reading catalogue %s.cat ..." %[x[0:len(listfields[0])-1] for x in listfields][count])
        start_time = time.time()
        catfield = Table.read('%s.cat' % [x[0:len(listfields[0])-1] for x in listfields][count], format='ascii')
        print("--- %s seconds ---" % (time.time() - start_time))

# Define the center of each cell as a matrix of SkyCoord and test which cells to discard because of large masked area (I ignore the area covered in the lens)
# Mask the lens catalogue using the field mask of each cell

        worldfield = WCS('%s' % [x[0:len(listmasks[0])-1] for x in listmasks][count])
        coordorigin = SkyCoord(ra=worldfield.wcs_pix2world(0,0,0)[0]*u.degree, dec=worldfield.wcs_pix2world(0,0,0)[1]*u.degree, frame='fk5')
        centerfields = [ [SkyCoord('00:0:0.0 0:0:0.0', frame='fk5', unit=(u.hourangle, u.deg))]*cells_on_a_side for i in range(cells_on_a_side) ] # matrix of SkyCoord
        maskedcell = np.zeros((cells_on_a_side,cells_on_a_side))
        catupdate0435 = [[ [Table(names=('X','Y','ra','dec','z','z_inf','z_sup','mass','mass_inf','mass_sup','dist'),dtype=('<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))]*cells_on_a_side for i in range(cells_on_a_side) ] for j in range(cells_on_a_side)] # matrix of tables; this is threedimensional because if I create a two dimensional array it will append to all columns
        print "Masking the lens catalogue..."
        start_time = time.time()
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                xlow = (240.0 * u.arcsec / pixCFHT) * i
                xhigh = (240.0 * u.arcsec / pixCFHT) + (240.0 * u.arcsec / pixCFHT) * i
                ylow = (240.0 * u.arcsec / pixCFHT) * j
                yhigh = (240.0 * u.arcsec / pixCFHT) + (240.0 * u.arcsec / pixCFHT) * j
                centerfields[i][j] = SkyCoord(ra=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2, ylow + (yhigh-ylow)/2, 0)[0]*u.degree, dec=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2, ylow + (yhigh-ylow)/2,0)[1]*u.degree, frame='fk5')
                maskedcell[i][j] = msk[0].data[msk[0].data[ylow : yhigh, xlow : xhigh] == 0].size * 1.0 / msk[0].data[ylow : yhigh, xlow : xhigh].size
        #print i,j
                for k in range(len(cat0435)):
                    if (cat0435['z'][k] <= z_s0435) and (msk[0].data[(ylow + cat0435['Y'][k] * pixlens / pixCFHT).value, (xlow + cat0435['X'][k] * pixlens / pixCFHT).value] == 0):
                        catupdate0435[i][j][0].add_row(cat0435[k])


        #print i,j
        #print xlow,xhigh,ylow,yhigh
        #print msk[0].data[xlow : xhigh, ylow : yhigh].size

#print centerfields
        print maskedcell
        msk.close()
        print("--- %s seconds ---" % (time.time() - start_time))
        print "Computing lens weights..."
        start_time = time.time()

# Compute the weights for the lens
# I'm interested in using the probability distribution (PSF) for the z and Mstar of each object; First, I use the most probable values, then I sample from the (symplified) PDFs

        gal_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        zweight_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2rms_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3rms_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        oneoverr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        zoverr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        massoverr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        zmassoverr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        zmass2overr_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overrrms_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overrrms_0435orig = np.zeros((cells_on_a_side,cells_on_a_side))
        # for each cell, I will have 100 or 1000 realizations of the weighted sums
        gal_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        oneoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        zweight_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass2_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass3_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass2rms_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass3rms_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        zoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        massoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass2overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass3overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        zmassoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        zmass2overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass2overrrms_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        mass3overrrms_0435 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                gal_0435orig[i][j]=len(catupdate0435[i][j][0])
                zweight_0435orig[i][j] = np.sum((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z']))
                mass_0435orig[i][j] = np.sum(10**catupdate0435[i][j][0]['mass'])
                mass2_0435orig[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']))
                mass3_0435orig[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']))
                mass2rms_0435orig[i][j] = np.sqrt(mass2_0435orig[i][j])
                mass3rms_0435orig[i][j] = scipy.special.cbrt(mass3_0435orig[i][j])
                oneoverr_0435orig[i][j] = np.sum(1. / catupdate0435[i][j][0]['dist'])
                zoverr_0435orig[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) / catupdate0435[i][j][0]['dist'])
                massoverr_0435orig[i][j] = np.sum(10**catupdate0435[i][j][0]['mass'] / catupdate0435[i][j][0]['dist'])
                mass2overr_0435orig[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                mass3overr_0435orig[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                zmassoverr_0435orig[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                zmass2overr_0435orig[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                mass2overrrms_0435orig[i][j] = np.sqrt(mass2overr_0435orig[i][j])
                mass3overrrms_0435orig[i][j] = scipy.special.cbrt(mass3overr_0435orig[i][j])
                # from here on I sample from the PDFs of z and Mstar. Here I'm using symplified PDFs, since I only have the most probable values and +/-1sigma
                sample0435z=np.zeros((len(catupdate0435[i][j][0]['z']),int(str(sys.argv[3]))))
                # for each object, create 100 or 1000 samples
                sample0435mass=np.zeros((len(catupdate0435[i][j][0]['mass']),int(str(sys.argv[3]))))
                frac0435z=1-1.0*(catupdate0435[i][j][0]['z_sup']-catupdate0435[i][j][0]['z'])/(catupdate0435[i][j][0]['z_sup']-catupdate0435[i][j][0]['z_inf'])
                frac0435mass=1-(1.0*(10**catupdate0435[i][j][0]['mass_sup']-10**catupdate0435[i][j][0]['mass'])/(10**catupdate0435[i][j][0]['mass_sup']-10**catupdate0435[i][j][0]['mass_inf']))
                for k in range(len(catupdate0435[i][j][0]['z'])):
                    sample0435sup=np.zeros(int(str(sys.argv[3])))
                    sample0435inf=np.zeros(int(str(sys.argv[3])))
                    if (frac0435z[k] > 0) and (frac0435z[k] < 1):
                        # in this case I sample from two gaussians, to the left and right of the most probable value
                        sample0435sup=catupdate0435[i][j][0]['z'][k]+abs(np.random.normal(0, catupdate0435[i][j][0]['z_sup'][k]-catupdate0435[i][j][0]['z'][k], int(str(sys.argv[3]))))
                        sample0435inf=catupdate0435[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate0435[i][j][0]['z'][k]-catupdate0435[i][j][0]['z_inf'][k], int(str(sys.argv[3]))))
                        # make sure redshifts are positive
                        while len(sample0435inf[sample0435inf<0]) > 0:
                            sample0435inf[sample0435inf<0]=catupdate0435[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate0435[i][j][0]['z'][k]-catupdate0435[i][j][0]['z_inf'][k], len(sample0435inf[sample0435inf<0])))
                        rand=np.random.random(int(str(sys.argv[3])))
                        sample0435z[k]=sample0435inf
                        sample0435z[k][np.where(rand>frac0435z[k])]=sample0435sup[np.where(rand>frac0435z[k])]
                        #sample0435z[k][sample0435z[k]<0]=0
                    if (frac0435z[k] <= 0) or (frac0435z[k] >= 1):
                        # in this case I ignore the most probable value and I sample from a single gaussian centered in the middle of +/-1sigma
                        sample0435z[k]=np.random.normal(catupdate0435[i][j][0]['z_inf'][k]+(catupdate0435[i][j][0]['z_sup'][k]-catupdate0435[i][j][0]['z_inf'][k])/2, (catupdate0435[i][j][0]['z_sup'][k]-catupdate0435[i][j][0]['z_inf'][k])/2, int(str(sys.argv[3])))
                        #sample0435z[k][sample0435z[k]<0]=0
                    # for mass, assume gaussian distribution in log, not normal space, so I don't get things like negative mass
                    if (frac0435mass[k] > 0) and (frac0435mass[k] < 1):
                        sample0435sup=10**(catupdate0435[i][j][0]['mass'][k]+abs(np.random.normal(0, catupdate0435[i][j][0]['mass_sup'][k]-catupdate0435[i][j][0]['mass'][k], int(str(sys.argv[3])))))
                        sample0435inf=10**(catupdate0435[i][j][0]['mass'][k]-abs(np.random.normal(0, catupdate0435[i][j][0]['mass'][k]-catupdate0435[i][j][0]['mass_inf'][k], int(str(sys.argv[3])))))
                        rand=np.random.random(int(str(sys.argv[3])))
                        sample0435mass[k]=sample0435inf
                        sample0435mass[k][np.where(rand>frac0435mass[k])]=sample0435sup[np.where(rand>frac0435mass[k])]
                    if (frac0435mass[k] <= 0) or (frac0435mass[k] >= 1):
                        sample0435mass[k]=10**np.random.normal(catupdate0435[i][j][0]['mass_inf'][k]+(catupdate0435[i][j][0]['mass_sup'][k]-catupdate0435[i][j][0]['mass_inf'][k])/2, (catupdate0435[i][j][0]['mass_sup'][k]-catupdate0435[i][j][0]['mass_inf'][k])/2, int(str(sys.argv[3])))
                for k in range(int(str(sys.argv[3]))):
                    # sum up the weights of different objects, corresponding to the same realization ID
                    # ignore objects with z>z_source for all weights when drawing from the z PDF
                    gal_0435[i][j][k]=len(sample0435z[:,k][sample0435z[:,k]<z_s0435])
                    oneoverr_0435[i][j][k] = np.sum(1. / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    zweight_0435[i][j][k] = np.sum((z_s0435 * sample0435z[:,k][sample0435z[:,k]<z_s0435]) - (sample0435z[:,k][sample0435z[:,k]<z_s0435]**2))
                    mass_0435[i][j][k] = np.sum(sample0435mass[:,k][sample0435z[:,k]<z_s0435])
                    mass2_0435[i][j][k] = np.sum(sample0435mass[:,k][sample0435z[:,k]<z_s0435]**2)
                    mass3_0435[i][j][k] = np.sum(sample0435mass[:,k][sample0435z[:,k]<z_s0435]**3)
                    mass2rms_0435[i][j][k] = np.sqrt(mass2_0435[i][j][k])
                    mass3rms_0435[i][j][k] = scipy.special.cbrt(mass3_0435[i][j][k])
                    zoverr_0435[i][j][k] = np.sum(((z_s0435 * sample0435z[:,k][sample0435z[:,k]<z_s0435]) - (sample0435z[:,k][sample0435z[:,k]<z_s0435]**2)) / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    massoverr_0435[i][j][k] = np.sum(sample0435mass[:,k][sample0435z[:,k]<z_s0435] / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    mass2overr_0435[i][j][k] = np.sum((sample0435mass[:,k][sample0435z[:,k]<z_s0435]**2) / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    mass3overr_0435[i][j][k] = np.sum((sample0435mass[:,k][sample0435z[:,k]<z_s0435]**3) / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    zmassoverr_0435[i][j][k] = np.sum(((z_s0435 * sample0435z[:,k][sample0435z[:,k]<z_s0435]) - (sample0435z[:,k][sample0435z[:,k]<z_s0435]**2)) * (sample0435mass[:,k][sample0435z[:,k]<z_s0435]) / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    zmass2overr_0435[i][j][k] = np.sum(((z_s0435 * sample0435z[:,k][sample0435z[:,k]<z_s0435]) - (sample0435z[:,k][sample0435z[:,k]<z_s0435]**2)) * (sample0435mass[:,k][sample0435z[:,k]<z_s0435]**2) / catupdate0435[i][j][0]['dist'][sample0435z[:,k]<z_s0435])
                    mass2overrrms_0435[i][j][k] = np.sqrt(mass2overr_0435[i][j][k])
                    mass3overrrms_0435[i][j][k] = scipy.special.cbrt(mass3overr_0435[i][j][k])

        print("--- %s seconds ---" % (time.time() - start_time))

#print gal_0435, sigma_mass2rms_0435

# I think RuntimeWarning: invalid value encountered in double_scalars can be ignored because it shows up only for rms weights, which means it's due to division by zero for the cells completely covered by the masks

# Initialize all field weight measurements as blank matrices

        origfield_gal = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_zweight = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_mass = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_mass2 = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_mass3 = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_zoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_massoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origfield_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side))

# Use objects from the field catalogue that pass the redshift, separation, and mask tests

        print 'Computing field weights...'
        start_time = time.time()
        for i in range(len(catfield)): #range(len(catfield)):   #range(1000):    JUST TO CHECK QUICKLY
            if i in [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000]:
                print i, "objects..."
            if (catfield['Z_B'][i] <= z_s0435) and (catfield['MASK'][i] == 0) and (catfield['star_flag'][i] == 0) and ((catfield['MAG_i'][i] <= 24 and catfield['MAG_i'][i] > 0) or (catfield['MAG_y'][i] <= 24 and catfield['MAG_y'][i] > 0)): # detection in either i or y bands
                coordfieldwcs = SkyCoord(ra=catfield[i]['ALPHA_J2000']*u.degree, dec=catfield[i]['DELTA_J2000']*u.degree, frame='fk5')
        #print coordfieldwcs
                coordfieldpix = worldfield.wcs_world2pix(coordfieldwcs.ra.deg, coordfieldwcs.dec.deg, 0)
                x = int((coordfieldpix[0] * pixCFHT) / (240.0 * u.arcsec))
                y = int((coordfieldpix[1] * pixCFHT) / (240.0 * u.arcsec))
                cellpix_x = coordfieldpix[0] - (x * 240.0 * u.arcsec / pixCFHT)
                cellpix_y = coordfieldpix[1] - (y * 240.0 * u.arcsec / pixCFHT)
                sep = coordfieldwcs.separation(centerfields[x][y]).arcsec
                if (maskedcell[x][y] >= 0.50):
                    if (sep >= 4) and (msk0435[0].data[(cellpix_y * pixCFHT / pixlens).value][(cellpix_x * pixCFHT / pixlens).value] == 0): # ignore objects too close to the center and masked in the lens catalogue
                        origfield_gal[x][y] = origfield_gal[x][y] + 1
                        origfield_zweight[x][y] = origfield_zweight[x][y] + (z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])
                        origfield_mass[x][y] = origfield_mass[x][y] + 10**catfield['LP_log10_SM_MED'][i]
                        origfield_mass2[x][y] = origfield_mass2[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]))
                        origfield_mass3[x][y] = origfield_mass3[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]))
                        # ignore objects with z>z_source for all weights when drawing from the z PDF
                #print (catfield['Z_B_MAX'][i] - catfield['Z_B_MIN'][i])/2
                        if (sep <= 10):
                            origfield_oneoverr[x][y] = origfield_oneoverr[x][y] + 0.1
                            origfield_zoverr[x][y] = origfield_zoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) / 10)
                            origfield_massoverr[x][y] = origfield_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) / 10)
                            origfield_mass2overr[x][y] = origfield_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / 10)
                            origfield_mass3overr[x][y] = origfield_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / 10)
                            origfield_zmassoverr[x][y] = origfield_zmassoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) / 10)
                            origfield_zmass2overr[x][y] = origfield_zmass2overr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) / 10)

                        else:
                            origfield_oneoverr[x][y] = origfield_oneoverr[x][y] + 1. / sep
                            origfield_zoverr[x][y] = origfield_zoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) / sep)
                            origfield_massoverr[x][y] = origfield_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) / sep)
                            origfield_mass2overr[x][y] = origfield_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / sep)
                            origfield_mass3overr[x][y] = origfield_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / sep)
                            origfield_zmassoverr[x][y] = origfield_zmassoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) / sep)
                            origfield_zmass2overr[x][y] = origfield_zmass2overr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) / sep)


        origfield_mass2rms = np.sqrt(origfield_mass2)
        origfield_mass3rms = scipy.special.cbrt(origfield_mass3)
        origfield_mass2overrrms = np.sqrt(origfield_mass2overr)
        origfield_mass3overrrms = scipy.special.cbrt(origfield_mass3overr)

        origq_gal = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_zweight = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass2 = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass3 = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_zoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_massoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass2rms = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass3rms = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass2overrrms = np.zeros((cells_on_a_side,cells_on_a_side))
        origq_mass3overrrms = np.zeros((cells_on_a_side,cells_on_a_side))
        q_gal = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_zweight = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass2 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass3 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_zoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_massoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass2rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass3rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass2overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))
        q_mass3overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[3]))))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50: # skipping the cells with large masked area
                    origq_gal[i][j] = gal_0435orig[i][j] * 1.0 / origfield_gal[i][j]
                    origq_zweight[i][j] = zweight_0435orig[i][j] * 1.0 / origfield_zweight[i][j]
                    origq_mass[i][j] = mass_0435orig[i][j] * 1.0 / origfield_mass[i][j]
                    origq_mass2[i][j] = mass2_0435orig[i][j] * 1.0 / origfield_mass2[i][j]
                    origq_mass3[i][j] = mass3_0435orig[i][j] * 1.0 / origfield_mass3[i][j]
                    origq_oneoverr[i][j] = oneoverr_0435orig[i][j] * 1.0 / origfield_oneoverr[i][j]
                    origq_zoverr[i][j] = zoverr_0435orig[i][j] * 1.0 / origfield_zoverr[i][j]
                    origq_massoverr[i][j] = massoverr_0435orig[i][j] * 1.0 / origfield_massoverr[i][j]
                    origq_mass2overr[i][j] = mass2overr_0435orig[i][j] * 1.0 / origfield_mass2overr[i][j]
                    origq_mass3overr[i][j] = mass3overr_0435orig[i][j] * 1.0 / origfield_mass3overr[i][j]
                    origq_zmassoverr[i][j] = zmassoverr_0435orig[i][j] * 1.0 / origfield_zmassoverr[i][j]
                    origq_zmass2overr[i][j] = zmass2overr_0435orig[i][j] * 1.0 / origfield_zmass2overr[i][j]
                    origq_mass2rms[i][j] = mass2rms_0435orig[i][j] * 1.0 / origfield_mass2rms[i][j]
                    origq_mass3rms[i][j] = mass3rms_0435orig[i][j] * 1.0 / origfield_mass3rms[i][j]
                    origq_mass2overrrms[i][j] = mass2overrrms_0435orig[i][j] * 1.0 / origfield_mass2overrrms[i][j]
                    origq_mass3overrrms[i][j] = mass3overrrms_0435orig[i][j] * 1.0 / origfield_mass3overrrms[i][j]
                    q_gal[i][j] = gal_0435[i][j] * 1.0 / origfield_gal[i][j]
                    q_oneoverr[i][j] = oneoverr_0435[i][j] * 1.0 / origfield_oneoverr[i][j]
                    q_zweight[i][j] = zweight_0435[i][j] * 1.0 / origfield_zweight[i][j]
                    q_mass[i][j] = mass_0435[i][j] * 1.0 / origfield_mass[i][j]
                    q_mass2[i][j] = mass2_0435[i][j] * 1.0 / origfield_mass2[i][j]
                    q_mass3[i][j] = mass3_0435[i][j] * 1.0 / origfield_mass3[i][j]
                    q_zoverr[i][j] = zoverr_0435[i][j] * 1.0 / origfield_zoverr[i][j]
                    q_massoverr[i][j] = massoverr_0435[i][j] * 1.0 / origfield_massoverr[i][j]
                    q_mass2overr[i][j] = mass2overr_0435[i][j] * 1.0 / origfield_mass2overr[i][j]
                    q_mass3overr[i][j] = mass3overr_0435[i][j] * 1.0 / origfield_mass3overr[i][j]
                    q_zmassoverr[i][j] = zmassoverr_0435[i][j] * 1.0 / origfield_zmassoverr[i][j]
                    q_zmass2overr[i][j] = zmass2overr_0435[i][j] * 1.0 / origfield_zmass2overr[i][j]
                    q_mass2rms[i][j] = mass2rms_0435[i][j] * 1.0 / origfield_mass2rms[i][j]
                    q_mass3rms[i][j] = mass3rms_0435[i][j] * 1.0 / origfield_mass3rms[i][j]
                    q_mass2overrrms[i][j] = mass2overrrms_0435[i][j] * 1.0 / origfield_mass2overrrms[i][j]
                    q_mass3overrrms[i][j] = mass3overrrms_0435[i][j] * 1.0 / origfield_mass3overrrms[i][j]


        print("--- %s seconds ---" % (time.time() - start_time))

        print ("No. of fields with area without mask > 75 percent & > 50 percent: %d/%d, %d/%d" % ((maskedcell>=0.75).sum(), cells_on_a_side ** 2, (maskedcell>=0.50).sum(), cells_on_a_side ** 2))

# write a file where the most probable value of z and Mstar has been used, and a file for each of the 100 or 1000 samples
        f = open('%s_q50_orig_noCFHTLENSsamp.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')
        g = open('%s_q50_samp_noCFHTLENSsamp.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (origq_gal[i][j], origq_oneoverr[i][j], origq_zweight[i][j], origq_mass[i][j], origq_mass2[i][j], origq_mass2rms[i][j], origq_mass3[i][j], origq_mass3rms[i][j], origq_zoverr[i][j], origq_massoverr[i][j], origq_mass2overr[i][j], origq_mass3overr[i][j], origq_mass2overrrms[i][j], origq_mass3overrrms[i][j], origq_zmassoverr[i][j], origq_zmass2overr[i][j]))
                    for k in range(int(str(sys.argv[3]))):
                        g.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j][k], q_oneoverr[i][j][k], q_zweight[i][j][k], q_mass[i][j][k], q_mass2[i][j][k], q_mass2rms[i][j][k], q_mass3[i][j][k], q_mass3rms[i][j][k], q_zoverr[i][j][k], q_massoverr[i][j][k], q_mass2overr[i][j][k], q_mass3overr[i][j][k], q_mass2overrrms[i][j][k], q_mass3overrrms[i][j][k], q_zmassoverr[i][j][k], q_zmass2overr[i][j][k]))

        f.close()
        g.close()

        f = open('%s_q75_orig_noCFHTLENSsamp.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')
        g = open('%s_q75_samp_noCFHTLENSsamp.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.75:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (origq_gal[i][j], origq_oneoverr[i][j], origq_zweight[i][j], origq_mass[i][j], origq_mass2[i][j], origq_mass2rms[i][j], origq_mass3[i][j], origq_mass3rms[i][j], origq_zoverr[i][j], origq_massoverr[i][j], origq_mass2overr[i][j], origq_mass3overr[i][j], origq_mass2overrrms[i][j], origq_mass3overrrms[i][j], origq_zmassoverr[i][j], origq_zmass2overr[i][j]))
                    for k in range(int(str(sys.argv[3]))):
                        g.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j][k], q_oneoverr[i][j][k], q_zweight[i][j][k], q_mass[i][j][k], q_mass2[i][j][k], q_mass2rms[i][j][k], q_mass3[i][j][k], q_mass3rms[i][j][k], q_zoverr[i][j][k], q_massoverr[i][j][k], q_mass2overr[i][j][k], q_mass3overr[i][j][k], q_mass2overrrms[i][j][k], q_mass3overrrms[i][j][k], q_zmassoverr[i][j][k], q_zmass2overr[i][j][k]))

        f.close()
        g.close()

        print("Total time for subfield: --- %s seconds ---" % (time.time() - start_timesubfield))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

