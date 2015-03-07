# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid

import numpy as np
import scipy
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

with open('fields1.lst') as f:
    listfields = f.readlines()

with open('masks1.lst') as f:
    listmasks = f.readlines()

msk0435 = fits.open('msk0435.fits')
cat0435 = Table.read('HE0435rugiJK_i24nostarXYlepharephotozstellarmass.cat', names=('X', 'Y', 'ra','dec','z','mass'), format='ascii')
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
        catupdate0435 = [[ [Table(names=('X','Y','ra','dec','z','mass','dist'),dtype=('<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))]*cells_on_a_side for i in range(cells_on_a_side) ] for j in range(cells_on_a_side)] # matrix of tables; this is threedimensional because if I create a two dimensional array it will append to all columns
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

# Compute the weights for the lens

        gal_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        zweight_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2rms_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3rms_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        oneoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        zoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        massoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        zmassoverr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        zmass2overr_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overrrms_0435 = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overrrms_0435 = np.zeros((cells_on_a_side,cells_on_a_side))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                gal_0435[i][j]=len(catupdate0435[i][j][0])
                zweight_0435[i][j] = np.sum((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z']))
                mass_0435[i][j] = np.sum(10**catupdate0435[i][j][0]['mass'])
                mass2_0435[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']))
                mass3_0435[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']))
                mass2rms_0435[i][j] = np.sqrt(mass2_0435[i][j])
                mass3rms_0435[i][j] = scipy.special.cbrt(mass3_0435[i][j])
                oneoverr_0435[i][j] = np.sum(1. / catupdate0435[i][j][0]['dist'])
                zoverr_0435[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) / catupdate0435[i][j][0]['dist'])
                massoverr_0435[i][j] = np.sum(10**catupdate0435[i][j][0]['mass'] / catupdate0435[i][j][0]['dist'])
                mass2overr_0435[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                mass3overr_0435[i][j] = np.sum((10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                zmassoverr_0435[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                zmass2overr_0435[i][j] = np.sum(((z_s0435 * catupdate0435[i][j][0]['z']) - (catupdate0435[i][j][0]['z'] * catupdate0435[i][j][0]['z'])) * (10**catupdate0435[i][j][0]['mass']) * (10**catupdate0435[i][j][0]['mass']) / catupdate0435[i][j][0]['dist'])
                mass2overrrms_0435[i][j] = np.sqrt(mass2overr_0435[i][j])
                mass3overrrms_0435[i][j] = scipy.special.cbrt(mass3overr_0435[i][j])

# Initialize all field weight measurements as blank matrices

        field_gal = np.zeros((cells_on_a_side,cells_on_a_side))
        field_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_zweight = np.zeros((cells_on_a_side,cells_on_a_side))
        field_mass = np.zeros((cells_on_a_side,cells_on_a_side))
        field_mass2 = np.zeros((cells_on_a_side,cells_on_a_side))
        field_mass3 = np.zeros((cells_on_a_side,cells_on_a_side))
        field_zoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_massoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        field_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side))

# Use objects from the field catalogue that pass the redshift, separation, and mask tests

        print 'Computing weights...'
        start_time = time.time()
        for i in range(len(catfield)): #range(len(catfield)):   #range(1000):    JUST TO CHECK QUICKLY
            if i in [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000]:
                print i, "objects..."
            if (catfield['Z_B'][i] <= z_s0435) and (catfield['MASK'][i] == 0) and (catfield['star_flag'][i] == 0) and ((catfield['MAG_i'][i] <= 24 and catfield['MAG_i'][i] > 0) or (catfield['MAG_y'][i] <= 24 and catfield['MAG_y'][i] > 0)):
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
                        field_gal[x][y] = field_gal[x][y] + 1
                        field_zweight[x][y] = field_zweight[x][y] + (z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])
                        field_mass[x][y] = field_mass[x][y] + 10**catfield['LP_log10_SM_MED'][i]
                        field_mass2[x][y] = field_mass2[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]))
                        field_mass3[x][y] = field_mass3[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]))
                        if (sep <= 10):
                            field_oneoverr[x][y] = field_oneoverr[x][y] + 0.1
                            field_zoverr[x][y] = field_zoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) / 10)
                            field_massoverr[x][y] = field_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) / 10)
                            field_mass2overr[x][y] = field_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / 10)
                            field_mass3overr[x][y] = field_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / 10)
                            field_zmassoverr[x][y] = field_zmassoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) / 10)
                            field_zmass2overr[x][y] = field_zmass2overr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) / 10)
                        else:
                            field_oneoverr[x][y] = field_oneoverr[x][y] + 1. / sep
                            field_zoverr[x][y] = field_zoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) / sep)
                            field_massoverr[x][y] = field_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][i]) / sep)
                            field_mass2overr[x][y] = field_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / sep)
                            field_mass3overr[x][y] = field_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i])) / sep)
                            field_zmassoverr[x][y] = field_zmassoverr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) / sep)
                            field_zmass2overr[x][y] = field_zmass2overr[x][y] + (((z_s0435 * catfield['Z_B'][i]) - (catfield['Z_B'][i] * catfield['Z_B'][i])) * (10**catfield['LP_log10_SM_MED'][i]) * (10**catfield['LP_log10_SM_MED'][i]) / sep)


        field_mass2rms = np.sqrt(field_mass2)
        field_mass3rms = scipy.special.cbrt(field_mass3)
        field_mass2overrrms = np.sqrt(field_mass2overr)
        field_mass3overrrms = scipy.special.cbrt(field_mass3overr)

        q_gal = np.zeros((cells_on_a_side,cells_on_a_side))
        q_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_zweight = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass2 = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass3 = np.zeros((cells_on_a_side,cells_on_a_side))
        q_zoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_massoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass2rms = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass3rms = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass2overrrms = np.zeros((cells_on_a_side,cells_on_a_side))
        q_mass3overrrms = np.zeros((cells_on_a_side,cells_on_a_side))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50: # skipping the cells with large masked area
                    q_gal[i][j] = gal_0435[i][j] * 1.0 / field_gal[i][j]
                    q_zweight[i][j] = zweight_0435[i][j] * 1.0 / field_zweight[i][j]
                    q_mass[i][j] = mass_0435[i][j] * 1.0 / field_mass[i][j]
                    q_mass2[i][j] = mass2_0435[i][j] * 1.0 / field_mass2[i][j]
                    q_mass3[i][j] = mass3_0435[i][j] * 1.0 / field_mass3[i][j]
                    q_oneoverr[i][j] = oneoverr_0435[i][j] * 1.0 / field_oneoverr[i][j]
                    q_zoverr[i][j] = zoverr_0435[i][j] * 1.0 / field_zoverr[i][j]
                    q_massoverr[i][j] = massoverr_0435[i][j] * 1.0 / field_massoverr[i][j]
                    q_mass2overr[i][j] = mass2overr_0435[i][j] * 1.0 / field_mass2overr[i][j]
                    q_mass3overr[i][j] = mass3overr_0435[i][j] * 1.0 / field_mass3overr[i][j]
                    q_zmassoverr[i][j] = zmassoverr_0435[i][j] * 1.0 / field_zmassoverr[i][j]
                    q_zmass2overr[i][j] = zmass2overr_0435[i][j] * 1.0 / field_zmass2overr[i][j]
                    q_mass2rms[i][j] = mass2rms_0435[i][j] * 1.0 / field_mass2rms[i][j]
                    q_mass3rms[i][j] = mass3rms_0435[i][j] * 1.0 / field_mass3rms[i][j]
                    q_mass2overrrms[i][j] = mass2overrrms_0435[i][j] * 1.0 / field_mass2overrrms[i][j]
                    q_mass3overrrms[i][j] = mass3overrrms_0435[i][j] * 1.0 / field_mass3overrrms[i][j]

#print field_gal
#print field_zweight
#print field_mass
#print field_mass2
#print field_mass3
#print field_oneoverr
#print field_zoverr
#print field_massoverr
#print field_mass2overr
#print field_mass3overr
#print field_zmassoverr
#print field_zmass2overr
#print field_mass2rms
#print field_mass3rms
#print field_mass2overrrms
#print field_mass3overrrms

        print("--- %s seconds ---" % (time.time() - start_time))

        print ("No. of fields with area without mask > 75 percent & > 50 percent: %d/%d, %d/%d" % ((maskedcell>=0.75).sum(), cells_on_a_side ** 2, (maskedcell>=0.50).sum(), cells_on_a_side ** 2))


        f = open('%s_q50.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j], q_oneoverr[i][j], q_zweight[i][j], q_mass[i][j], q_mass2[i][j], q_mass2rms[i][j], q_mass3[i][j], q_mass3rms[i][j], q_zoverr[i][j], q_massoverr[i][j], q_mass2overr[i][j], q_mass3overr[i][j], q_mass2overrrms[i][j], q_mass3overrrms[i][j], q_zmassoverr[i][j], q_zmass2overr[i][j]))

        f.close()

        f = open('%s_q75.lst' % [x[0:len(listfields[0])-1] for x in listfields][count],'w')

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.75:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j], q_oneoverr[i][j], q_zweight[i][j], q_mass[i][j], q_mass2[i][j], q_mass2rms[i][j], q_mass3[i][j], q_mass3rms[i][j], q_zoverr[i][j], q_massoverr[i][j], q_mass2overr[i][j], q_mass3overr[i][j], q_mass2overrrms[i][j], q_mass3overrrms[i][j], q_zmassoverr[i][j], q_zmass2overr[i][j]))

        f.close()

        print("Total time for subfield: --- %s seconds ---" % (time.time() - start_timesubfield))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

