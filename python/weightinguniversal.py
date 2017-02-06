# fields CFHTLenS W1-4
# subfields: 171 1deg^2 throughout W1-4
# cells: 4x4arcmin covering each subfield, in a grid
# run as: python weightinguniversal.py lens fields#.lst masks#.lst samplesize# msk_lenssize maglimit classification mode
# where lens is B1608,HE0435,HE1104 or RX1131, samplesize number is 0, 100 or 1000 and msk_lenssize is 45, 60, 90 or 120; maglimit is 23 23.5 or 24; classification is "old" or "new", where "old" means original CFHTLENS, and "new" means extended beyond i=23;

import numpy as np
import scipy
import sys
from scipy import special
from scipy import stats
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column
import time


# Open the list of fields and associated masks; open the final HE_lens mask, with correct header; open the catalogue, define coordinates and redshift

start_timefield = time.time()

with open(sys.argv[2]) as f:  # fields#.lst
    listfields = f.readlines()

with open(sys.argv[3]) as f:   # masks#.lst
    listmasks = f.readlines()

print("Arguments: \n Lens field: %s \n List of fields to work on: %s \n List of masks associated with the fields: %s \n Number of samples to be drawn from P(z) and P(Mstar): %s \n Radius of each cell: %s \n Limiting magnitude: %s \n Classification: %s" % (str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]), float(str(sys.argv[6])), str(sys.argv[7])))

msk_lens = fits.open('msk%s_asecrad%s.fits' % (str(sys.argv[1]),str(sys.argv[5])))
maskedlens = np.where( msk_lens[0].data == 0 ) # will be used later in computing the mask cover fraction
#print maskedlens[0]
#print maskedlens[1]
cat_lens = Table.read('%s.cat' % (str(sys.argv[1])), names=('ID', 'X', 'Y', 'ra','dec','mag','z','z_inf','z_sup','mass_best','mass_inf','mass_med','mass_sup', 'specz?', 'correction'), format='ascii')
pdzcat_lens = np.loadtxt('%s_pdz.cat' % str(sys.argv[1])) # since the catalogue is small, I just load it into memory
mstarcat_lens = np.loadtxt('%s_mstar.cat' % str(sys.argv[1])) # since the catalogue is small, I just load it into memory
zgrid=np.linspace(0.05,3.5,70)
zgridint = np.arange(70) # because stats.rv_discrete only works with integer points
if str(sys.argv[7])=="old":
    cat_lens_reclassif=cat_lens
else:
    cat_lens_reclassif=Table(names=('ID','X','Y','ra','dec','mag','z','z_inf','z_sup','mass_best','mass_inf','mass_med','mass_sup','specz?','correction'),dtype=('<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))
    for i in range(len(cat_lens)):
        if cat_lens['mag'][i]<=23:
            cat_lens_reclassif.add_row(cat_lens[i])
        else:
            if cat_lens['specz?'][i]==1:
                cat_lens_reclassif.add_row(cat_lens[i])
            elif mstarcat_lens[mstarcat_lens[:,0]==cat_lens['ID'][i]][0][5]/mstarcat_lens[mstarcat_lens[:,0]==cat_lens['ID'][i]][0][6]<0.5:
                cat_lens_reclassif.add_row(cat_lens[i])
#print cat_lens_reclassif['ID']
dist_lens = Column(np.arange(len(cat_lens_reclassif)), name='dist', dtype=('<f8'))  # this is done in order to implement the 1/r weight convention for < 10 arcsec and also in order to use the mask of appropriate radius
cat_lens_reclassif.add_column(dist_lens)
if str(sys.argv[1]) == "B1608":
    center_lens = SkyCoord('16:09:13.956 +65:32:28.00', frame='fk5', unit=(u.hourangle, u.deg))
if str(sys.argv[1]) == "HE0435":
    center_lens = SkyCoord('04:38:14.871 -12:17:14.96', frame='fk5', unit=(u.hourangle, u.deg))
if str(sys.argv[1]) == "HE1104":
    center_lens = SkyCoord('11:06:33.450 -18:21:24.20', frame='fk5', unit=(u.hourangle, u.deg))
if str(sys.argv[1]) == "RX1131":
    center_lens = SkyCoord('11:31:51.435 -12:31:58.24', frame='fk5', unit=(u.hourangle, u.deg))

coord_lens=SkyCoord(ra=cat_lens_reclassif['ra']*u.degree, dec=cat_lens_reclassif['dec']*u.degree, frame='fk5')
cat_lens_reclassif['dist'] = coord_lens.separation(center_lens).arcsec
for i in range(len(cat_lens_reclassif)):
    if cat_lens_reclassif['dist'][i] < 10:
        cat_lens_reclassif['dist'][i] = 10

if str(sys.argv[1]) == "B1608":
    z_s_lens = 1.39
if str(sys.argv[1]) == "HE0435":
    z_s_lens = 1.69
if str(sys.argv[1]) == "HE1104":
    z_s_lens = 2.32
if str(sys.argv[1]) == "RX1131":
    z_s_lens = 0.66
pixlens = 0.200 * u.arcsec

# Open the subfield mask and catalogue, and set up the cell grid

pixCFHT = 0.187 * u.arcsec

for count in range(len(listfields)):
    #if "W1m" in [x[0:len(listmasks[0])] for x in listmasks][count]:
    
        start_timesubfield = time.time()
        msk = fits.open('%s' % [a[0:len(listmasks[0])-1] for a in listmasks][count])
        cells_on_a_side = int((len(msk[0].data[1]) * pixCFHT) / (1200 * pixlens))
        #print ("Reading catalogue %s.cat ..." %[x[0:len(listfields[0])-1] for x in listfields][count])
        #start_time = time.time()
        #catfield = Table.read('%s.cat' % [x[0:len(listfields[0])-1] for x in listfields][count], format='ascii')
        #print("--- %s seconds ---" % (time.time() - start_time))

# Define the center of each cell as a matrix of SkyCoord and test which cells to discard because of large masked area (I ignore the area covered in the lens)
# Mask the lens catalogue using the field mask of each cell

        worldfield = WCS('%s' % [b[0:len(listmasks[0])-1] for b in listmasks][count])
        coordorigin = SkyCoord(ra=worldfield.wcs_pix2world(0,0,0)[0]*u.degree, dec=worldfield.wcs_pix2world(0,0,0)[1]*u.degree, frame='fk5')
        centerfields = [ [SkyCoord('00:0:0.0 0:0:0.0', frame='fk5', unit=(u.hourangle, u.deg))]*cells_on_a_side for i in range(cells_on_a_side) ] # matrix of SkyCoord
        maskedcell = np.zeros((cells_on_a_side,cells_on_a_side))
        origtry = np.zeros((cells_on_a_side,cells_on_a_side))
        origtry1 = np.zeros((cells_on_a_side,cells_on_a_side))
        catupdate_lens = [[ [Table(names=('ID','X','Y','ra','dec','mag','z','z_inf','z_sup','mass_best','mass_inf','mass_med','mass_sup','specz?','correction','dist'),dtype=('<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))]*cells_on_a_side for i in range(cells_on_a_side) ] for j in range(cells_on_a_side)] # matrix of tables; this is threedimensional because if I create a two dimensional array it will append to all columns
        print "Masking the lens catalogue..."
        start_time = time.time()
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                xlow = (240.0 * u.arcsec / pixCFHT) * i
                xhigh = (240.0 * u.arcsec / pixCFHT) + (240.0 * u.arcsec / pixCFHT) * i
                ylow = (240.0 * u.arcsec / pixCFHT) * j
                yhigh = (240.0 * u.arcsec / pixCFHT) + (240.0 * u.arcsec / pixCFHT) * j
                centerfields[i][j] = SkyCoord(ra=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2, ylow + (yhigh-ylow)/2, 0)[0]*u.degree, dec=worldfield.wcs_pix2world(xlow + (xhigh-xlow)/2, ylow + (yhigh-ylow)/2,0)[1]*u.degree, frame='fk5')
                maskedfieldx = float(xlow)+(maskedlens[0]/1200.0)*(float(xhigh)-float(xlow))
                maskedfieldy = float(ylow)+(maskedlens[1]/1200.0)*(float(yhigh)-float(ylow))
                # the pixels in the field mask that correspond to the unmasked pixels in the lens mask
                maskedfield = msk[0].data[maskedfieldy.astype(int),maskedfieldx.astype(int)]
                maskedcell[i][j] = maskedfield[maskedfield == 0].size / (np.pi*((int((str(sys.argv[5])))*u.arcsec/pixlens)**2))
        #print i,j
                for k in range(len(cat_lens_reclassif)):
                    if (msk_lens[0].data[cat_lens_reclassif['Y'][k],cat_lens_reclassif['X'][k]] == 0) and (cat_lens_reclassif['z'][k] <= z_s_lens) and (cat_lens_reclassif['mag'][k] <= float(str(sys.argv[6]))) and (msk[0].data[(ylow + cat_lens_reclassif['Y'][k] * pixlens / pixCFHT).value, (xlow + cat_lens_reclassif['X'][k] * pixlens / pixCFHT).value] == 0):
                        #print i,j,k,len(cat_lens_reclassif),cat_lens[k]
                        catupdate_lens[i][j][0].add_row(cat_lens_reclassif[k])
                        # the first condition is necessary because, the way the code is written now, the lens field mask is not used to actually mask the catalogue

        #print catupdate_lens[5][10][0]['dist']
        #print i,j
        #print xlow,xhigh,ylow,yhigh
        #print msk[0].data[xlow : xhigh, ylow : yhigh].size

#print centerfields
        #print maskedcell
        msk.close()
        print("--- %s seconds ---" % (time.time() - start_time))
        print "Computing lens weights..."
        start_time = time.time()

# Compute the weights for the lens
# I'm interested in using the probability distribution (PSF) for the z and Mstar of each object; First, I use the most probable values, then I sample from the (symplified) PDFs

        #if str(sys.argv[8]) == "orig":
        gal_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        zweight_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2rms_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3rms_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        oneoverr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        zoverr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        massoverr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        zmassoverr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        zmass2overr_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass2overrrms_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        mass3overrrms_lensorig = np.zeros((cells_on_a_side,cells_on_a_side))
        #if str(sys.argv[8]) == "samp":
            # for each cell, I will have 100 or 1000 realizations of the weighted sums
        gal_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        oneoverr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zweight_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2rms_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3rms_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zoverr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        massoverr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2overr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3overr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zmassoverr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zmass2overr_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2overrrms_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3overrrms_lens = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        #if str(sys.argv[8]) == "tab":
            # tabulated PDF for z and Mstar
        gal_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        oneoverr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zweight_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2rms_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3rms_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zoverr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        massoverr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2overr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3overr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zmassoverr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        zmass2overr_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass2overrrms_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        mass3overrrms_lens_tab = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                #if str(sys.argv[8]) == "orig":
                    gal_lensorig[i][j]=len(catupdate_lens[i][j][0])
                    zweight_lensorig[i][j] = np.sum((z_s_lens * catupdate_lens[i][j][0]['z']) - (catupdate_lens[i][j][0]['z'] * catupdate_lens[i][j][0]['z']))
                    mass_lensorig[i][j] = np.sum(10**catupdate_lens[i][j][0]['mass_med'])
                    mass2_lensorig[i][j] = np.sum((10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']))
                    mass3_lensorig[i][j] = np.sum((10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']))
                    mass2rms_lensorig[i][j] = np.sqrt(mass2_lensorig[i][j])
                    mass3rms_lensorig[i][j] = scipy.special.cbrt(mass3_lensorig[i][j])
                    oneoverr_lensorig[i][j] = np.sum(1. / catupdate_lens[i][j][0]['dist'])
                    zoverr_lensorig[i][j] = np.sum(((z_s_lens * catupdate_lens[i][j][0]['z']) - (catupdate_lens[i][j][0]['z'] * catupdate_lens[i][j][0]['z'])) / catupdate_lens[i][j][0]['dist'])
                    massoverr_lensorig[i][j] = np.sum(10**catupdate_lens[i][j][0]['mass_med'] / catupdate_lens[i][j][0]['dist'])
                    mass2overr_lensorig[i][j] = np.sum((10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']) / catupdate_lens[i][j][0]['dist'])
                    mass3overr_lensorig[i][j] = np.sum((10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']) / catupdate_lens[i][j][0]['dist'])
                    zmassoverr_lensorig[i][j] = np.sum(((z_s_lens * catupdate_lens[i][j][0]['z']) - (catupdate_lens[i][j][0]['z'] * catupdate_lens[i][j][0]['z'])) * (10**catupdate_lens[i][j][0]['mass_med']) / catupdate_lens[i][j][0]['dist'])
                    zmass2overr_lensorig[i][j] = np.sum(((z_s_lens * catupdate_lens[i][j][0]['z']) - (catupdate_lens[i][j][0]['z'] * catupdate_lens[i][j][0]['z'])) * (10**catupdate_lens[i][j][0]['mass_med']) * (10**catupdate_lens[i][j][0]['mass_med']) / catupdate_lens[i][j][0]['dist'])
                    mass2overrrms_lensorig[i][j] = np.sqrt(mass2overr_lensorig[i][j])
                    mass3overrrms_lensorig[i][j] = scipy.special.cbrt(mass3overr_lensorig[i][j])

                # from here on I sample from the PDFs of z and Mstar. Here I'm using symplified PDFs, since I only have the most probable values and +/-1sigma
                #if str(sys.argv[8]) == "samp":
                    sample_lensz=np.zeros((len(catupdate_lens[i][j][0]['z']),int(str(sys.argv[4]))))
                    # for each object, create 100 or 1000 samples
                    sample_lensmass=np.zeros((len(catupdate_lens[i][j][0]['mass_med']),int(str(sys.argv[4]))))
                    frac_lensz=1-1.0*(catupdate_lens[i][j][0]['z_sup']-catupdate_lens[i][j][0]['z'])/(catupdate_lens[i][j][0]['z_sup']-catupdate_lens[i][j][0]['z_inf'])
                    frac_lensmass=1-(1.0*(10**catupdate_lens[i][j][0]['mass_sup']-10**catupdate_lens[i][j][0]['mass_med'])/(10**catupdate_lens[i][j][0]['mass_sup']-10**catupdate_lens[i][j][0]['mass_inf']))
                    for k in range(len(catupdate_lens[i][j][0]['z'])): # for each object, create its samples
                        sample_lenssup=np.zeros(int(str(sys.argv[4])))
                        sample_lensinf=np.zeros(int(str(sys.argv[4])))
                        if (frac_lensz[k] > 0) and (frac_lensz[k] < 1):
                            # in this case I sample from two gaussians, to the left and right of the most probable value
                            sample_lenssup=catupdate_lens[i][j][0]['z'][k]+abs(np.random.normal(0, catupdate_lens[i][j][0]['z_sup'][k]-catupdate_lens[i][j][0]['z'][k], int(str(sys.argv[4]))))
                            sample_lensinf=catupdate_lens[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['z'][k]-catupdate_lens[i][j][0]['z_inf'][k], int(str(sys.argv[4]))))
                            # make sure redshifts are positive
                            while len(sample_lensinf[sample_lensinf<0]) > 0:
                                sample_lensinf[sample_lensinf<0]=catupdate_lens[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['z'][k]-catupdate_lens[i][j][0]['z_inf'][k], len(sample_lensinf[sample_lensinf<0])))
                            rand=np.random.random(int(str(sys.argv[4])))
                            sample_lensz[k]=sample_lensinf
                            sample_lensz[k][np.where(rand>frac_lensz[k])]=sample_lenssup[np.where(rand>frac_lensz[k])]
                            #sample_lensz[k][sample_lensz[k]<0]=0
                        if (frac_lensz[k] <= 0) or (frac_lensz[k] >= 1):
                            # in this case I ignore the most probable value and I sample from a single gaussian centered in the middle of +/-1sigma
                            sample_lensz[k]=np.random.normal(catupdate_lens[i][j][0]['z_inf'][k]+(catupdate_lens[i][j][0]['z_sup'][k]-catupdate_lens[i][j][0]['z_inf'][k])/2, (catupdate_lens[i][j][0]['z_sup'][k]-catupdate_lens[i][j][0]['z_inf'][k])/2, int(str(sys.argv[4])))
                            #sample_lensz[k][sample_lensz[k]<0]=0
                        # for mass, assume gaussian distribution in log, not normal space, so I don't get things like negative mass
                        if (frac_lensmass[k] > 0) and (frac_lensmass[k] < 1):
                            sample_lenssup=10**(catupdate_lens[i][j][0]['mass_med'][k]+abs(np.random.normal(0, catupdate_lens[i][j][0]['mass_sup'][k]-catupdate_lens[i][j][0]['mass_med'][k], int(str(sys.argv[4])))))
                            sample_lensinf=10**(catupdate_lens[i][j][0]['mass_med'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['mass_med'][k]-catupdate_lens[i][j][0]['mass_inf'][k], int(str(sys.argv[4])))))
                            rand=np.random.random(int(str(sys.argv[4])))
                            sample_lensmass[k]=sample_lensinf
                            sample_lensmass[k][np.where(rand>frac_lensmass[k])]=sample_lenssup[np.where(rand>frac_lensmass[k])] # sampling to the left or right
                        if (frac_lensmass[k] <= 0) or (frac_lensmass[k] >= 1):
                            sample_lensmass[k]=10**np.random.normal(catupdate_lens[i][j][0]['mass_inf'][k]+(catupdate_lens[i][j][0]['mass_sup'][k]-catupdate_lens[i][j][0]['mass_inf'][k])/2, (catupdate_lens[i][j][0]['mass_sup'][k]-catupdate_lens[i][j][0]['mass_inf'][k])/2, int(str(sys.argv[4])))
                    for k in range(int(str(sys.argv[4]))):
                        # for each sample ID, sum up the weights of different objects
                        # ignore objects with z>z_source for all weights when drawing from the z PDF
                        gal_lens[i][j][k]=len(sample_lensz[:,k][sample_lensz[:,k]<z_s_lens])
                        oneoverr_lens[i][j][k] = np.sum(1. / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zweight_lens[i][j][k] = np.sum((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2))
                        mass_lens[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens])
                        mass2_lens[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2)
                        mass3_lens[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**3)
                        mass2rms_lens[i][j][k] = np.sqrt(mass2_lens[i][j][k])
                        mass3rms_lens[i][j][k] = scipy.special.cbrt(mass3_lens[i][j][k])
                        zoverr_lens[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        massoverr_lens[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens] / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass2overr_lens[i][j][k] = np.sum((sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass3overr_lens[i][j][k] = np.sum((sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**3) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zmassoverr_lens[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) * (sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zmass2overr_lens[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) * (sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass2overrrms_lens[i][j][k] = np.sqrt(mass2overr_lens[i][j][k])
                        mass3overrrms_lens[i][j][k] = scipy.special.cbrt(mass3overr_lens[i][j][k])

                #if str(sys.argv[8]) == "tab":
                    sample_lensz=np.zeros((len(catupdate_lens[i][j][0]['z']),int(str(sys.argv[4]))))
                    # for each object, create 100 or 1000 samples
                    sample_lensmass=np.zeros((len(catupdate_lens[i][j][0]['mass_med']),int(str(sys.argv[4]))))
                    frac_lensz=1-1.0*(catupdate_lens[i][j][0]['z_sup']-catupdate_lens[i][j][0]['z'])/(catupdate_lens[i][j][0]['z_sup']-catupdate_lens[i][j][0]['z_inf'])
                    frac_lensmass=1-(1.0*(10**catupdate_lens[i][j][0]['mass_sup']-10**catupdate_lens[i][j][0]['mass_med'])/(10**catupdate_lens[i][j][0]['mass_sup']-10**catupdate_lens[i][j][0]['mass_inf']))
                    for k in range(len(catupdate_lens[i][j][0]['z'])): # for each object, create its samples
                        #print "cat",k
                        if catupdate_lens[i][j][0]['specz?'][k]==0: # tabulated values for non-specz objects
                            dummy=pdzcat_lens[pdzcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][1:]
                            dummy[dummy<0.001]=0
                            z_tab=dummy
                            #print "z_tab", catupdate_lens[i][j][0]['ID'][k], z_tab
                            massbest_tab=np.zeros(70)
                            massinf_tab=np.zeros(70)
                            massmed_tab=np.zeros(70)
                            masssup_tab=np.zeros(70)
                            for l in range(70): # for each of the tabulated points
                                if z_tab[l]!=0:
                                    massbest_tab[l]=mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][7+l*4]+catupdate_lens[i][j][0]['correction'][k]
                                    if mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][9+l*4]>=0:
                                        massmed_tab[l]=mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][9+l*4]+catupdate_lens[i][j][0]['correction'][k]
                                    else:
                                        massmed_tab[l]=massbest_tab[l]
                                    if mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][8+l*4]>=0:
                                        massinf_tab[l]=mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][8+l*4]+catupdate_lens[i][j][0]['correction'][k]
                                    else:
                                        massinf_tab[l]=massmed_tab[l]-0.1
                                    if mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][10+l*4]>=0:
                                        masssup_tab[l]=mstarcat_lens[mstarcat_lens[:,0]==catupdate_lens[i][j][0]['ID'][k]][0][10+l*4]+catupdate_lens[i][j][0]['correction'][k]
                                    else:
                                        masssup_tab[l]=massmed_tab[l]+0.1
                            #print "massbest_tab", massbest_tab
                            #print "massinf_tab", massinf_tab
                            #print "massmed_tab", massmed_tab
                            #print "masssup_tab", masssup_tab
                            custm = stats.rv_discrete(name='custm', values=(zgridint, z_tab)) # sample from the tabulated distribution
                            sample=custm.rvs(size=int(str(sys.argv[4]))) # Weird, sometimes custm samples from places where the distribution is 0. Below I ensure this doesn't happen. The cause is probably due to the probabilities not summing exactly to 1
                            while len(sample[z_tab[sample]==0]) != 0:
                                sample=custm.rvs(size=int(str(sys.argv[4])))
                            sample_lensz[k]=zgrid[sample]
                            #print "sample",sample
                            #print "sample_lensz1[k]", k, sample_lensz[k]
                            sample_massinf_tab=massinf_tab[sample] # since sample is constant, this insures that Mstar corresponds to z
                            sample_massmed_tab=massmed_tab[sample]
                            sample_masssup_tab=masssup_tab[sample]
                            #print "sample_massmed_tab", sample_massmed_tab
                            #print "sample_masssup_tab", sample_masssup_tab
                            sample_lenssup=np.zeros(int(str(sys.argv[4])))
                            sample_lensinf=np.zeros(int(str(sys.argv[4])))
                            for l in range(int(str(sys.argv[4]))):
                                sample_lenssup[l]=10**(sample_massmed_tab[l]+abs(np.random.normal(0, sample_masssup_tab[l]-sample_massmed_tab[l], 1)))
                                sample_lensinf[l]=10**(sample_massmed_tab[l]-abs(np.random.normal(0, sample_massmed_tab[l]-sample_massinf_tab[l], 1)))
                                rand=np.random.random(1)
                                sample_lensmass[k][l]=sample_lensinf[l]
                                frac_lensmass_samp=1-(1.0*(10**sample_masssup_tab[l]-10**sample_massmed_tab[l])/(10**sample_masssup_tab[l]-10**sample_massinf_tab[l]))
                                if rand>frac_lensmass_samp:
                                    sample_lensmass[k][l]=sample_lenssup[l]
                            #print "sample_lensmass1[k]", k, sample_lensmass[k]
                        else: # for the zpecz objects, since the input catalogue is properly edited so that "frac" is well-behaved, this section is a stripped-down version of "samp"
                            sample_lenssup=np.zeros(int(str(sys.argv[4])))
                            sample_lensinf=np.zeros(int(str(sys.argv[4])))
                            # in this case I sample from two gaussians, to the left and right of the most probable value
                            sample_lenssup=catupdate_lens[i][j][0]['z'][k]+abs(np.random.normal(0, catupdate_lens[i][j][0]['z_sup'][k]-catupdate_lens[i][j][0]['z'][k], int(str(sys.argv[4]))))
                            sample_lensinf=catupdate_lens[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['z'][k]-catupdate_lens[i][j][0]['z_inf'][k], int(str(sys.argv[4]))))
                            # make sure redshifts are positive
                            while len(sample_lensinf[sample_lensinf<0]) > 0:
                                sample_lensinf[sample_lensinf<0]=catupdate_lens[i][j][0]['z'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['z'][k]-catupdate_lens[i][j][0]['z_inf'][k], len(sample_lensinf[sample_lensinf<0])))
                            rand=np.random.random(int(str(sys.argv[4])))
                            sample_lensz[k]=sample_lensinf
                            sample_lensz[k][np.where(rand>frac_lensz[k])]=sample_lenssup[np.where(rand>frac_lensz[k])]
                            #print "sample_lensz2[k]", k, sample_lensz[k]
                            #sample_lensz[k][sample_lensz[k]<0]=0
                            sample_lenssup=10**(catupdate_lens[i][j][0]['mass_med'][k]+abs(np.random.normal(0, catupdate_lens[i][j][0]['mass_sup'][k]-catupdate_lens[i][j][0]['mass_med'][k], int(str(sys.argv[4])))))
                            sample_lensinf=10**(catupdate_lens[i][j][0]['mass_med'][k]-abs(np.random.normal(0, catupdate_lens[i][j][0]['mass_med'][k]-catupdate_lens[i][j][0]['mass_inf'][k], int(str(sys.argv[4])))))
                            rand=np.random.random(int(str(sys.argv[4])))
                            sample_lensmass[k]=sample_lensinf
                            sample_lensmass[k][np.where(rand>frac_lensmass[k])]=sample_lenssup[np.where(rand>frac_lensmass[k])]
                            #print "sample_lensmass2[k]", k, sample_lensmass[k]
                    for k in range(int(str(sys.argv[4]))):
                        # sum up the weights of different objects, corresponding to the same realization ID
                        # ignore objects with z>z_source for all weights when drawing from the z PDF
                        gal_lens_tab[i][j][k]=len(sample_lensz[:,k][sample_lensz[:,k]<z_s_lens])
                        oneoverr_lens_tab[i][j][k] = np.sum(1. / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zweight_lens_tab[i][j][k] = np.sum((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2))
                        mass_lens_tab[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens])
                        mass2_lens_tab[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2)
                        mass3_lens_tab[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**3)
                        mass2rms_lens_tab[i][j][k] = np.sqrt(mass2_lens[i][j][k])
                        mass3rms_lens_tab[i][j][k] = scipy.special.cbrt(mass3_lens[i][j][k])
                        zoverr_lens_tab[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        massoverr_lens_tab[i][j][k] = np.sum(sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens] / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass2overr_lens_tab[i][j][k] = np.sum((sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass3overr_lens_tab[i][j][k] = np.sum((sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**3) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zmassoverr_lens_tab[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) * (sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        zmass2overr_lens_tab[i][j][k] = np.sum(((z_s_lens * sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]) - (sample_lensz[:,k][sample_lensz[:,k]<z_s_lens]**2)) * (sample_lensmass[:,k][sample_lensz[:,k]<z_s_lens]**2) / catupdate_lens[i][j][0]['dist'][sample_lensz[:,k]<z_s_lens])
                        mass2overrrms_lens_tab[i][j][k] = np.sqrt(mass2overr_lens[i][j][k])
                        mass3overrrms_lens_tab[i][j][k] = scipy.special.cbrt(mass3overr_lens[i][j][k])

        print("--- %s seconds ---" % (time.time() - start_time))

#print gal_lens, sigma_mass2rms_lens

# I think RuntimeWarning: invalid value encountered in double_scalars can be ignored because it shows up only for rms weights, which means it's due to division by zero for the cells completely covered by the masks

# Initialize all field weight measurements as blank matrices

        #if str(sys.argv[8]) == "orig":
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
        #if str(sys.argv[8]) == "samp":
        field_gal = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_zweight = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_mass = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_mass2 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_mass3 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_zoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_massoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        field_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        #if str(sys.argv[8]) == "tab":
        tab_field_gal = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_zweight = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_mass = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_mass2 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_mass3 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_zoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_massoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_field_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))

# Use objects from the field catalogue that pass the redshift, separation, and mask tests

        print 'Computing field weights...'
        start_time = time.time()
        i=0
        pos_mstar=0
        #pos_mstar_remember=0
        with open('%s.cat' % [c[0:len(listfields[0])-1] for c in listfields][count]) as fieldcat:
            for gal in fieldcat:
                if (gal!="\n") and (gal.split()[0]!="#id"):
                    i=i+1
                    #print i
                    catfield=Table(names=('ID','ALPHA_J2000','DELTA_J2000','MASK','star_flag','MAG_i','MAG_y','Z_B','Z_B_MIN','Z_B_MAX','LP_log10_SM_MED','LP_log10_SM_INF','LP_log10_SM_SUP'),dtype=('<S13', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))
                    catfield.add_row([str(gal.split()[0]),gal.split()[4],gal.split()[5],gal.split()[36],gal.split()[60],gal.split()[79],gal.split()[84],gal.split()[39],gal.split()[40],gal.split()[41],gal.split()[61],gal.split()[62],gal.split()[63]])
                    #print catfield[0]
                    if i in [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000]:
                        print i, "objects..."
                    if (catfield['MASK'][0] == 0) and (catfield['star_flag'][0] == 0) and ((catfield['MAG_i'][0] <= float(str(sys.argv[6])) and catfield['MAG_i'][0] > 0) or (catfield['MAG_y'][0] <= float(str(sys.argv[6])) and catfield['MAG_y'][0] > 0)): # above threshold in either i or y bands
                      coordfieldwcs = SkyCoord(ra=catfield[0]['ALPHA_J2000']*u.degree, dec=catfield[0]['DELTA_J2000']*u.degree, frame='fk5')
                      #print coordfieldwcs.ra.deg, coordfieldwcs.dec.deg
                      #print i
                      with open('%s_pdz.cat' % [z[0:len(listfields[0])-1] for z in listfields][count]) as pdzcat:
                          galpdz=pdzcat.readlines()
                          linepdz=galpdz[i]
                      with open('%s_mstar.cat' % [z[0:len(listfields[0])-1] for z in listfields][count]) as mstarcat:
                          for m in xrange(pos_mstar):
                              mstarcat.next()
                          for galmstar in mstarcat:
                              pos_mstar=pos_mstar+1
                              if galmstar!="\n":
                                  if galmstar.split()[0]==catfield['ID'][0]:
                                      linemstar=galmstar
                                      break
                      # classification criteria:
                      where = (int(np.max([catfield['Z_B'][0]/0.05, 1])) * 6) # find the position in mstar.cat corresponding to z_b, in order to read chi_star and chi_gal
                      #print catfield['ID'][0],float(linemstar.split()[where-1]),float(linemstar.split()[where]),where
                      cond1 = str(sys.argv[7])=="old"
                      cond2 = (str(sys.argv[7])=="new") and (((catfield['MAG_i'][0] <=23) and (catfield['MAG_i'][0] >0)) or ((catfield['MAG_y'][0] <=23) and (catfield['MAG_y'][0] >0)))
                      if (where < 420) and (float(linemstar.split()[where])!=0):
                          cond3 = (str(sys.argv[7])=="new") and (float(linemstar.split()[where-1])/float(linemstar.split()[where])<0.5)
                      else:
                          cond3 = 1>2
                      if (where < 420) and (cond1 or cond2 or cond3):
                        #print linemstar.split()[0], linepdz.split()[0], catfield['ID'][0]
                        coordfieldpix = worldfield.wcs_world2pix(coordfieldwcs.ra.deg, coordfieldwcs.dec.deg, 0)
                        x = int((coordfieldpix[0] * pixCFHT) / (240.0 * u.arcsec))
                        y = int((coordfieldpix[1] * pixCFHT) / (240.0 * u.arcsec))
                        #print "x,y",x,y
                        cellpix_x = coordfieldpix[0] - (x * 240.0 * u.arcsec / pixCFHT)
                        cellpix_y = coordfieldpix[1] - (y * 240.0 * u.arcsec / pixCFHT)
                        sep = coordfieldwcs.separation(centerfields[x][y]).arcsec
                        if (maskedcell[x][y] >= 0.50):
                            if (sep >= 4) and (msk_lens[0].data[(cellpix_y * pixCFHT / pixlens).value][(cellpix_x * pixCFHT / pixlens).value] == 0): # ignore objects too close to the center and masked in the lens catalogue
                                #print sep
                                # fix anomalouss mass values in CFHTLENS. For "orig" and "samp" the point is to avoid unnecessary reading of the catalogues
                                #if (str(sys.argv[8]) == "orig") or (str(sys.argv[8]) == "samp"):
                                if catfield['LP_log10_SM_MED'][0] < 0:
                                    if float(linemstar.split()[where-3])>0: #mass_med
                                        catfield['LP_log10_SM_MED'][0]=linemstar.split()[where-3]
                                    elif float(linemstar.split()[where-5])>0: #mass_best
                                        catfield['LP_log10_SM_MED'][0]=linemstar.split()[where-5]
                                    else:
                                        pos=int(np.max([catfield['Z_B'][0]/0.05, 1]))
                                        while (pos<69) and (float(linemstar.split()[(pos * 6) - 5])<=0):
                                            pos=pos+1
                                        catfield['LP_log10_SM_MED'][0]=linemstar.split()[(pos * 6) - 5]
                                        if catfield['LP_log10_SM_MED'][0]<=0:
                                            pos=int(np.max([catfield['Z_B'][0]/0.05, 1]))
                                            while (pos>1) and (float(linemstar.split()[(pos * 6) - 5])<=0):
                                                pos=pos-1
                                            catfield['LP_log10_SM_MED'][0]=linemstar.split()[(pos * 6) - 5]
                                    if catfield['LP_log10_SM_MED'][0]<=0:
                                        catfield['LP_log10_SM_MED'][0]=9
                                #print "replaced mass: ", catfield['LP_log10_SM_MED'][0]
                                if (catfield['LP_log10_SM_INF'][0] < 0) or (catfield['LP_log10_SM_SUP'][0] < 0):
                                    catfield['LP_log10_SM_INF'][0]=catfield['LP_log10_SM_MED'][0]-0.1
                                    catfield['LP_log10_SM_SUP'][0]=catfield['LP_log10_SM_MED'][0]+0.1
                                if catfield['Z_B'][0] <= z_s_lens:
                                    origfield_gal[x][y] = origfield_gal[x][y] + 1
                                    origfield_zweight[x][y] = origfield_zweight[x][y] + (z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])
                                    origfield_mass[x][y] = origfield_mass[x][y] + 10**catfield['LP_log10_SM_MED'][0]
                                    origfield_mass2[x][y] = origfield_mass2[x][y] + ((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]))
                                    origfield_mass3[x][y] = origfield_mass3[x][y] + ((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]))
                                    if (sep <= 10):
                                        origfield_oneoverr[x][y] = origfield_oneoverr[x][y] + 0.1
                                        origfield_zoverr[x][y] = origfield_zoverr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) / 10)
                                        origfield_massoverr[x][y] = origfield_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][0]) / 10)
                                        origfield_mass2overr[x][y] = origfield_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0])) / 10)
                                        origfield_mass3overr[x][y] = origfield_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0])) / 10)
                                        origfield_zmassoverr[x][y] = origfield_zmassoverr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) * (10**catfield['LP_log10_SM_MED'][0]) / 10)
                                        origfield_zmass2overr[x][y] = origfield_zmass2overr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) * (10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]) / 10)
                                    else:
                                        origfield_oneoverr[x][y] = origfield_oneoverr[x][y] + 1. / sep
                                        origfield_zoverr[x][y] = origfield_zoverr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) / sep)
                                        origfield_massoverr[x][y] = origfield_massoverr[x][y] + ((10**catfield['LP_log10_SM_MED'][0]) / sep)
                                        origfield_mass2overr[x][y] = origfield_mass2overr[x][y] + (((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0])) / sep)
                                        origfield_mass3overr[x][y] = origfield_mass3overr[x][y] + (((10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0])) / sep)
                                        origfield_zmassoverr[x][y] = origfield_zmassoverr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) * (10**catfield['LP_log10_SM_MED'][0]) / sep)
                                        origfield_zmass2overr[x][y] = origfield_zmass2overr[x][y] + (((z_s_lens * catfield['Z_B'][0]) - (catfield['Z_B'][0] * catfield['Z_B'][0])) * (10**catfield['LP_log10_SM_MED'][0]) * (10**catfield['LP_log10_SM_MED'][0]) / sep)
                                    
                                samplesup=np.zeros(int(str(sys.argv[4])))
                                sampleinf=np.zeros(int(str(sys.argv[4])))
                                fracz=1-1.0*(catfield['Z_B_MAX'][0]-catfield['Z_B'][0])/(catfield['Z_B_MAX'][0]-catfield['Z_B_MIN'][0])
                                fracmass=1-(1.0*(10**catfield['LP_log10_SM_SUP'][0]-10**catfield['LP_log10_SM_MED'][0])/(10**catfield['LP_log10_SM_SUP'][0]-10**catfield['LP_log10_SM_INF'][0]))
                                #print "fracmass"
                                #print fracmass
                                #print catfield['LP_log10_SM_MED'][i],catfield['LP_log10_SM_SUP'][i],catfield['LP_log10_SM_INF'][i]

                                #if str(sys.argv[8]) == "samp":
                                if (fracz > 0) and (fracz < 1):
                                    samplesup=catfield['Z_B'][0]+abs(np.random.normal(0, catfield['Z_B_MAX'][0]-catfield['Z_B'][0], int(str(sys.argv[4]))))
                                    sampleinf=catfield['Z_B'][0]-abs(np.random.normal(0, catfield['Z_B'][0]-catfield['Z_B_MIN'][0], int(str(sys.argv[4]))))
                                        # no negative redshifts
                                    while len(sampleinf[sampleinf<0]) > 0:
                                        sampleinf[sampleinf<0]=catfield['Z_B'][0]-abs(np.random.normal(0, catfield['Z_B'][0]-catfield['Z_B_MIN'][0], len(sampleinf[sampleinf<0])))
                                    rand=np.random.random(int(str(sys.argv[4])))
                                    samplez=sampleinf
                                    samplez[np.where(rand>fracz)]=samplesup[np.where(rand>fracz)]
                                    #samplez[samplez<0]=0
                                if (fracz <= 0) or (fracz >= 1):
                                    samplez=np.random.normal(catfield['Z_B_MIN'][0]+(catfield['Z_B_MAX'][0]-catfield['Z_B_MIN'][0])/2, (catfield['Z_B_MAX'][0]-catfield['Z_B_MIN'][0])/2, int(str(sys.argv[4])))
                                    #samplez[samplez<0]=0
                                # for mass, assume gaussian distribution in log, not normal space, so I don't get things like negative mass
                                if (fracmass > 0) and (fracmass < 1):
                                    samplesup=10**(catfield['LP_log10_SM_MED'][0]+abs(np.random.normal(0, catfield['LP_log10_SM_SUP'][0]-catfield['LP_log10_SM_MED'][0], int(str(sys.argv[4])))))
                                    sampleinf=10**(catfield['LP_log10_SM_MED'][0]-abs(np.random.normal(0, catfield['LP_log10_SM_MED'][0]-catfield['LP_log10_SM_INF'][0], int(str(sys.argv[4])))))
                                    rand=np.random.random(int(str(sys.argv[4])))
                                    samplemass=sampleinf
                                    samplemass[np.where(rand>fracmass)]=samplesup[np.where(rand>fracmass)]
                                if (fracmass <= 0) or (fracmass >= 1):
                                    samplemass=10**np.random.normal(catfield['LP_log10_SM_INF'][0]+(catfield['LP_log10_SM_SUP'][0]-catfield['LP_log10_SM_INF'][0])/2, (catfield['LP_log10_SM_SUP'][0]-catfield['LP_log10_SM_INF'][0])/2, int(str(sys.argv[4])))
                                print samplemass
                                #ignore objects with z>z_source for all weights when drawing from the z PDF
                                field_gal[x][y][samplez<z_s_lens] = field_gal[x][y][samplez<z_s_lens] + 1
                                field_zweight[x][y][samplez<z_s_lens] = field_zweight[x][y][samplez<z_s_lens] + (z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)
                                field_mass[x][y][samplez<z_s_lens] = field_mass[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]
                                field_mass2[x][y][samplez<z_s_lens] = field_mass2[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]**2
                                field_mass3[x][y][samplez<z_s_lens] = field_mass3[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]**3
                                #print (catfield['Z_B_MAX'][i] - catfield['Z_B_MIN'][i])/2
                                if (sep <= 10):
                                    field_oneoverr[x][y][samplez<z_s_lens] = field_oneoverr[x][y][samplez<z_s_lens] + 0.1
                                    field_zoverr[x][y][samplez<z_s_lens] = field_zoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) / 10)
                                    field_massoverr[x][y][samplez<z_s_lens] = field_massoverr[x][y][samplez<z_s_lens] + (samplemass[samplez<z_s_lens] / 10)
                                    field_mass2overr[x][y][samplez<z_s_lens] = field_mass2overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**2) / 10)
                                    field_mass3overr[x][y][samplez<z_s_lens] = field_mass3overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**3) / 10)
                                    field_zmassoverr[x][y][samplez<z_s_lens] = field_zmassoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * samplemass[samplez<z_s_lens] / 10)
                                    field_zmass2overr[x][y][samplez<z_s_lens] = field_zmass2overr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * (samplemass[samplez<z_s_lens]**2) / 10)
                                else:
                                    field_oneoverr[x][y][samplez<z_s_lens] = field_oneoverr[x][y][samplez<z_s_lens] + 1/sep
                                    field_zoverr[x][y][samplez<z_s_lens] = field_zoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) / sep)
                                    field_massoverr[x][y][samplez<z_s_lens] = field_massoverr[x][y][samplez<z_s_lens] + (samplemass[samplez<z_s_lens] / sep)
                                    field_mass2overr[x][y][samplez<z_s_lens] = field_mass2overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**2) / sep)
                                    field_mass3overr[x][y][samplez<z_s_lens] = field_mass3overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**3) / sep)
                                    field_zmassoverr[x][y][samplez<z_s_lens] = field_zmassoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * samplemass[samplez<z_s_lens] / sep)
                                    field_zmass2overr[x][y][samplez<z_s_lens] = field_zmass2overr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * (samplemass[samplez<z_s_lens]**2) / sep)
                        
                                #if str(sys.argv[8]) == "tab":
                                massbest_tab=np.zeros(70)
                                massinf_tab=np.zeros(70)
                                massmed_tab=np.zeros(70)
                                masssup_tab=np.zeros(70)
                                pdz_tab=np.zeros(70)
                                for m in range(69):
                                    pdz_tab[m]=float(linepdz.split()[2+m][:-1]) # because the last character is comma
                                    if pdz_tab[m]<0.001:
                                        pdz_tab[m]=0
                                pdz_tab[69]=float(linepdz.split()[71])
                                if pdz_tab[69]<0.001:
                                    pdz_tab[69]=0
                                for m in range(70):
                                    if pdz_tab[m]!=0:
                                        massbest_tab[m]=float(linemstar.split()[1+6*m])
                                    if massbest_tab[m]<0:
                                        massbest_tab[m]=9 # very small number of exceptions
                                for m in range(70):
                                    if pdz_tab[m]!=0:
                                        massmed_tab[m]=float(linemstar.split()[3+6*m])
                                        if massmed_tab[m]<0:
                                            massmed_tab[m]=massbest_tab[m]
                                        massinf_tab[m]=float(linemstar.split()[2+6*m])
                                        masssup_tab[m]=float(linemstar.split()[4+6*m])
                                        if massinf_tab[m]<0:
                                            massinf_tab[m]=massmed_tab[m]-0.1
                                        if masssup_tab[m]<0:
                                            masssup_tab[m]=massmed_tab[m]+0.1
                                    #print ID,m,massbest_tab[m],massinf_tab[m],masssup_tab[m],masssup_tab[m]
                                custm = stats.rv_discrete(name='custm', values=(zgridint, pdz_tab))
                                sample=custm.rvs(size=int(str(sys.argv[4])))
                                iter=0
                                while len(sample[pdz_tab[sample]==0]) != 0: # happens because the probabilities do not sum exactly to 1; first reshuffle 10 times; if this does not solve the problem replace with the value having maximum probability; happens because the probabilities do not sum exactly to 1
                                    iter=iter+1
                                    sample=custm.rvs(size=int(str(sys.argv[4])))
                                    if iter==10:
                                        print pdz_tab,catfield['ID'][0],sample,str(sys.argv[1])
                                    sample[pdz_tab[sample]==0]=np.where(pdz_tab==np.max(pdz_tab[sample[pdz_tab[sample]!=0]]))[0][0]
                                    #print zgrid[sample],catfield['ID'][0]
                                samplez=zgrid[sample]
                                sample_massinf_tab=massinf_tab[sample] # since "sample" is constant, this insures that Mstar corresponds to z
                                sample_massmed_tab=massmed_tab[sample]
                                sample_masssup_tab=masssup_tab[sample]
                                sample_lenssup=np.zeros(int(str(sys.argv[4])))
                                sample_lensinf=np.zeros(int(str(sys.argv[4])))
                                #print sample_massmed_tab, samplez, catfield['ID'][0]
                                for l in range(int(str(sys.argv[4]))):
                                    if sample_massmed_tab[l]==0:     #I SHOULD NOT HAVE TO DO THIS, THERE IS A BUG
                                        sample_massmed_tab[l]=9
                                        sample_massinf_tab[l]=8.9
                                        sample_masssup_tab[l]=9.1
                                        print "Exception!", catfield['ID'][0]
                                    sample_lenssup[l]=10**(sample_massmed_tab[l]+abs(np.random.normal(0, sample_masssup_tab[l]-sample_massmed_tab[l], 1)))
                                    sample_lensinf[l]=10**(sample_massmed_tab[l]-abs(np.random.normal(0, sample_massmed_tab[l]-sample_massinf_tab[l], 1)))
                                    rand=np.random.random(1)
                                    samplemass[l]=sample_lensinf[l]
                                    fracmass=1-(1.0*(10**sample_masssup_tab[l]-10**sample_massmed_tab[l])/(10**sample_masssup_tab[l]-10**sample_massinf_tab[l]))
                                    if rand>fracmass:
                                        samplemass[l]=sample_lenssup[l]
                                #print catfield['ID'][0]
                                #print zgrid[sample]
                                #print samplemass[sample]
                                #print massinf_tab[sample]
                                #print massmed_tab[sample]
                                #print masssup_tab[sample]
                                #ignore objects with z>z_source for all weights when drawing from the z PDF
                                tab_field_gal[x][y][samplez<z_s_lens] = tab_field_gal[x][y][samplez<z_s_lens] + 1
                                tab_field_zweight[x][y][samplez<z_s_lens] = tab_field_zweight[x][y][samplez<z_s_lens] + (z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)
                                tab_field_mass[x][y][samplez<z_s_lens] = tab_field_mass[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]
                                tab_field_mass2[x][y][samplez<z_s_lens] = tab_field_mass2[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]**2
                                tab_field_mass3[x][y][samplez<z_s_lens] = tab_field_mass3[x][y][samplez<z_s_lens] + samplemass[samplez<z_s_lens]**3
                                #print (catfield['Z_B_MAX'][i] - catfield['Z_B_MIN'][i])/2
                                if (sep <= 10):
                                    tab_field_oneoverr[x][y][samplez<z_s_lens] = tab_field_oneoverr[x][y][samplez<z_s_lens] + 0.1
                                    tab_field_zoverr[x][y][samplez<z_s_lens] = tab_field_zoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) / 10)
                                    tab_field_massoverr[x][y][samplez<z_s_lens] = tab_field_massoverr[x][y][samplez<z_s_lens] + (samplemass[samplez<z_s_lens] / 10)
                                    tab_field_mass2overr[x][y][samplez<z_s_lens] = tab_field_mass2overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**2) / 10)
                                    tab_field_mass3overr[x][y][samplez<z_s_lens] = tab_field_mass3overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**3) / 10)
                                    tab_field_zmassoverr[x][y][samplez<z_s_lens] = tab_field_zmassoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * samplemass[samplez<z_s_lens] / 10)
                                    tab_field_zmass2overr[x][y][samplez<z_s_lens] = tab_field_zmass2overr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * (samplemass[samplez<z_s_lens]**2) / 10)
                                else:
                                    tab_field_oneoverr[x][y][samplez<z_s_lens] = tab_field_oneoverr[x][y][samplez<z_s_lens] + 1/sep
                                    tab_field_zoverr[x][y][samplez<z_s_lens] = tab_field_zoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) / sep)
                                    tab_field_massoverr[x][y][samplez<z_s_lens] = tab_field_massoverr[x][y][samplez<z_s_lens] + (samplemass[samplez<z_s_lens] / sep)
                                    tab_field_mass2overr[x][y][samplez<z_s_lens] = tab_field_mass2overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**2) / sep)
                                    tab_field_mass3overr[x][y][samplez<z_s_lens] = tab_field_mass3overr[x][y][samplez<z_s_lens] + ((samplemass[samplez<z_s_lens]**3) / sep)
                                    tab_field_zmassoverr[x][y][samplez<z_s_lens] = tab_field_zmassoverr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * samplemass[samplez<z_s_lens] / sep)
                                    tab_field_zmass2overr[x][y][samplez<z_s_lens] = tab_field_zmass2overr[x][y][samplez<z_s_lens] + (((z_s_lens * samplez[samplez<z_s_lens]) - (samplez[samplez<z_s_lens]**2)) * (samplemass[samplez<z_s_lens]**2) / sep)

        #if str(sys.argv[8]) == "orig":
        origfield_mass2rms = np.sqrt(origfield_mass2)
        origfield_mass3rms = scipy.special.cbrt(origfield_mass3)
        origfield_mass2overrrms = np.sqrt(origfield_mass2overr)
        origfield_mass3overrrms = scipy.special.cbrt(origfield_mass3overr)
        #if str(sys.argv[8]) == "samp":
        field_mass2rms = np.sqrt(field_mass2)
        field_mass3rms = scipy.special.cbrt(field_mass3)
        field_mass2overrrms = np.sqrt(field_mass2overr)
        field_mass3overrrms = scipy.special.cbrt(field_mass3overr)
        #if str(sys.argv[8]) == "tab":
        tab_field_mass2rms = np.sqrt(tab_field_mass2)
        tab_field_mass3rms = scipy.special.cbrt(tab_field_mass3)
        tab_field_mass2overrrms = np.sqrt(tab_field_mass2overr)
        tab_field_mass3overrrms = scipy.special.cbrt(tab_field_mass3overr)

        #if str(sys.argv[8]) == "orig":
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
        #if str(sys.argv[8]) == "samp":
        q_gal = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_zweight = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass2 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass3 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_zoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_massoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass2rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass3rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass2overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        q_mass3overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        #if str(sys.argv[8]) == "tab":
        tab_q_gal = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_oneoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_zweight = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass2 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass3 = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_zoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_massoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass3overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_zmassoverr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_zmass2overr = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass2rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass3rms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass2overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))
        tab_q_mass3overrrms = np.zeros((cells_on_a_side,cells_on_a_side,int(str(sys.argv[4]))))

        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50: # skipping the cells with large masked area
                    #if str(sys.argv[8]) == "orig":
                    origq_gal[i][j] = gal_lensorig[i][j] * 1.0 / origfield_gal[i][j]
                    origq_zweight[i][j] = zweight_lensorig[i][j] * 1.0 / origfield_zweight[i][j]
                    origq_mass[i][j] = mass_lensorig[i][j] * 1.0 / origfield_mass[i][j]
                    origq_mass2[i][j] = mass2_lensorig[i][j] * 1.0 / origfield_mass2[i][j]
                    origq_mass3[i][j] = mass3_lensorig[i][j] * 1.0 / origfield_mass3[i][j]
                    origq_oneoverr[i][j] = oneoverr_lensorig[i][j] * 1.0 / origfield_oneoverr[i][j]
                    origq_zoverr[i][j] = zoverr_lensorig[i][j] * 1.0 / origfield_zoverr[i][j]
                    origq_massoverr[i][j] = massoverr_lensorig[i][j] * 1.0 / origfield_massoverr[i][j]
                    origq_mass2overr[i][j] = mass2overr_lensorig[i][j] * 1.0 / origfield_mass2overr[i][j]
                    origq_mass3overr[i][j] = mass3overr_lensorig[i][j] * 1.0 / origfield_mass3overr[i][j]
                    origq_zmassoverr[i][j] = zmassoverr_lensorig[i][j] * 1.0 / origfield_zmassoverr[i][j]
                    origq_zmass2overr[i][j] = zmass2overr_lensorig[i][j] * 1.0 / origfield_zmass2overr[i][j]
                    origq_mass2rms[i][j] = mass2rms_lensorig[i][j] * 1.0 / origfield_mass2rms[i][j]
                    origq_mass3rms[i][j] = mass3rms_lensorig[i][j] * 1.0 / origfield_mass3rms[i][j]
                    origq_mass2overrrms[i][j] = mass2overrrms_lensorig[i][j] * 1.0 / origfield_mass2overrrms[i][j]
                    origq_mass3overrrms[i][j] = mass3overrrms_lensorig[i][j] * 1.0 / origfield_mass3overrrms[i][j]
                    #if str(sys.argv[8]) == "samp":
                    q_gal[i][j] = gal_lens[i][j] * 1.0 / field_gal[i][j]
                    q_oneoverr[i][j] = oneoverr_lens[i][j] * 1.0 / field_oneoverr[i][j]
                    q_zweight[i][j] = zweight_lens[i][j] * 1.0 / field_zweight[i][j]
                    q_mass[i][j] = mass_lens[i][j] * 1.0 / field_mass[i][j]
                    q_mass2[i][j] = mass2_lens[i][j] * 1.0 / field_mass2[i][j]
                    q_mass3[i][j] = mass3_lens[i][j] * 1.0 / field_mass3[i][j]
                    q_zoverr[i][j] = zoverr_lens[i][j] * 1.0 / field_zoverr[i][j]
                    q_massoverr[i][j] = massoverr_lens[i][j] * 1.0 / field_massoverr[i][j]
                    q_mass2overr[i][j] = mass2overr_lens[i][j] * 1.0 / field_mass2overr[i][j]
                    q_mass3overr[i][j] = mass3overr_lens[i][j] * 1.0 / field_mass3overr[i][j]
                    q_zmassoverr[i][j] = zmassoverr_lens[i][j] * 1.0 / field_zmassoverr[i][j]
                    q_zmass2overr[i][j] = zmass2overr_lens[i][j] * 1.0 / field_zmass2overr[i][j]
                    q_mass2rms[i][j] = mass2rms_lens[i][j] * 1.0 / field_mass2rms[i][j]
                    q_mass3rms[i][j] = mass3rms_lens[i][j] * 1.0 / field_mass3rms[i][j]
                    q_mass2overrrms[i][j] = mass2overrrms_lens[i][j] * 1.0 / field_mass2overrrms[i][j]
                    q_mass3overrrms[i][j] = mass3overrrms_lens[i][j] * 1.0 / field_mass3overrrms[i][j]
                    #if str(sys.argv[8]) == "tab":
                    tab_q_gal[i][j] = gal_lens_tab[i][j] * 1.0 / tab_field_gal[i][j]
                    tab_q_oneoverr[i][j] = oneoverr_lens_tab[i][j] * 1.0 / tab_field_oneoverr[i][j]
                    tab_q_zweight[i][j] = zweight_lens_tab[i][j] * 1.0 / tab_field_zweight[i][j]
                    tab_q_mass[i][j] = mass_lens_tab[i][j] * 1.0 / tab_field_mass[i][j]
                    tab_q_mass2[i][j] = mass2_lens_tab[i][j] * 1.0 / tab_field_mass2[i][j]
                    tab_q_mass3[i][j] = mass3_lens_tab[i][j] * 1.0 / tab_field_mass3[i][j]
                    tab_q_zoverr[i][j] = zoverr_lens_tab[i][j] * 1.0 / tab_field_zoverr[i][j]
                    tab_q_massoverr[i][j] = massoverr_lens_tab[i][j] * 1.0 / tab_field_massoverr[i][j]
                    tab_q_mass2overr[i][j] = mass2overr_lens_tab[i][j] * 1.0 / tab_field_mass2overr[i][j]
                    tab_q_mass3overr[i][j] = mass3overr_lens_tab[i][j] * 1.0 / tab_field_mass3overr[i][j]
                    tab_q_zmassoverr[i][j] = zmassoverr_lens_tab[i][j] * 1.0 / tab_field_zmassoverr[i][j]
                    tab_q_zmass2overr[i][j] = zmass2overr_lens_tab[i][j] * 1.0 / tab_field_zmass2overr[i][j]
                    tab_q_mass2rms[i][j] = mass2rms_lens_tab[i][j] * 1.0 / tab_field_mass2rms[i][j]
                    tab_q_mass3rms[i][j] = mass3rms_lens_tab[i][j] * 1.0 / tab_field_mass3rms[i][j]
                    tab_q_mass2overrrms[i][j] = mass2overrrms_lens_tab[i][j] * 1.0 / tab_field_mass2overrrms[i][j]
                    tab_q_mass3overrrms[i][j] = mass3overrrms_lens_tab[i][j] * 1.0 / tab_field_mass3overrrms[i][j]


        print("--- %s seconds ---" % (time.time() - start_time))

        print ("No. of fields with area without mask > 75 percent & > 50 percent: %d/%d, %d/%d" % ((maskedcell>=0.75).sum(), cells_on_a_side ** 2, (maskedcell>=0.50).sum(), cells_on_a_side ** 2))

# write a file where the most probable value of z and Mstar has been used, and a file for each of the 100 or 1000 samples
        address="/Volumes/G-RAIDStudio/CFHTlens/forhistograms"
        f = open('%s/%s_%s_q50_orig_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        #if str(sys.argv[8]) == "samp":
        g = open('%s/%s_%s_q50_samp_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        #if str(sys.argv[8]) == "tab":
        h = open('%s/%s_%s_q50_tab_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.50:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (origq_gal[i][j], origq_oneoverr[i][j], origq_zweight[i][j], origq_mass[i][j], origq_mass2[i][j], origq_mass2rms[i][j], origq_mass3[i][j], origq_mass3rms[i][j], origq_zoverr[i][j], origq_massoverr[i][j], origq_mass2overr[i][j], origq_mass3overr[i][j], origq_mass2overrrms[i][j], origq_mass3overrrms[i][j], origq_zmassoverr[i][j], origq_zmass2overr[i][j]))
                    for k in range(int(str(sys.argv[4]))):
                        g.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j][k], q_oneoverr[i][j][k], q_zweight[i][j][k], q_mass[i][j][k], q_mass2[i][j][k], q_mass2rms[i][j][k], q_mass3[i][j][k], q_mass3rms[i][j][k], q_zoverr[i][j][k], q_massoverr[i][j][k], q_mass2overr[i][j][k], q_mass3overr[i][j][k], q_mass2overrrms[i][j][k], q_mass3overrrms[i][j][k], q_zmassoverr[i][j][k], q_zmass2overr[i][j][k]))
                    for k in range(int(str(sys.argv[4]))):
                        h.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (tab_q_gal[i][j][k], tab_q_oneoverr[i][j][k], tab_q_zweight[i][j][k], tab_q_mass[i][j][k], tab_q_mass2[i][j][k], tab_q_mass2rms[i][j][k], tab_q_mass3[i][j][k], tab_q_mass3rms[i][j][k], tab_q_zoverr[i][j][k], tab_q_massoverr[i][j][k], tab_q_mass2overr[i][j][k], tab_q_mass3overr[i][j][k], tab_q_mass2overrrms[i][j][k], tab_q_mass3overrrms[i][j][k], tab_q_zmassoverr[i][j][k], tab_q_zmass2overr[i][j][k]))
        f.close()
        #if str(sys.argv[8]) == "samp":
        g.close()
        #if str(sys.argv[8]) == "tab":
        h.close()

        f = open('%s/%s_%s_q75_orig_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        #if str(sys.argv[8]) == "samp":
        g = open('%s/%s_%s_q75_samp_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        #if str(sys.argv[8]) == "tab":
        h = open('%s/%s_%s_q75_tab_size%s_i%s_%s.lst' % (address,[d[38:len(listfields[0])-1] for d in listfields][count],str(sys.argv[1]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7])),'w')
        for i in range(cells_on_a_side):
            for j in range(cells_on_a_side):
                if maskedcell[i][j] >= 0.75:
                    f.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (origq_gal[i][j], origq_oneoverr[i][j], origq_zweight[i][j], origq_mass[i][j], origq_mass2[i][j], origq_mass2rms[i][j], origq_mass3[i][j], origq_mass3rms[i][j], origq_zoverr[i][j], origq_massoverr[i][j], origq_mass2overr[i][j], origq_mass3overr[i][j], origq_mass2overrrms[i][j], origq_mass3overrrms[i][j], origq_zmassoverr[i][j], origq_zmass2overr[i][j]))
                    for k in range(int(str(sys.argv[4]))):
                        g.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (q_gal[i][j][k], q_oneoverr[i][j][k], q_zweight[i][j][k], q_mass[i][j][k], q_mass2[i][j][k], q_mass2rms[i][j][k], q_mass3[i][j][k], q_mass3rms[i][j][k], q_zoverr[i][j][k], q_massoverr[i][j][k], q_mass2overr[i][j][k], q_mass3overr[i][j][k], q_mass2overrrms[i][j][k], q_mass3overrrms[i][j][k], q_zmassoverr[i][j][k], q_zmass2overr[i][j][k]))
                    for k in range(int(str(sys.argv[4]))):
                        h.write('q_gal= %.4e q_oneoverr= %.4e q_zweight= %.4e q_mass= %.4e q_mass2= %.4e q_mass2rms= %.4e q_mass3= %.4e q_mass3rms= %.4e q_zoverr= %.4e q_massoverr= %.4e q_mass2overr= %.4e q_mass3overr= %.4e q_mass2overrrms= %.4e q_mass3overrrms= %.4e q_zmassoverr= %.4e q_zmass2overr= %.4e \n' % (tab_q_gal[i][j][k], tab_q_oneoverr[i][j][k], tab_q_zweight[i][j][k], tab_q_mass[i][j][k], tab_q_mass2[i][j][k], tab_q_mass2rms[i][j][k], tab_q_mass3[i][j][k], tab_q_mass3rms[i][j][k], tab_q_zoverr[i][j][k], tab_q_massoverr[i][j][k], tab_q_mass2overr[i][j][k], tab_q_mass3overr[i][j][k], tab_q_mass2overrrms[i][j][k], tab_q_mass3overrrms[i][j][k], tab_q_zmassoverr[i][j][k], tab_q_zmass2overr[i][j][k]))

        f.close()
        #mstarcat_lens.close()
        #if str(sys.argv[8]) == "samp":
        g.close()
        #if str(sys.argv[8]) == "tab":
        h.close()

        print("Total time for subfield: --- %s seconds ---" % (time.time() - start_timesubfield))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'

