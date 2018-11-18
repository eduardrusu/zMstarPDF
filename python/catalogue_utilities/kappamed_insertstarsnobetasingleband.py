# CE Rusu, Feb 12 2018
# This code uses the Millenium Sumilation convergence and shear maps as well as the associated SA catalogue of galaxies, in order to compute the weighted counts for fields centered around each kappa and gamma point. This is done for a variety of limiting magnitudes, aperture radii, and weights.
# run with the following arguments: lens name, field name, limiting mag, outer mask radius, type, inner mask radius, zinf, zsup (in case I remove redshift slices); e.g.: python /lfs08/rusucs/code/kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
# the code can oly be used for limmag 23 or 24 currently

import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
import pandas as pd
import time
#import distances
from scipy.interpolate import griddata
import astropy.table as table
from astropy.io import fits # for tables

############################
# function definitions
############################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...") # use for debugging

def readbinary(replacestr):
    replace = plane + replacestr
    os.system("sed \"11s/.*/  const char kappa_file_name[]   = \\\"\%s\\\";/\" readKappaBinary.c > readKappaBinary_%s_%s_%s_%s_%s_%sinner.c_" % (replace,lens,type,plane,float(limmag),radius,innermsk))
    os.system("sed \"35s/.*/  fpt = fopen (\\\"kappa_values_%s_%s_%s_%s_%s_%sinner.dat\\\", \\\"w\\\");/\"  readKappaBinary_%s_%s_%s_%s_%s_%sinner.c_ > readKappaBinary_%s_%s_%s_%s_%s_%sinner.c" % (lens,type,plane,float(limmag),radius,innermsk,lens,type,plane,float(limmag),radius,innermsk,lens,type,plane,float(limmag),radius,innermsk))
    os.system("rm -f readKappaBinary_%s_%s_%s_%s_%s_%sinner.c_" % (lens,type,plane,float(limmag),radius,innermsk))
    os.system("gcc readKappaBinary_%s_%s_%s_%s_%s_%sinner.c -o compiled_%s_%s_%s_%s_%s_%sinner.out" % (lens,type,plane,float(limmag),radius,innermsk,lens,type,plane,float(limmag),radius,innermsk))
    os.system("./compiled_%s_%s_%s_%s_%s_%sinner.out" % (lens,type,plane,float(limmag),radius,innermsk))
    os.system("rm -f readKappaBinary_%s_%s_%s_%s_%s_%sinner.c" % (lens,type,plane,float(limmag),radius,innermsk))
    os.system("rm -f compiled_%s_%s_%s_%s_%s_%sinner.out" % (lens,type,plane,float(limmag),radius,innermsk))

def contaminants(count,cont_ugr,posxmin,posxmax,posymin,posymax,star_rmag):
    cont = np.random.random_integers(0,499,int(count * cont_ugr)) # randomly select from the star catalogues which contain 500 stars each
    cont_posx = np.random.uniform(posxmin,posxmax,int(count * cont_ugr))
    cont_posy = np.random.uniform(posymin,posymax,int(count * cont_ugr))
    cont_rmag = star_rmag[cont]
    return cont, cont_posx, cont_posy, cont_rmag

def weightedcounts(cat,spacing,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,bands):
    initialized = 0
    for i in range(spacing):
        for j in range(spacing):
            print "(%s,%s)/(%s,%s) for radius %s" % (i+1,j+1,spacing,spacing,radius)
            grid_x, grid_y = np.mgrid[lim1D + i:4096 - lim1D - (4096 - 2 * lim1D) % spacing - spacing + i:complex(0,cells_on_a_side), lim1D + j:4096 - lim1D - (4096 - 2 * lim1D) % spacing - spacing + j:complex(0,cells_on_a_side)] # the grid containing the kappa pixel at the center of cells
            cellx = grid_x.flatten()
            celly = grid_y.flatten()
            cellxy = np.array([cellx,celly])
            posxy = np.array([-0.5 * L_field  + (1 + cellxy[0] + 0.5) * L_pix, -0.5 * L_field  + (1 + cellxy[1] + 0.5) * L_pix])
            cellkappagamma = np.c_[cells,kappagamma[:,][(cellxy[0] * 4096 + cellxy[1]).astype(int)]]
            index = griddata(posxy.T, cells, (cat[:,index_posx], cat[:,index_posy]), method='nearest')
            sep = np.sqrt((posxy[0][index.astype(int)]-cat[:,index_posx])**2 + (posxy[1][index.astype(int)]-cat[:,index_posy])**2)*1/degree*3600
            cat_msk = np.c_[cat,index,sep]
            catinner = cat_msk[cat_msk[:,index_sep] <= innermsk] # so that I can count how many objects are inside the inner mask
            cat_msk = cat_msk[cat_msk[:,index_sep] <= radius] # mask objects at distance larger than the aperture from the center
            cat_msk = cat_msk[cat_msk[:,index_sep] > innermsk] # uses the inner mask
            cat_msk[:,index_sep][cat_msk[:,index_sep] < 10] = 10
            #cat = np.c_[z,posx,posy,mstar,rmag...]

            w_gal_2X = np.bincount(cat_msk[:,index_index].astype(int)) # 2X stands for 23 or 24 limmag
            galinner = np.bincount(catinner[:,index_index].astype(int)) # counts objects inside the inner mask
            index_all = np.unique(index.astype(int)) # sorts all unique entries; galinner and w_gal_2X may contain fewer elements, because in some indexed cells (it seems largest index only) there may be no galaxies, so I need to expand them to include all indices
            #print len(index_all),len(w_gal_2X),len(galinner)
            try:
                w_gal_2X = np.append(w_gal_2X,np.zeros(len(index_all)-len(w_gal_2X))) # in rare cases, for the 45" aperture, index will have missing fields (it will miss one of the integers from 0 to len(index)). So this assignment will not work, because np.zeros(len(index_all)-len(w_gal_2X))is a negative number
                galinner = np.append(galinner,np.zeros(len(index_all)-len(galinner)))
            except:
                missing = []
                for k in range(len(index_all)):
                    if k not in index_all: missing = np.append(missing,k)
                for k in range(len(missing)):
                    insert = np.copy(cat_msk[0]) # any entry would do
                    insert[index_index] = missing[k]
                    cat_msk = np.append(cat_msk,insert.reshape(1,5),axis = 0)
                index_all = np.append(index_all,missing)
                w_gal_2X = np.append(w_gal_2X,np.zeros(len(index_all)-len(w_gal_2X)))
                galinner = np.append(galinner,np.zeros(len(index_all)-len(galinner)))

            try:
                p_oneoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'oneoverr':1.0 / cat_msk[:,index_sep]})
                w_oneoverr_2X = p_oneoverr.groupby(['cell']).median().values[:,0] * w_gal_2X # this might fail for radius=45, where there are the larger fluctuations in the number of galaxies, and as I remove galaxies from cat_msk there might be an index for which all cells contain zero galaxies. In that case the index is removed from cat_msk, but _zweight.groupby(['cell']).median().values[:,0] needs all indices to be present. The solution is to insert a ghost line into cat_msk for each missing index
            except:
                #print len(w_gal_2X[w_gal_2X==0])
                missing = np.array([])
                cat_msk_unique = np.unique(cat_msk[:,index_index]).astype(int) # to speed up the search
                for k in range(int(np.max(index_all))):
                    if k not in cat_msk_unique:
                        missing = np.append(missing,np.array([k]))
                for k in range(len(missing)):
                    insert = np.copy(cat_msk[0]) # any entry would do
                    insert[index_index] = missing[k]
                    cat_msk=np.append(cat_msk,insert.reshape(1,5),axis = 0)

            p_oneoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'oneoverr':1.0 / cat_msk[:,index_sep]})
            #for k in range(len(missing)):
                #cat_msk = np.delete(cat_msk,-1,axis = 0) # delete the last line I inserted above; actually this is not necessary because cat_msk is no longer used
            w_oneoverr_2X = p_oneoverr.groupby(['cell']).median().values[:,0] * w_gal_2X

            cellkappagamma = np.c_[cellkappagamma,w_gal_2X,w_oneoverr_2X,galinner]
            cellkappagammastyle = np.c_[cellkappagamma[:,1].astype(int),np.around(cellkappagamma[:,2],decimals=5),np.around(cellkappagamma[:,3],decimals=5),np.around(cellkappagamma[:,4],decimals=5),cellkappagamma[:,5].astype(int),np.around(cellkappagamma[:,6],decimals=4),cellkappagamma[:,7].astype(int)]
            if initialized != 0:
            	cellkappagammafinal = np.r_[cellkappagammafinal,cellkappagammastyle]
            else:
                f = '%snobeta%s%smedinject_%s_%s_%s_%s_%s_%sarcsecinner.fits' % (rootwghtratios,pln,type,bands,lens,plane[0:13],limmag,radius,innermsk)
                os.system('rm -f %s' % f)
                cellkappagammafinal = cellkappagammastyle
                initialized = 1
            if (i == spacing - 1) and (j == spacing - 1):
                tableout = table.Table(cellkappagammafinal, names=('ID', 'kappa', 'gamma1', 'gamma2', 'w_gal_%s' % limmag, 'w_oneoverr_%s' % limmag, 'galinner_%s' % limmag), dtype=(np.int32,np.float32,np.float32,np.float32,np.int32,np.float32,np.float32))
                #fits.append(f, tableout.as_array())
                tableout.write(f)
                del tableout

############################
# lens information
############################

start_time = time.time()

lens = str(sys.argv[1])
plane = str(sys.argv[2])
limmag = float(str(sys.argv[3]))
radiusstr = str(sys.argv[4])
type = str(sys.argv[5]) # computed or measured
innermsk = int(str(sys.argv[6])) # inner mask in arcsec

if lens == "PG1115":
    brightmag = 15.3
    #limmag = 22.5
    pln = 34
    if (radiusstr == "45"):
        hstcoverage = 0
        radius = 45
        fracspec20 = 1 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.73
        fracspec23 = 0.15
        fracspec24 = 0.00
    if (radiusstr == "120"):
        hstcoverage = 0
        radius = 120
        fracspec20 = 0.69
        fracspec21 = 0.81
        fracspec22 = 0.52
        fracspec23 = 0.07
        fracspec24 = 0.00

rootwghtratios = "/lfs08/rusucs/%s/MSwghtratios/" % lens
#rootwghtratios = "/u/flashscratch/c/cerusu/MSwghtratios/"
#rootwghtratios = "/mnt/scratch/rusucs/%s/MSwghtratios/" % lens
#rootwghtratios = "/Volumes/LaCieSubaru/MSweights/"
rootgals = "/lfs08/rusucs/%s/MSgals/" % lens
#rootgals = "/u/flashscratch/c/cerusu/MSgals/"
#rootgals = "/mnt/scratch/rusucs/%s/MSgals/" % lens
#rootgals = "/Volumes/LaCieSubaru/MSgals/"
rootkappaplanes = "/lfs08/rusucs/kappaplanes/"
#rootkappaplanes = "/u/flashscratch/c/cerusu/kappaplanes/"
#rootkappaplanes = "/mnt/scratch/rusucs/kappaplanes/"
#rootkappaplanes = "/Volumes/LaCieSubaru/kappaplanes/"
rootstars = "/lfs08/rusucs/insertstars/"
#rootstars = "/u/flashscratch/c/cerusu/insertstars/"
#rootstars = "/Volumes/LaCieSubaru/insertstars/"
#rootstars = "/mnt/scratch/rusucs/insertstars/"

# contamination and incompleteness based on Figure 9 W1 from Hildebrandt 2012

cont_h12_18 = 0.00
cont_h12_185 = 0.12
cont_h12_19 = 0.08
cont_h12_195 = 0.03
cont_h12_20 = 0.04
cont_h12_205 = 0.06
cont_h12_21 = 0.05
cont_h12_215 = 0.02
cont_h12_22 = 0.01
cont_h12_225 = 0.02
cont_h12_23 = 0.03
cont_h12_235 = 0.01
cont_h12_24 = 0.01

inc_h12_18 = 0.20
inc_h12_185 = 0.13
inc_h12_19 = 0.20
inc_h12_195 = 0.00
inc_h12_20 = 0.03
inc_h12_205 = 0.02
inc_h12_21 = 0.01
inc_h12_215 = 0.07
inc_h12_22 = 0.05
inc_h12_225 = 0.05
inc_h12_23 = 0.03
inc_h12_235 = 0.02
inc_h12_24 = 0.01

nospec_18 = 1 - fracspec20
nospec_185 = 1 - fracspec20
nospec_19 = 1 - fracspec20
nospec_195 = 1 - fracspec20
nospec_20 = 1 - fracspec20
nospec_205 = 1 - fracspec21
nospec_21 = 1 - fracspec21
nospec_215 = 1 - fracspec22
nospec_22 = 1 - fracspec22
nospec_225 = 1 - fracspec23
nospec_23 = 1 - fracspec23
nospec_235 = 1 - fracspec24
nospec_24 = 1 - fracspec24

cont_ugrizJHK_18 = cont_h12_18 * nospec_18 * (1 - hstcoverage)
cont_ugrizJHK_185 = cont_h12_185 * nospec_185 * (1 - hstcoverage)
cont_ugrizJHK_19 = cont_h12_19 * nospec_19 * (1 - hstcoverage)
cont_ugrizJHK_195 = cont_h12_195 * nospec_195 * (1 - hstcoverage)
cont_ugrizJHK_20 = cont_h12_20 * nospec_20 * (1 - hstcoverage)
cont_ugrizJHK_205 = cont_h12_205 * nospec_205 * (1 - hstcoverage)
cont_ugrizJHK_21 = cont_h12_21 * nospec_21 * (1 - hstcoverage)
cont_ugrizJHK_215 = cont_h12_215 * nospec_215 * (1 - hstcoverage)
cont_ugrizJHK_22 = cont_h12_22 * nospec_22 * (1 - hstcoverage)
cont_ugrizJHK_225 = cont_h12_225 * nospec_225 * (1 - hstcoverage)
cont_ugrizJHK_23 = cont_h12_23 * nospec_23 * (1 - hstcoverage)
cont_ugrizJHK_235 = cont_h12_235 * nospec_235 * (1 - hstcoverage)
cont_ugrizJHK_24 = cont_h12_24 * nospec_24 * (1 - hstcoverage)

cont_ugriz_18 = cont_h12_18
cont_ugriz_185 = cont_h12_185
cont_ugriz_19 = cont_h12_19
cont_ugriz_195 = cont_h12_195
cont_ugriz_20 = cont_h12_20
cont_ugriz_205 = cont_h12_205
cont_ugriz_21 = cont_h12_21
cont_ugriz_215 = cont_h12_215
cont_ugriz_22 = cont_h12_22
cont_ugriz_225 = cont_h12_225
cont_ugriz_23 = cont_h12_23
cont_ugriz_235 = cont_h12_235
cont_ugriz_24 = cont_h12_24

inc_ugrizJHK_18 = inc_h12_18 * nospec_18 * (1 - hstcoverage)
inc_ugrizJHK_185 = inc_h12_185 * nospec_185 * (1 - hstcoverage)
inc_ugrizJHK_19 = inc_h12_19 * nospec_19 * (1 - hstcoverage)
inc_ugrizJHK_195 = inc_h12_195 * nospec_195 * (1 - hstcoverage)
inc_ugrizJHK_20 = inc_h12_20 * nospec_20 * (1 - hstcoverage)
inc_ugrizJHK_205 = inc_h12_205 * nospec_205 * (1 - hstcoverage)
inc_ugrizJHK_21 = inc_h12_21 * nospec_21 * (1 - hstcoverage)
inc_ugrizJHK_215 = inc_h12_215 * nospec_215 * (1 - hstcoverage)
inc_ugrizJHK_22 = inc_h12_22 * nospec_22 * (1 - hstcoverage)
inc_ugrizJHK_225 = inc_h12_225 * nospec_225 * (1 - hstcoverage)
inc_ugrizJHK_23 = inc_h12_23 * nospec_23 * (1 - hstcoverage)
inc_ugrizJHK_235 = inc_h12_235 * nospec_235 * (1 - hstcoverage)
inc_ugrizJHK_24 = inc_h12_24 * nospec_24 * (1 - hstcoverage)

inc_ugriz_18 = inc_h12_18
inc_ugriz_185 = inc_h12_185
inc_ugriz_19 = inc_h12_19
inc_ugriz_195 = inc_h12_195
inc_ugriz_20 = inc_h12_20
inc_ugriz_205 = inc_h12_205
inc_ugriz_21 = inc_h12_21
inc_ugriz_215 = inc_h12_215
inc_ugriz_22 = inc_h12_22
inc_ugriz_225 = inc_h12_225
inc_ugriz_23 = inc_h12_23
inc_ugriz_235 = inc_h12_235
inc_ugriz_24 = inc_h12_24

# read the stellar contaminants I will insert

# !!!!!!!!!!! the mags are actually i-band, not r-band, but I will ignore that since I'm not 'estimating' z and Mstar
star_rmag_18 = np.loadtxt("%sstar018zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_185 = np.loadtxt("%sstar18185zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_19 = np.loadtxt("%sstar18519zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_195 = np.loadtxt("%sstar19195zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_20 = np.loadtxt("%sstar19520zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_205 = np.loadtxt("%sstar20205zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_21 = np.loadtxt("%sstar20521zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_215 = np.loadtxt("%sstar21215zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_22 = np.loadtxt("%sstar21522zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_225 = np.loadtxt("%sstar22225zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_23 = np.loadtxt("%sstar22523zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_235 = np.loadtxt("%sstar23235zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)
star_rmag_24 = np.loadtxt("%sstar23524zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0], unpack=True)

############################
# read the binary kappa file
############################

start_readkappa = time.time()

if str(pln) in plane:
    os.chdir(rootkappaplanes)
    readbinary(".kappa")
    pos1D,kappa = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner.dat" % (lens,type,plane,float(limmag),radius,innermsk), unpack=True)
    readbinary(".gamma_1")
    gamma1 = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner.dat" % (lens,type,plane,float(limmag),radius,innermsk), usecols = [1], unpack=True)
    readbinary(".gamma_2")
    gamma2 = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner.dat" % (lens,type,plane,float(limmag),radius,innermsk), usecols = [1], unpack=True)
    kappagamma = np.c_[pos1D,kappa,gamma1,gamma2]
    os.system("rm -f kappa_values_%s_%s_%s_%s_%s_%sinner.dat" % (lens,type,plane,float(limmag),radius,innermsk))
else: sys.exit('Wrong MS plane for this lens!!!')

############################
# prepare the galaxy catalogues
############################

start_readcat = time.time()

posx_ugrizJHK = np.array([])
posy_ugrizJHK = np.array([])
rmag_ugrizJHK = np.array([])

posx_ugriz = np.array([])
posy_ugriz = np.array([])
rmag_ugriz = np.array([])

root = plane[0:13]

for i in range(4):
    for j in range(4):
        file_ugrizJHK = '%s%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_sampledr.images.txt' % (rootgals,root,i,j)
        file_ugriz = '%s%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_originalr.images.txt' % (rootgals,root,i,j)
        if "measured" in type:
            posx__ugrizJHK, posy__ugrizJHK, rmag__ugrizJHK = np.loadtxt(file_ugrizJHK, usecols = (2,3,4), unpack=True)
            posx__ugriz, posy__ugriz, rmag__ugriz = np.loadtxt(file_ugriz, usecols = (2,3,4), unpack=True)
        elif "computed" in type:
            posx__ugrizJHK, posy__ugrizJHK, rmag__ugrizJHK = np.loadtxt(file_ugriz, usecols = (2,3,4), unpack=True)
            posx__ugriz, posy__ugriz, rmag__ugriz = np.loadtxt(file_ugriz, usecols = (2,3,4), unpack=True)
        else: sys.exit('Only \"measured\" and \"computed\" are accepted as the fourth argument of the code. Execution stopped.')

        if "measured" in type:
            # count the galaxies in order to implement the fraction of stars
            count18 = posx__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= 18)].size
            count185 = posx__ugrizJHK[(rmag__ugrizJHK > 18) & (rmag__ugrizJHK <= 18.5)].size
            count19 = posx__ugrizJHK[(rmag__ugrizJHK > 18.5) & (rmag__ugrizJHK <= 19)].size
            count195 = posx__ugrizJHK[(rmag__ugrizJHK > 19) & (rmag__ugrizJHK <= 19.5)].size
            count20 = posx__ugrizJHK[(rmag__ugrizJHK > 19.5) & (rmag__ugrizJHK <= 20)].size
            count205 = posx__ugrizJHK[(rmag__ugrizJHK > 20) & (rmag__ugrizJHK <= 20.5)].size
            count21 = posx__ugrizJHK[(rmag__ugrizJHK > 20.5) & (rmag__ugrizJHK <= 21)].size
            count215 = posx__ugrizJHK[(rmag__ugrizJHK > 21) & (rmag__ugrizJHK <= 21.5)].size
            count22 = posx__ugrizJHK[(rmag__ugrizJHK > 21.5) & (rmag__ugrizJHK <= 22)].size
            count225 = posx__ugrizJHK[(rmag__ugrizJHK > 22) & (rmag__ugrizJHK <= 22.5)].size
            count23 = posx__ugrizJHK[(rmag__ugrizJHK > 22.5) & (rmag__ugrizJHK <= 23)].size
            count235 = posx__ugrizJHK[(rmag__ugrizJHK > 23) & (rmag__ugrizJHK <= 23.5)].size
            count24 = posx__ugrizJHK[(rmag__ugrizJHK > 23.5) & (rmag__ugrizJHK <= 24)].size
            #print count18,count185,count19,count195,count20,count205,count21,count215,count22,count225,count23,count235,count24

            posxmin = np.min(posx__ugrizJHK) # used to insert random stellar contaminants
            posxmax = np.max(posx__ugrizJHK)
            posymin = np.min(posy__ugrizJHK)
            posymax = np.max(posy__ugrizJHK)

            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_rmag_18 = contaminants(count18,cont_ugrizJHK_18,posxmin,posxmax,posymin,posymax,star_rmag_18)
            cont_185,cont_posx_185,cont_posy_185,cont_rmag_185 = contaminants(count185,cont_ugrizJHK_185,posxmin,posxmax,posymin,posymax,star_rmag_185)
            cont_19,cont_posx_19,cont_posy_19,cont_rmag_19 = contaminants(count19,cont_ugrizJHK_19,posxmin,posxmax,posymin,posymax,star_rmag_19)
            cont_195,cont_posx_195,cont_posy_195,cont_rmag_195 = contaminants(count195,cont_ugrizJHK_195,posxmin,posxmax,posymin,posymax,star_rmag_195)
            cont_20,cont_posx_20,cont_posy_20,cont_rmag_20 = contaminants(count20,cont_ugrizJHK_20,posxmin,posxmax,posymin,posymax,star_rmag_20)
            cont_205,cont_posx_205,cont_posy_205,cont_rmag_205 = contaminants(count205,cont_ugrizJHK_205,posxmin,posxmax,posymin,posymax,star_rmag_205)
            cont_21,cont_posx_21,cont_posy_21,cont_rmag_21 = contaminants(count21,cont_ugrizJHK_21,posxmin,posxmax,posymin,posymax,star_rmag_21)
            cont_215,cont_posx_215,cont_posy_215,cont_rmag_215 = contaminants(count215,cont_ugrizJHK_215,posxmin,posxmax,posymin,posymax,star_rmag_215)
            cont_22,cont_posx_22,cont_posy_22,cont_rmag_22 = contaminants(count22,cont_ugrizJHK_22,posxmin,posxmax,posymin,posymax,star_rmag_22)
            cont_225,cont_posx_225,cont_posy_225,cont_rmag_225 = contaminants(count225,cont_ugrizJHK_225,posxmin,posxmax,posymin,posymax,star_rmag_225)
            cont_23,cont_posx_23,cont_posy_23,cont_rmag_23 = contaminants(count23,cont_ugrizJHK_23,posxmin,posxmax,posymin,posymax,star_rmag_23)
            cont_235,cont_posx_235,cont_posy_235,cont_rmag_235 = contaminants(count235,cont_ugrizJHK_235,posxmin,posxmax,posymin,posymax,star_rmag_235)
            cont_24,cont_posx_24,cont_posy_24,cont_rmag_24 = contaminants(count24,cont_ugrizJHK_24,posxmin,posxmax,posymin,posymax,star_rmag_24)
            #print cont_18.size,cont_185.size,cont_19.size,cont_195.size,cont_20.size,cont_205.size,cont_21.size,cont_215.size,cont_22.size,cont_225.size,cont_23.size,cont_235.size,cont_24.size

            # masking the fraction of galaxies expected to not be detected as galaxies, due to incompleteness; here also apply the brightmag, limmag, z_s and z gap cuts
            posx__ugrizJHK = posx__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag) & (((rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= 18)) | ((rmag__ugrizJHK > 18) & (rmag__ugrizJHK <= 18.5)) | ((rmag__ugrizJHK > 18.5) & (rmag__ugrizJHK <= 19)) | ((rmag__ugrizJHK > 19) & (rmag__ugrizJHK <= 19.5)) | ((rmag__ugrizJHK > 19.5) & (rmag__ugrizJHK <= 20)) | ((rmag__ugrizJHK > 20) & (rmag__ugrizJHK <= 20.5)) | ((rmag__ugrizJHK > 20.5) & (rmag__ugrizJHK <= 21)) | ((rmag__ugrizJHK > 21) & (rmag__ugrizJHK <= 21.5)) | ((rmag__ugrizJHK > 21.5) & (rmag__ugrizJHK <= 22)) | ((rmag__ugrizJHK > 22) & (rmag__ugrizJHK <= 22.5)) | ((rmag__ugrizJHK > 22.5) & (rmag__ugrizJHK <= 23)) | ((rmag__ugrizJHK > 23) & (rmag__ugrizJHK <= 23.5)) | ((rmag__ugrizJHK > 23.5) & (rmag__ugrizJHK <= 24)))]
            posy__ugrizJHK = posy__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag) & (((rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= 18)) | ((rmag__ugrizJHK > 18) & (rmag__ugrizJHK <= 18.5)) | ((rmag__ugrizJHK > 18.5) & (rmag__ugrizJHK <= 19)) | ((rmag__ugrizJHK > 19) & (rmag__ugrizJHK <= 19.5)) | ((rmag__ugrizJHK > 19.5) & (rmag__ugrizJHK <= 20)) | ((rmag__ugrizJHK > 20) & (rmag__ugrizJHK <= 20.5)) | ((rmag__ugrizJHK > 20.5) & (rmag__ugrizJHK <= 21)) | ((rmag__ugrizJHK > 21) & (rmag__ugrizJHK <= 21.5)) | ((rmag__ugrizJHK > 21.5) & (rmag__ugrizJHK <= 22)) | ((rmag__ugrizJHK > 22) & (rmag__ugrizJHK <= 22.5)) | ((rmag__ugrizJHK > 22.5) & (rmag__ugrizJHK <= 23)) | ((rmag__ugrizJHK > 23) & (rmag__ugrizJHK <= 23.5)) | ((rmag__ugrizJHK > 23.5) & (rmag__ugrizJHK <= 24)))]
            rmag__ugrizJHK = rmag__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag) & (((rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= 18)) | ((rmag__ugrizJHK > 18) & (rmag__ugrizJHK <= 18.5)) | ((rmag__ugrizJHK > 18.5) & (rmag__ugrizJHK <= 19)) | ((rmag__ugrizJHK > 19) & (rmag__ugrizJHK <= 19.5)) | ((rmag__ugrizJHK > 19.5) & (rmag__ugrizJHK <= 20)) | ((rmag__ugrizJHK > 20) & (rmag__ugrizJHK <= 20.5)) | ((rmag__ugrizJHK > 20.5) & (rmag__ugrizJHK <= 21)) | ((rmag__ugrizJHK > 21) & (rmag__ugrizJHK <= 21.5)) | ((rmag__ugrizJHK > 21.5) & (rmag__ugrizJHK <= 22)) | ((rmag__ugrizJHK > 22) & (rmag__ugrizJHK <= 22.5)) | ((rmag__ugrizJHK > 22.5) & (rmag__ugrizJHK <= 23)) | ((rmag__ugrizJHK > 23) & (rmag__ugrizJHK <= 23.5)) | ((rmag__ugrizJHK > 23.5) & (rmag__ugrizJHK <= 24)))]

            # inserting the stellar contaminants
            if limmag == 22.5:
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225))
                rmag_ugrizJHK = np.concatenate((rmag_ugrizJHK,rmag__ugrizJHK,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225))
            if limmag == 23:
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23))
                rmag_ugrizJHK = np.concatenate((rmag_ugrizJHK,rmag__ugrizJHK,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23))
            if limmag == 23.5:
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235))
                rmag_ugrizJHK = np.concatenate((rmag_ugrizJHK,rmag__ugrizJHK,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23,cont_rmag_235))
            if limmag == 24:
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235,cont_posx_24))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235,cont_posy_24))
                rmag_ugrizJHK = np.concatenate((rmag_ugrizJHK,rmag__ugrizJHK,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23,cont_rmag_235,cont_rmag_24))
            #if len(rmag_ugrizJHK) > 0: print np.max(rmag_ugrizJHK)

            # repeat for ugriz catalogues

            count18 = posx__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= 18)].size
            count185 = posx__ugriz[(rmag__ugriz > 18) & (rmag__ugriz <= 18.5)].size
            count19 = posx__ugriz[(rmag__ugriz > 18.5) & (rmag__ugriz <= 19)].size
            count195 = posx__ugriz[(rmag__ugriz > 19) & (rmag__ugriz <= 19.5)].size
            count20 = posx__ugriz[(rmag__ugriz > 19.5) & (rmag__ugriz <= 20)].size
            count205 = posx__ugriz[(rmag__ugriz > 20) & (rmag__ugriz <= 20.5)].size
            count21 = posx__ugriz[(rmag__ugriz > 20.5) & (rmag__ugriz <= 21)].size
            count215 = posx__ugriz[(rmag__ugriz > 21) & (rmag__ugriz <= 21.5)].size
            count22 = posx__ugriz[(rmag__ugriz > 21.5) & (rmag__ugriz <= 22)].size
            count225 = posx__ugriz[(rmag__ugriz > 22) & (rmag__ugriz <= 22.5)].size
            count23 = posx__ugriz[(rmag__ugriz > 22.5) & (rmag__ugriz <= 23)].size
            count235 = posx__ugriz[(rmag__ugriz > 23) & (rmag__ugriz <= 23.5)].size
            count24 = posx__ugriz[(rmag__ugriz > 23.5) & (rmag__ugriz <= 24)].size

            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_rmag_18 = contaminants(count18,cont_ugriz_18,posxmin,posxmax,posymin,posymax,star_rmag_18)
            cont_185,cont_posx_185,cont_posy_185,cont_rmag_185 = contaminants(count185,cont_ugriz_185,posxmin,posxmax,posymin,posymax,star_rmag_185)
            cont_19,cont_posx_19,cont_posy_19,cont_rmag_19 = contaminants(count19,cont_ugriz_19,posxmin,posxmax,posymin,posymax,star_rmag_19)
            cont_195,cont_posx_195,cont_posy_195,cont_rmag_195 = contaminants(count195,cont_ugriz_195,posxmin,posxmax,posymin,posymax,star_rmag_195)
            cont_20,cont_posx_20,cont_posy_20,cont_rmag_20 = contaminants(count20,cont_ugriz_20,posxmin,posxmax,posymin,posymax,star_rmag_20)
            cont_205,cont_posx_205,cont_posy_205,cont_rmag_205 = contaminants(count205,cont_ugriz_205,posxmin,posxmax,posymin,posymax,star_rmag_205)
            cont_21,cont_posx_21,cont_posy_21,cont_rmag_21 = contaminants(count21,cont_ugriz_21,posxmin,posxmax,posymin,posymax,star_rmag_21)
            cont_215,cont_posx_215,cont_posy_215,cont_rmag_215 = contaminants(count215,cont_ugriz_215,posxmin,posxmax,posymin,posymax,star_rmag_215)
            cont_22,cont_posx_22,cont_posy_22,cont_rmag_22 = contaminants(count22,cont_ugriz_22,posxmin,posxmax,posymin,posymax,star_rmag_22)
            cont_225,cont_posx_225,cont_posy_225,cont_rmag_225 = contaminants(count225,cont_ugriz_225,posxmin,posxmax,posymin,posymax,star_rmag_225)
            cont_23,cont_posx_23,cont_posy_23,cont_rmag_23 = contaminants(count23,cont_ugriz_23,posxmin,posxmax,posymin,posymax,star_rmag_23)
            cont_235,cont_posx_235,cont_posy_235,cont_rmag_235 = contaminants(count235,cont_ugriz_235,posxmin,posxmax,posymin,posymax,star_rmag_235)
            cont_24,cont_posx_24,cont_posy_24,cont_rmag_24 = contaminants(count24,cont_ugriz_24,posxmin,posxmax,posymin,posymax,star_rmag_24)

            posx__ugriz = posx__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag) & (((rmag__ugriz > brightmag) & (rmag__ugriz <= 18)) | ((rmag__ugriz > 18) & (rmag__ugriz <= 18.5)) | ((rmag__ugriz > 18.5) & (rmag__ugriz <= 19)) | ((rmag__ugriz > 19) & (rmag__ugriz <= 19.5)) | ((rmag__ugriz > 19.5) & (rmag__ugriz <= 20)) | ((rmag__ugriz > 20) & (rmag__ugriz <= 20.5)) | ((rmag__ugriz > 20.5) & (rmag__ugriz <= 21)) | ((rmag__ugriz > 21) & (rmag__ugriz <= 21.5)) | ((rmag__ugriz > 21.5) & (rmag__ugriz <= 22)) | ((rmag__ugriz > 22) & (rmag__ugriz <= 22.5)) | ((rmag__ugriz > 22.5) & (rmag__ugriz <= 23)) | ((rmag__ugriz > 23) & (rmag__ugriz <= 23.5)) | ((rmag__ugriz > 23.5) & (rmag__ugriz <= 24)))]
            posy__ugriz = posy__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag) & (((rmag__ugriz > brightmag) & (rmag__ugriz <= 18)) | ((rmag__ugriz > 18) & (rmag__ugriz <= 18.5)) | ((rmag__ugriz > 18.5) & (rmag__ugriz <= 19)) | ((rmag__ugriz > 19) & (rmag__ugriz <= 19.5)) | ((rmag__ugriz > 19.5) & (rmag__ugriz <= 20)) | ((rmag__ugriz > 20) & (rmag__ugriz <= 20.5)) | ((rmag__ugriz > 20.5) & (rmag__ugriz <= 21)) | ((rmag__ugriz > 21) & (rmag__ugriz <= 21.5)) | ((rmag__ugriz > 21.5) & (rmag__ugriz <= 22)) | ((rmag__ugriz > 22) & (rmag__ugriz <= 22.5)) | ((rmag__ugriz > 22.5) & (rmag__ugriz <= 23)) | ((rmag__ugriz > 23) & (rmag__ugriz <= 23.5)) | ((rmag__ugriz > 23.5) & (rmag__ugriz <= 24)))]
            rmag__ugriz = rmag__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag) & (((rmag__ugriz > brightmag) & (rmag__ugriz <= 18)) | ((rmag__ugriz > 18) & (rmag__ugriz <= 18.5)) | ((rmag__ugriz > 18.5) & (rmag__ugriz <= 19)) | ((rmag__ugriz > 19) & (rmag__ugriz <= 19.5)) | ((rmag__ugriz > 19.5) & (rmag__ugriz <= 20)) | ((rmag__ugriz > 20) & (rmag__ugriz <= 20.5)) | ((rmag__ugriz > 20.5) & (rmag__ugriz <= 21)) | ((rmag__ugriz > 21) & (rmag__ugriz <= 21.5)) | ((rmag__ugriz > 21.5) & (rmag__ugriz <= 22)) | ((rmag__ugriz > 22) & (rmag__ugriz <= 22.5)) | ((rmag__ugriz > 22.5) & (rmag__ugriz <= 23)) | ((rmag__ugriz > 23) & (rmag__ugriz <= 23.5)) | ((rmag__ugriz > 23.5) & (rmag__ugriz <= 24)))]

            if limmag == 22.5:
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225))
                rmag_ugriz = np.concatenate((rmag_ugriz,rmag__ugriz,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225))
            if limmag == 23:
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23))
                rmag_ugriz = np.concatenate((rmag_ugriz,rmag__ugriz,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23))
            if limmag == 23.5:
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235))
                rmag_ugriz = np.concatenate((rmag_ugriz,rmag__ugriz,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23,cont_rmag_235))
            if limmag == 24:
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235,cont_posx_24))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235,cont_posy_24))
                rmag_ugriz = np.concatenate((rmag_ugriz,rmag__ugriz,cont_rmag_18,cont_rmag_185,cont_rmag_19,cont_rmag_195,cont_rmag_20,cont_rmag_205,cont_rmag_21,cont_rmag_215,cont_rmag_22,cont_rmag_225,cont_rmag_23,cont_rmag_235,cont_rmag_24))
            #if len(rmag_ugriz) > 0: print np.max(rmag_ugriz)

        else:
            posx__ugrizJHK = posx__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag)]
            posy__ugrizJHK = posy__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag)]
            rmag__ugrizJHK = rmag__ugrizJHK[(rmag__ugrizJHK > brightmag) & (rmag__ugrizJHK <= limmag)]
            posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK))
            posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK))
            rmag_ugrizJHK = np.concatenate((rmag_ugrizJHK,rmag__ugrizJHK))

            posx__ugriz = posx__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag)]
            posy__ugriz = posy__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag)]
            rmag__ugriz = rmag__ugriz[(rmag__ugriz > brightmag) & (rmag__ugriz <= limmag)]
            posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz))
            posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz))
            rmag_ugriz = np.concatenate((rmag_ugriz,rmag__ugriz))

cat_ugrizJHK = np.c_[posx_ugrizJHK,posy_ugrizJHK,rmag_ugrizJHK]
del posx_ugrizJHK
del posy_ugrizJHK
del rmag_ugrizJHK
del posx__ugrizJHK
del posy__ugrizJHK
del rmag__ugrizJHK
cat_ugriz = np.c_[posx_ugriz,posy_ugriz,rmag_ugriz]
del posx_ugriz
del posy_ugriz
del rmag_ugriz
del posx__ugriz
del posy__ugriz
del rmag__ugriz

index_posx = 0
index_posy = 1
index_rmag = 2
index_sep = -1
index_index = -2

print "Read galaxy catalogues in ", time.time() - start_readcat, "seconds"

############################
# compute and write weighted counts
############################

degree = np.pi / 180
L_field = 4.0 * degree
N_pix_per_dim = 4096
L_pix = L_field / N_pix_per_dim
pixscl_asec = 4.0 * 3600 / 4096 # kappa pixel size in arcsec
pixscl_rad = pixscl_asec * degree / 3600 # kappa pixel size in radian

# divide the field into cells and compute weighted counts
start_weights = time.time()

lim1D = int(radius / pixscl_asec) + 1 # number of kappa pixels that should be ignored at the edge of the field so that the field can be covered by full cells
spacing = int(2.0 * radius / pixscl_asec + 1) # number of pixels between each pixels considered as cell center; in order for the cells not to overlap, the grid spacing should be at least this large
cells_on_a_side = int(1.0 * (4096 - 2 * lim1D) / spacing)
cells = np.linspace(0,cells_on_a_side**2 - 1,cells_on_a_side**2)
start_radius = time.time()

cat = cat_ugrizJHK
bands = "ugrizJHK"
weightedcounts(cat,spacing,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,bands)

cat = cat_ugriz
bands = "ugriz"
weightedcounts(cat,spacing,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,bands)

print "Computed weights in ", time.time() - start_weights, "seconds"

print(" Field done in --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
