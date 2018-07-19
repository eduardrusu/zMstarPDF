# CE Rusu, Feb 12 2018
# This code uses the Millenium Sumilation convergence and shear maps as well as the associated SA catalogue of galaxies, in order to compute the weighted counts for fields centered around each kappa and gamma point. This is done for a variety of limiting magnitudes, aperture radii, and weights.
# run with the following arguments: lens name, field name, limiting mag, outer mask radius, type, inner mask radius, zinf, zsup (in case I remove redshift slices); e.g.: python /lfs08/rusucs/code/kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1.0 -1.0
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
    os.system("sed \"11s/.*/  const char kappa_file_name[]   = \\\"\%s\\\";/\" readKappaBinary.c > readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c_" % (replace,lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("sed \"35s/.*/  fpt = fopen (\\\"kappa_values_%s_%s_%s_%s_%s_%sinner_%s.dat\\\", \\\"w\\\");/\"  readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c_ > readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c" % (lens,type,plane,int(limmag),radius,innermsk,gap,lens,type,plane,int(limmag),radius,innermsk,gap,lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("rm -f readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c_" % (lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("gcc readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c -o compiled_%s_%s_%s_%s_%s_%sinner_%s.out" % (lens,type,plane,int(limmag),radius,innermsk,gap,lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("./compiled_%s_%s_%s_%s_%s_%sinner_%s.out" % (lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("rm -f readKappaBinary_%s_%s_%s_%s_%s_%sinner_%s.c" % (lens,type,plane,int(limmag),radius,innermsk,gap))
    os.system("rm -f compiled_%s_%s_%s_%s_%s_%sinner_%s.out" % (lens,type,plane,int(limmag),radius,innermsk,gap))

def contaminants(count,cont_ugr,posxmin,posxmax,posymin,posymax,star_imag,star_z):
    cont = np.random.random_integers(0,499,int(count * cont_ugr)) # randomly select from the star catalogues which contain 500 stars each
    cont_posx = np.random.uniform(posxmin,posxmax,int(count * cont_ugr))
    cont_posy = np.random.uniform(posymin,posymax,int(count * cont_ugr))
    cont_imag = star_imag[cont]
    cont_z = star_z[cont]
    cont = cont[((cont_z >= zsup) | (cont_z <= zinf))] # remove the redshift gap
    cont_posx = cont_posx[((cont_z >= zsup) | (cont_z <= zinf))]
    cont_posy = cont_posy[((cont_z >= zsup) | (cont_z <= zinf))]
    cont_imag = cont_imag[((cont_z >= zsup) | (cont_z <= zinf))]
    cont_z = cont_z[((cont_z >= zsup) | (cont_z <= zinf))]
    return cont, cont_posx, cont_posy, cont_imag, cont_z

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
            #cat = np.c_[z,posx,posy,mstar,imag...]
            
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
                    cat_msk = np.append(cat_msk,insert.reshape(1,6),axis = 0)
                index_all = np.append(index_all,missing)
                w_gal_2X = np.append(w_gal_2X,np.zeros(len(index_all)-len(w_gal_2X)))
                galinner = np.append(galinner,np.zeros(len(index_all)-len(galinner)))
            
            try:
                p_zweight = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'zweight':1.0 * (z_s * cat_msk[:,index_z]) - cat_msk[:,index_z]**2})
                w_zweight_2X = p_zweight.groupby(['cell']).median().values[:,0] * w_gal_2X # this might fail for radius=45, where there are the larger fluctuations in the number of galaxies, and as I remove galaxies from cat_msk there might be an index for which all cells contain zero galaxies. In that case the index is removed from cat_msk, but _zweight.groupby(['cell']).median().values[:,0] needs all indices to be present. The solution is to insert a ghost line into cat_msk for each missing index
            except:
                missing = np.array([])
                cat_msk_unique = np.unique(cat_msk[:,index_index]).astype(int) # to speed up the search
                for k in range(int(np.max(index_all))):
                    if k not in cat_msk_unique:
                        missing = np.append(missing,np.array([k]))
                for k in range(len(missing)):
                    insert = np.copy(cat_msk[0]) # any entry would do
                    insert[index_index] = missing[k]
                    cat_msk=np.append(cat_msk,insert.reshape(1,6),axis = 0)
            
            p_zweight = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'zweight':1.0 * (z_s * cat_msk[:,index_z]) - cat_msk[:,index_z]**2})
            p_zoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'zoverr':1.0 * ((z_s * cat_msk[:,index_z]) - cat_msk[:,index_z]**2) / cat_msk[:,index_sep]})
            p_oneoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'oneoverr':1.0 / cat_msk[:,index_sep]})
            
            w_zweight_2X = p_zweight.groupby(['cell']).median().values[:,0] * w_gal_2X
            w_oneoverr_2X = p_oneoverr.groupby(['cell']).median().values[:,0] * w_gal_2X
            w_zoverr_2X = p_zoverr.groupby(['cell']).median().values[:,0] * w_gal_2X
            
            cellkappagamma = np.c_[cellkappagamma,w_gal_2X,w_zweight_2X,w_oneoverr_2X,w_zoverr_2X,galinner]
            cellkappagammastyle = np.c_[cellkappagamma[:,1].astype(int),np.around(cellkappagamma[:,2],decimals=5),np.around(cellkappagamma[:,3],decimals=5),np.around(cellkappagamma[:,4],decimals=5),cellkappagamma[:,5].astype(int),np.around(cellkappagamma[:,6],decimals=4),np.around(cellkappagamma[:,7],decimals=4),np.around(cellkappagamma[:,8],decimals=4),cellkappagamma[:,9].astype(int)]
            if initialized != 0:
                cellkappagammafinal = np.r_[cellkappagammafinal,cellkappagammastyle]
            else:
                f = '%snobeta%s%smedinject_%s_%s_%s_%s_%s_%sarcsecinner_%s.fits' % (rootwghtratios,pln,type,bands,lens,plane[0:13],int(limmag),radius,innermsk,gap)
                os.system('rm -f %s' % f)
                cellkappagammafinal = cellkappagammastyle
                initialized = 1
            if (i == spacing - 1) and (j == spacing - 1):
                tableout = table.Table(cellkappagammafinal, names=('ID', 'kappa', 'gamma1', 'gamma2', 'w_gal_%s' % limmag, 'w_zweight_%s' % limmag, 'w_oneoverr_%s' % limmag, 'w_zoverr_%s' % limmag, 'galinner_%s' % limmag), dtype=(np.int32,np.float32,np.float32,np.float32,np.int32,np.float32,np.float32,np.float32,np.int32))
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
zinf = float(str(sys.argv[7]))
zsup = float(str(sys.argv[8]))
gap = 'gap_%s_%s' % (zinf,zsup)

if lens == "B1608":
    z_s = 1.39
    pln = 37 # MS simulation plane
if lens == "HE0435":
    z_s = 1.69
    z_l = 0.455
    brightmag = 17.50
    #pln = 34 & 35
    if (radiusstr == "45"):
        radius = 45
        fracspec20 = 1 # gal+stars
        fracspec21 = 0.75
        fracspec22 = 0.83
        fracspec23 = 0.11
        fracspec24 = 0.00
    if (radiusstr == "120"):
        radius = 120
        fracspec20 = 1
        fracspec21 = 0.87
        fracspec22 = 0.44
        fracspec23 = 0.08
        fracspec24 = 0.00
    if (radiusstr == "60"):
        radius = 60
        fracspec20 = 1 # gal+stars
        fracspec21 = 0.8
        fracspec22 = 0.78
        fracspec23 = 0.15
        fracspec24 = 0.00
    if (radiusstr == "90"):
        radius = 90
        fracspec20 = 1
        fracspec21 = 0.85
        fracspec22 = 0.5
        fracspec23 = 0.1
        fracspec24 = 0.00
if lens == "WFI2033":
    z_s = 1.66
    z_l = 0.66
    brightmag = 16.90
    limmag = 23
    pln = 35
    if (radiusstr == "45"):
        hstcoverage = 1
        radius = 45
        fracspec20 = 1 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.73
        fracspec23 = 0.15
        fracspec24 = 0.00
    if (radiusstr == "120"):
        hstcoverage = 1 * 0.47
        radius = 120
        fracspec20 = 0.69
        fracspec21 = 0.81
        fracspec22 = 0.52
        fracspec23 = 0.07
        fracspec24 = 0.00
    if (radiusstr == "60"):
        hstcoverage = 1
        radius = 60
        fracspec20 = 1 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.70
        fracspec23 = 0.11
        fracspec24 = 0.00
    if (radiusstr == "90"):
        hstcoverage = 1 * 0.47 * (90.0**2)/(120.0**2)
        radius = 90
        fracspec20 = 0.71
        fracspec21 = 0.88
        fracspec22 = 0.6
        fracspec23 = 0.06
        fracspec24 = 0.00
if lens == "HE1104":
    z_s = 2.32
#pln = 30 & 31
if lens == "RX1131":
    z_s = 0.66
if lens == "J1206":
    z_s = 1.79
    z_l = 0.745
    pln = 34
    #limmag = 24
    brightmag = 18.00
    if (radiusstr == "45"):
        hstcoverage = 1
        radius = 45
        fracspec20 = 0.40 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.83
        fracspec23 = 0.29
        fracspec24 = 0.13
    if (radiusstr == "120"):
        hstcoverage = 1 * 0.44
        radius = 120
        fracspec20 = 0.40
        fracspec21 = 0.81
        fracspec22 = 0.65
        fracspec23 = 0.19
        fracspec24 = 0.04

#rootwghtratios = "/lfs08/rusucs/%s/MSwghtratios/" % lens
#rootwghtratios = "/u/flashscratch/c/cerusu/MSwghtratios/"
#rootwghtratios = "/mnt/scratch/rusucs/%s/MSwghtratios/" % lens
rootwghtratios = "/Volumes/LaCieSubaru/MSweights/"
#rootgals = "/lfs08/rusucs/%s/MSgals/" % lens
#rootgals = "/u/flashscratch/c/cerusu/MSgals/"
#rootgals = "/mnt/scratch/rusucs/%s/MSgals/" % lens
rootgals = "/Volumes/LaCieSubaru/MSgals/"
#rootkappaplanes = "/lfs08/rusucs/kappaplanes/"
#rootkappaplanes = "/u/flashscratch/c/cerusu/kappaplanes/"
#rootkappaplanes = "/mnt/scratch/rusucs/kappaplanes/"
rootkappaplanes = "/Volumes/LaCieSubaru/kappaplanes/"
#rootstars = "/lfs08/rusucs/insertstars/"
#rootstars = "/u/flashscratch/c/cerusu/insertstars/"
rootstars = "/Volumes/LaCieSubaru/insertstars/"
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

star_imag_18, star_z_18 = np.loadtxt("%sstar018zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_185, star_z_185 = np.loadtxt("%sstar18185zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_19, star_z_19 = np.loadtxt("%sstar18519zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_195, star_z_195 = np.loadtxt("%sstar19195zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_20, star_z_20 = np.loadtxt("%sstar19520zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_205, star_z_205 = np.loadtxt("%sstar20205zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_21, star_z_21 = np.loadtxt("%sstar20521zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_215, star_z_215 = np.loadtxt("%sstar21215zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_22, star_z_22 = np.loadtxt("%sstar21522zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_225, star_z_225 = np.loadtxt("%sstar22225zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_23, star_z_23 = np.loadtxt("%sstar22523zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_235, star_z_235 = np.loadtxt("%sstar23235zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)
star_imag_24, star_z_24 = np.loadtxt("%sstar23524zcut_catpdzmstar_magrepaired.cat" % rootstars, usecols = [0,2], unpack=True)

############################
# read the binary kappa file
############################

start_readkappa = time.time()

if str(pln) in plane:
    os.chdir(rootkappaplanes)
    readbinary(".kappa")
    pos1D,kappa = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner_%s.dat" % (lens,type,plane,int(limmag),radius,innermsk,gap), unpack=True)
    readbinary(".gamma_1")
    gamma1 = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner_%s.dat" % (lens,type,plane,int(limmag),radius,innermsk,gap), usecols = [1], unpack=True)
    readbinary(".gamma_2")
    gamma2 = np.loadtxt("kappa_values_%s_%s_%s_%s_%s_%sinner_%s.dat" % (lens,type,plane,int(limmag),radius,innermsk,gap), usecols = [1], unpack=True)
    kappagamma = np.c_[pos1D,kappa,gamma1,gamma2]
    os.system("rm -f kappa_values_%s_%s_%s_%s_%s_%sinner_%s.dat" % (lens,type,plane,int(limmag),radius,innermsk,gap))
else: sys.exit('Wrong MS plane for this lens!!!')

############################
# prepare the galaxy catalogues
############################

start_readcat = time.time()

z_ugrizJHK = np.array([])
posx_ugrizJHK = np.array([])
posy_ugrizJHK = np.array([])
imag_ugrizJHK = np.array([])

z_ugriz = np.array([])
posx_ugriz = np.array([])
posy_ugriz = np.array([])
imag_ugriz = np.array([])

root = plane[0:13]

for i in range(4):
    for j in range(4):
        #file_ugrizJHK = '%s%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_%s.images_forNAOJ.txt' % (rootgals,root,i,j,lens)
        file_ugrizJHK = '/Volumes/LaCieDavis/lensing_simulations/SA_galaxies/original/J1206/%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_%s.images_forNAOJ.txt' % (root,i,j,lens)
        #file_ugriz = '/lfs08/rusucs/WFI2033/MSgals/%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt' % (root,i,j)
        #file_ugriz = '/mnt/scratch/rusucs/CFHTLenS/%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt' % (root,i,j)
        file_ugriz = '%s%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_ugriz.images_forNAOJ.txt' % (rootgals,root,i,j)
        if "measured" in type:
            posx__ugrizJHK, posy__ugrizJHK, imag__ugrizJHK, z__ugrizJHK, zspec__ugrizJHK = np.loadtxt(file_ugrizJHK, usecols = (2,3,7,8,1), unpack=True)
            posx__ugriz, posy__ugriz, imag__ugriz, z__ugriz = np.loadtxt(file_ugriz, usecols = (2,3,7,8), unpack=True)
        elif "computed" in type:
            posx__ugrizJHK, posy__ugrizJHK, imag__ugrizJHK, z__ugrizJHK, mstar__ugrizJHK, mhalo__ugrizJHK = np.loadtxt(file_ugrizJHK, usecols = (2,3,6,1,5,4), unpack=True)
            posx__ugriz, posy__ugriz, imag__ugriz, z__ugriz, mstar__ugriz, mhalo__ugriz = np.loadtxt(file_ugriz, usecols = (2,3,6,1,5,4), unpack=True)
        else: sys.exit('Only \"measured\" and \"computed\" are accepted as the fourth argument of the code. Execution stopped.')

        if "measured" in type:
            spec = np.random.uniform(0,1,len(z__ugrizJHK)) # is used to indicate the objects which have "spectroscopic" redshifts
            z__ugrizJHK[(spec < fracspec20) & (imag__ugrizJHK <= 20)] = zspec__ugrizJHK[(spec < fracspec20) & (imag__ugrizJHK <= 20)] # use "spectroscopic" redshifts for the randomy selected objects
            z__ugrizJHK[(spec < fracspec21) & (imag__ugrizJHK > 20) & (imag__ugrizJHK <= 21)] = zspec__ugrizJHK[(spec < fracspec21) & (imag__ugrizJHK > 20) & (imag__ugrizJHK <= 21)]
            z__ugrizJHK[(spec < fracspec22) & (imag__ugrizJHK > 21) & (imag__ugrizJHK <= 22)] = zspec__ugrizJHK[(spec < fracspec22) & (imag__ugrizJHK > 21) & (imag__ugrizJHK <= 22)]
            z__ugrizJHK[(spec < fracspec23) & (imag__ugrizJHK > 22) & (imag__ugrizJHK <= 23)] = zspec__ugrizJHK[(spec < fracspec23) & (imag__ugrizJHK > 22) & (imag__ugrizJHK <= 23)]
            z__ugrizJHK[(spec < fracspec24) & (imag__ugrizJHK > 23) & (imag__ugrizJHK <= 24)] = zspec__ugrizJHK[(spec < fracspec24) & (imag__ugrizJHK > 23) & (imag__ugrizJHK <= 24)]

            # count the galaxies in order to implement the fraction of stars; counting without the redshift gap, which I implement inside the function
            count18 = z__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= 18) & (z__ugrizJHK <= z_s)].size
            count185 = z__ugrizJHK[(imag__ugrizJHK > 18) & (imag__ugrizJHK <= 18.5) & (z__ugrizJHK <= z_s)].size
            count19 = z__ugrizJHK[(imag__ugrizJHK > 18.5) & (imag__ugrizJHK <= 19) & (z__ugrizJHK <= z_s)].size
            count195 = z__ugrizJHK[(imag__ugrizJHK > 19) & (imag__ugrizJHK <= 19.5) & (z__ugrizJHK <= z_s)].size
            count20 = z__ugrizJHK[(imag__ugrizJHK > 19.5) & (imag__ugrizJHK <= 20) & (z__ugrizJHK <= z_s)].size
            count205 = z__ugrizJHK[(imag__ugrizJHK > 20) & (imag__ugrizJHK <= 20.5) & (z__ugrizJHK <= z_s)].size
            count21 = z__ugrizJHK[(imag__ugrizJHK > 20.5) & (imag__ugrizJHK <= 21) & (z__ugrizJHK <= z_s)].size
            count215 = z__ugrizJHK[(imag__ugrizJHK > 21) & (imag__ugrizJHK <= 21.5) & (z__ugrizJHK <= z_s)].size
            count22 = z__ugrizJHK[(imag__ugrizJHK > 21.5) & (imag__ugrizJHK <= 22) & (z__ugrizJHK <= z_s)].size
            count225 = z__ugrizJHK[(imag__ugrizJHK > 22) & (imag__ugrizJHK <= 22.5) & (z__ugrizJHK <= z_s)].size
            count23 = z__ugrizJHK[(imag__ugrizJHK > 22.5) & (imag__ugrizJHK <= 23) & (z__ugrizJHK <= z_s)].size
            count235 = z__ugrizJHK[(imag__ugrizJHK > 23) & (imag__ugrizJHK <= 23.5) & (z__ugrizJHK <= z_s)].size
            count24 = z__ugrizJHK[(imag__ugrizJHK > 23.5) & (imag__ugrizJHK <= 24) & (z__ugrizJHK <= z_s)].size
            #print count18,count185,count19,count195,count20,count205,count21,count215,count22,count225,count23,count235,count24
        
            posxmin = np.min(posx__ugrizJHK) # used to insert random stellar contaminants
            posxmax = np.max(posx__ugrizJHK)
            posymin = np.min(posy__ugrizJHK)
            posymax = np.max(posy__ugrizJHK)

            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_imag_18,cont_z_18 = contaminants(count18,cont_ugrizJHK_18,posxmin,posxmax,posymin,posymax,star_imag_18,star_z_18)
            cont_185,cont_posx_185,cont_posy_185,cont_imag_185,cont_z_185 = contaminants(count185,cont_ugrizJHK_185,posxmin,posxmax,posymin,posymax,star_imag_185,star_z_185)
            cont_19,cont_posx_19,cont_posy_19,cont_imag_19,cont_z_19 = contaminants(count19,cont_ugrizJHK_19,posxmin,posxmax,posymin,posymax,star_imag_19,star_z_19)
            cont_195,cont_posx_195,cont_posy_195,cont_imag_195,cont_z_195 = contaminants(count195,cont_ugrizJHK_195,posxmin,posxmax,posymin,posymax,star_imag_195,star_z_195)
            cont_20,cont_posx_20,cont_posy_20,cont_imag_20,cont_z_20 = contaminants(count20,cont_ugrizJHK_20,posxmin,posxmax,posymin,posymax,star_imag_20,star_z_20)
            cont_205,cont_posx_205,cont_posy_205,cont_imag_205,cont_z_205 = contaminants(count205,cont_ugrizJHK_205,posxmin,posxmax,posymin,posymax,star_imag_205,star_z_205)
            cont_21,cont_posx_21,cont_posy_21,cont_imag_21,cont_z_21 = contaminants(count21,cont_ugrizJHK_21,posxmin,posxmax,posymin,posymax,star_imag_21,star_z_21)
            cont_215,cont_posx_215,cont_posy_215,cont_imag_215,cont_z_215 = contaminants(count215,cont_ugrizJHK_215,posxmin,posxmax,posymin,posymax,star_imag_215,star_z_215)
            cont_22,cont_posx_22,cont_posy_22,cont_imag_22,cont_z_22 = contaminants(count22,cont_ugrizJHK_22,posxmin,posxmax,posymin,posymax,star_imag_22,star_z_22)
            cont_225,cont_posx_225,cont_posy_225,cont_imag_225,cont_z_225 = contaminants(count225,cont_ugrizJHK_225,posxmin,posxmax,posymin,posymax,star_imag_225,star_z_225)
            cont_23,cont_posx_23,cont_posy_23,cont_imag_23,cont_z_23 = contaminants(count23,cont_ugrizJHK_23,posxmin,posxmax,posymin,posymax,star_imag_23,star_z_23)
            cont_235,cont_posx_235,cont_posy_235,cont_imag_235,cont_z_235 = contaminants(count235,cont_ugrizJHK_235,posxmin,posxmax,posymin,posymax,star_imag_235,star_z_235)
            cont_24,cont_posx_24,cont_posy_24,cont_imag_24,cont_z_24 = contaminants(count24,cont_ugrizJHK_24,posxmin,posxmax,posymin,posymax,star_imag_24,star_z_24)
            #print cont_18.size,cont_185.size,cont_19.size,cont_195.size,cont_20.size,cont_205.size,cont_21.size,cont_215.size,cont_22.size,cont_225.size,cont_23.size,cont_235.size,cont_24.size
        
            # masking the fraction of galaxies expected to not be detected as galaxies, due to incompleteness; here also apply the brightmag, limmag, z_s and z gap cuts
            z___ugrizJHK = z__ugrizJHK
            z__ugrizJHK = z___ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s) & ((z___ugrizJHK >= zsup) | (z___ugrizJHK <= zinf)) & (((imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= 18) & (spec < 1 - inc_ugrizJHK_18)) | ((imag__ugrizJHK > 18) & (imag__ugrizJHK <= 18.5) & (spec < 1 - inc_ugrizJHK_185)) | ((imag__ugrizJHK > 18.5) & (imag__ugrizJHK <= 19) & (spec < 1 - inc_ugrizJHK_19)) | ((imag__ugrizJHK > 19) & (imag__ugrizJHK <= 19.5) & (spec < 1 - inc_ugrizJHK_195)) | ((imag__ugrizJHK > 19.5) & (imag__ugrizJHK <= 20) & (spec < 1 - inc_ugrizJHK_20)) | ((imag__ugrizJHK > 20) & (imag__ugrizJHK <= 20.5) & (spec < 1 - inc_ugrizJHK_205)) | ((imag__ugrizJHK > 20.5) & (imag__ugrizJHK <= 21) & (spec < 1 - inc_ugrizJHK_21)) | ((imag__ugrizJHK > 21) & (imag__ugrizJHK <= 21.5) & (spec < 1 - inc_ugrizJHK_215)) | ((imag__ugrizJHK > 21.5) & (imag__ugrizJHK <= 22) & (spec < 1 - inc_ugrizJHK_22)) | ((imag__ugrizJHK > 22) & (imag__ugrizJHK <= 22.5) & (spec < 1 - inc_ugrizJHK_225)) | ((imag__ugrizJHK > 22.5) & (imag__ugrizJHK <= 23) & (spec < 1 - inc_ugrizJHK_23)) | ((imag__ugrizJHK > 23) & (imag__ugrizJHK <= 23.5) & (spec < 1 - inc_ugrizJHK_235)) | ((imag__ugrizJHK > 23.5) & (imag__ugrizJHK <= 24) & (spec < 1 - inc_ugrizJHK_24)))]
            posx__ugrizJHK = posx__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s) & ((z___ugrizJHK >= zsup) | (z___ugrizJHK <= zinf)) & (((imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= 18) & (spec < 1 - inc_ugrizJHK_18)) | ((imag__ugrizJHK > 18) & (imag__ugrizJHK <= 18.5) & (spec < 1 - inc_ugrizJHK_185)) | ((imag__ugrizJHK > 18.5) & (imag__ugrizJHK <= 19) & (spec < 1 - inc_ugrizJHK_19)) | ((imag__ugrizJHK > 19) & (imag__ugrizJHK <= 19.5) & (spec < 1 - inc_ugrizJHK_195)) | ((imag__ugrizJHK > 19.5) & (imag__ugrizJHK <= 20) & (spec < 1 - inc_ugrizJHK_20)) | ((imag__ugrizJHK > 20) & (imag__ugrizJHK <= 20.5) & (spec < 1 - inc_ugrizJHK_205)) | ((imag__ugrizJHK > 20.5) & (imag__ugrizJHK <= 21) & (spec < 1 - inc_ugrizJHK_21)) | ((imag__ugrizJHK > 21) & (imag__ugrizJHK <= 21.5) & (spec < 1 - inc_ugrizJHK_215)) | ((imag__ugrizJHK > 21.5) & (imag__ugrizJHK <= 22) & (spec < 1 - inc_ugrizJHK_22)) | ((imag__ugrizJHK > 22) & (imag__ugrizJHK <= 22.5) & (spec < 1 - inc_ugrizJHK_225)) | ((imag__ugrizJHK > 22.5) & (imag__ugrizJHK <= 23) & (spec < 1 - inc_ugrizJHK_23)) | ((imag__ugrizJHK > 23) & (imag__ugrizJHK <= 23.5) & (spec < 1 - inc_ugrizJHK_235)) | ((imag__ugrizJHK > 23.5) & (imag__ugrizJHK <= 24) & (spec < 1 - inc_ugrizJHK_24)))]
            posy__ugrizJHK = posy__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s) & ((z___ugrizJHK >= zsup) | (z___ugrizJHK <= zinf)) & (((imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= 18) & (spec < 1 - inc_ugrizJHK_18)) | ((imag__ugrizJHK > 18) & (imag__ugrizJHK <= 18.5) & (spec < 1 - inc_ugrizJHK_185)) | ((imag__ugrizJHK > 18.5) & (imag__ugrizJHK <= 19) & (spec < 1 - inc_ugrizJHK_19)) | ((imag__ugrizJHK > 19) & (imag__ugrizJHK <= 19.5) & (spec < 1 - inc_ugrizJHK_195)) | ((imag__ugrizJHK > 19.5) & (imag__ugrizJHK <= 20) & (spec < 1 - inc_ugrizJHK_20)) | ((imag__ugrizJHK > 20) & (imag__ugrizJHK <= 20.5) & (spec < 1 - inc_ugrizJHK_205)) | ((imag__ugrizJHK > 20.5) & (imag__ugrizJHK <= 21) & (spec < 1 - inc_ugrizJHK_21)) | ((imag__ugrizJHK > 21) & (imag__ugrizJHK <= 21.5) & (spec < 1 - inc_ugrizJHK_215)) | ((imag__ugrizJHK > 21.5) & (imag__ugrizJHK <= 22) & (spec < 1 - inc_ugrizJHK_22)) | ((imag__ugrizJHK > 22) & (imag__ugrizJHK <= 22.5) & (spec < 1 - inc_ugrizJHK_225)) | ((imag__ugrizJHK > 22.5) & (imag__ugrizJHK <= 23) & (spec < 1 - inc_ugrizJHK_23)) | ((imag__ugrizJHK > 23) & (imag__ugrizJHK <= 23.5) & (spec < 1 - inc_ugrizJHK_235)) | ((imag__ugrizJHK > 23.5) & (imag__ugrizJHK <= 24) & (spec < 1 - inc_ugrizJHK_24)))]
            imag__ugrizJHK = imag__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s) & ((z___ugrizJHK >= zsup) | (z___ugrizJHK <= zinf)) & (((imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= 18) & (spec < 1 - inc_ugrizJHK_18)) | ((imag__ugrizJHK > 18) & (imag__ugrizJHK <= 18.5) & (spec < 1 - inc_ugrizJHK_185)) | ((imag__ugrizJHK > 18.5) & (imag__ugrizJHK <= 19) & (spec < 1 - inc_ugrizJHK_19)) | ((imag__ugrizJHK > 19) & (imag__ugrizJHK <= 19.5) & (spec < 1 - inc_ugrizJHK_195)) | ((imag__ugrizJHK > 19.5) & (imag__ugrizJHK <= 20) & (spec < 1 - inc_ugrizJHK_20)) | ((imag__ugrizJHK > 20) & (imag__ugrizJHK <= 20.5) & (spec < 1 - inc_ugrizJHK_205)) | ((imag__ugrizJHK > 20.5) & (imag__ugrizJHK <= 21) & (spec < 1 - inc_ugrizJHK_21)) | ((imag__ugrizJHK > 21) & (imag__ugrizJHK <= 21.5) & (spec < 1 - inc_ugrizJHK_215)) | ((imag__ugrizJHK > 21.5) & (imag__ugrizJHK <= 22) & (spec < 1 - inc_ugrizJHK_22)) | ((imag__ugrizJHK > 22) & (imag__ugrizJHK <= 22.5) & (spec < 1 - inc_ugrizJHK_225)) | ((imag__ugrizJHK > 22.5) & (imag__ugrizJHK <= 23) & (spec < 1 - inc_ugrizJHK_23)) | ((imag__ugrizJHK > 23) & (imag__ugrizJHK <= 23.5) & (spec < 1 - inc_ugrizJHK_235)) | ((imag__ugrizJHK > 23.5) & (imag__ugrizJHK <= 24) & (spec < 1 - inc_ugrizJHK_24)))]

            # inserting the stellar contaminants
            if limmag == 23:
                z_ugrizJHK = np.concatenate((z_ugrizJHK,z__ugrizJHK,cont_z_18,cont_z_185,cont_z_19,cont_z_195,cont_z_20,cont_z_205,cont_z_21,cont_z_215,cont_z_22,cont_z_225,cont_z_23))
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23))
                imag_ugrizJHK = np.concatenate((imag_ugrizJHK,imag__ugrizJHK,cont_imag_18,cont_imag_185,cont_imag_19,cont_imag_195,cont_imag_20,cont_imag_205,cont_imag_21,cont_imag_215,cont_imag_22,cont_imag_225,cont_imag_23))
            if limmag == 24:
                z_ugrizJHK = np.concatenate((z_ugrizJHK,z__ugrizJHK,cont_z_18,cont_z_185,cont_z_19,cont_z_195,cont_z_20,cont_z_205,cont_z_21,cont_z_215,cont_z_22,cont_z_225,cont_z_23,cont_z_235,cont_z_24))
                posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235,cont_posx_24))
                posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235,cont_posy_24))
                imag_ugrizJHK = np.concatenate((imag_ugrizJHK,imag__ugrizJHK,cont_imag_18,cont_imag_185,cont_imag_19,cont_imag_195,cont_imag_20,cont_imag_205,cont_imag_21,cont_imag_215,cont_imag_22,cont_imag_225,cont_imag_23,cont_imag_235,cont_imag_24))
            #if len(imag_ugrizJHK) > 0: print np.max(imag_ugrizJHK)
            
            # repeat for ugriz catalogues

            count18 = z__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= 18) & (z__ugriz <= z_s)].size
            count185 = z__ugriz[(imag__ugriz > 18) & (imag__ugriz <= 18.5) & (z__ugriz <= z_s)].size
            count19 = z__ugriz[(imag__ugriz > 18.5) & (imag__ugriz <= 19) & (z__ugriz <= z_s)].size
            count195 = z__ugriz[(imag__ugriz > 19) & (imag__ugriz <= 19.5) & (z__ugriz <= z_s)].size
            count20 = z__ugriz[(imag__ugriz > 19.5) & (imag__ugriz <= 20) & (z__ugriz <= z_s)].size
            count205 = z__ugriz[(imag__ugriz > 20) & (imag__ugriz <= 20.5) & (z__ugriz <= z_s)].size
            count21 = z__ugriz[(imag__ugriz > 20.5) & (imag__ugriz <= 21) & (z__ugriz <= z_s)].size
            count215 = z__ugriz[(imag__ugriz > 21) & (imag__ugriz <= 21.5) & (z__ugriz <= z_s)].size
            count22 = z__ugriz[(imag__ugriz > 21.5) & (imag__ugriz <= 22) & (z__ugriz <= z_s)].size
            count225 = z__ugriz[(imag__ugriz > 22) & (imag__ugriz <= 22.5) & (z__ugriz <= z_s)].size
            count23 = z__ugriz[(imag__ugriz > 22.5) & (imag__ugriz <= 23) & (z__ugriz <= z_s)].size
            count235 = z__ugriz[(imag__ugriz > 23) & (imag__ugriz <= 23.5) & (z__ugriz <= z_s)].size
            count24 = z__ugriz[(imag__ugriz > 23.5) & (imag__ugriz <= 24) & (z__ugriz <= z_s)].size
        
            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_imag_18,cont_z_18 = contaminants(count18,cont_ugriz_18,posxmin,posxmax,posymin,posymax,star_imag_18,star_z_18)
            cont_185,cont_posx_185,cont_posy_185,cont_imag_185,cont_z_185 = contaminants(count185,cont_ugriz_185,posxmin,posxmax,posymin,posymax,star_imag_185,star_z_185)
            cont_19,cont_posx_19,cont_posy_19,cont_imag_19,cont_z_19 = contaminants(count19,cont_ugriz_19,posxmin,posxmax,posymin,posymax,star_imag_19,star_z_19)
            cont_195,cont_posx_195,cont_posy_195,cont_imag_195,cont_z_195 = contaminants(count195,cont_ugriz_195,posxmin,posxmax,posymin,posymax,star_imag_195,star_z_195)
            cont_20,cont_posx_20,cont_posy_20,cont_imag_20,cont_z_20 = contaminants(count20,cont_ugriz_20,posxmin,posxmax,posymin,posymax,star_imag_20,star_z_20)
            cont_205,cont_posx_205,cont_posy_205,cont_imag_205,cont_z_205 = contaminants(count205,cont_ugriz_205,posxmin,posxmax,posymin,posymax,star_imag_205,star_z_205)
            cont_21,cont_posx_21,cont_posy_21,cont_imag_21,cont_z_21 = contaminants(count21,cont_ugriz_21,posxmin,posxmax,posymin,posymax,star_imag_21,star_z_21)
            cont_215,cont_posx_215,cont_posy_215,cont_imag_215,cont_z_215 = contaminants(count215,cont_ugriz_215,posxmin,posxmax,posymin,posymax,star_imag_215,star_z_215)
            cont_22,cont_posx_22,cont_posy_22,cont_imag_22,cont_z_22 = contaminants(count22,cont_ugriz_22,posxmin,posxmax,posymin,posymax,star_imag_22,star_z_22)
            cont_225,cont_posx_225,cont_posy_225,cont_imag_225,cont_z_225 = contaminants(count225,cont_ugriz_225,posxmin,posxmax,posymin,posymax,star_imag_225,star_z_225)
            cont_23,cont_posx_23,cont_posy_23,cont_imag_23,cont_z_23 = contaminants(count23,cont_ugriz_23,posxmin,posxmax,posymin,posymax,star_imag_23,star_z_23)
            cont_235,cont_posx_235,cont_posy_235,cont_imag_235,cont_z_235 = contaminants(count235,cont_ugriz_235,posxmin,posxmax,posymin,posymax,star_imag_235,star_z_235)
            cont_24,cont_posx_24,cont_posy_24,cont_imag_24,cont_z_24 = contaminants(count24,cont_ugriz_24,posxmin,posxmax,posymin,posymax,star_imag_24,star_z_24)

            spec = np.random.uniform(0,1,len(z__ugriz))
            z___ugriz = z__ugriz
            z__ugriz = z___ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf)) & (((imag__ugriz > brightmag) & (imag__ugriz <= 18) & (spec < 1 - inc_ugriz_18)) | ((imag__ugriz > 18) & (imag__ugriz <= 18.5) & (spec < 1 - inc_ugriz_185)) | ((imag__ugriz > 18.5) & (imag__ugriz <= 19) & (spec < 1 - inc_ugriz_19)) | ((imag__ugriz > 19) & (imag__ugriz <= 19.5) & (spec < 1 - inc_ugriz_195)) | ((imag__ugriz > 19.5) & (imag__ugriz <= 20) & (spec < 1 - inc_ugriz_20)) | ((imag__ugriz > 20) & (imag__ugriz <= 20.5) & (spec < 1 - inc_ugriz_205)) | ((imag__ugriz > 20.5) & (imag__ugriz <= 21) & (spec < 1 - inc_ugriz_21)) | ((imag__ugriz > 21) & (imag__ugriz <= 21.5) & (spec < 1 - inc_ugriz_215)) | ((imag__ugriz > 21.5) & (imag__ugriz <= 22) & (spec < 1 - inc_ugriz_22)) | ((imag__ugriz > 22) & (imag__ugriz <= 22.5) & (spec < 1 - inc_ugriz_225)) | ((imag__ugriz > 22.5) & (imag__ugriz <= 23) & (spec < 1 - inc_ugriz_23)) | ((imag__ugriz > 23) & (imag__ugriz <= 23.5) & (spec < 1 - inc_ugriz_235)) | ((imag__ugriz > 23.5) & (imag__ugriz <= 24) & (spec < 1 - inc_ugriz_24)))]
            posx__ugriz = posx__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf)) & (((imag__ugriz > brightmag) & (imag__ugriz <= 18) & (spec < 1 - inc_ugriz_18)) | ((imag__ugriz > 18) & (imag__ugriz <= 18.5) & (spec < 1 - inc_ugriz_185)) | ((imag__ugriz > 18.5) & (imag__ugriz <= 19) & (spec < 1 - inc_ugriz_19)) | ((imag__ugriz > 19) & (imag__ugriz <= 19.5) & (spec < 1 - inc_ugriz_195)) | ((imag__ugriz > 19.5) & (imag__ugriz <= 20) & (spec < 1 - inc_ugriz_20)) | ((imag__ugriz > 20) & (imag__ugriz <= 20.5) & (spec < 1 - inc_ugriz_205)) | ((imag__ugriz > 20.5) & (imag__ugriz <= 21) & (spec < 1 - inc_ugriz_21)) | ((imag__ugriz > 21) & (imag__ugriz <= 21.5) & (spec < 1 - inc_ugriz_215)) | ((imag__ugriz > 21.5) & (imag__ugriz <= 22) & (spec < 1 - inc_ugriz_22)) | ((imag__ugriz > 22) & (imag__ugriz <= 22.5) & (spec < 1 - inc_ugriz_225)) | ((imag__ugriz > 22.5) & (imag__ugriz <= 23) & (spec < 1 - inc_ugriz_23)) | ((imag__ugriz > 23) & (imag__ugriz <= 23.5) & (spec < 1 - inc_ugriz_235)) | ((imag__ugriz > 23.5) & (imag__ugriz <= 24) & (spec < 1 - inc_ugriz_24)))]
            posy__ugriz = posy__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf)) & (((imag__ugriz > brightmag) & (imag__ugriz <= 18) & (spec < 1 - inc_ugriz_18)) | ((imag__ugriz > 18) & (imag__ugriz <= 18.5) & (spec < 1 - inc_ugriz_185)) | ((imag__ugriz > 18.5) & (imag__ugriz <= 19) & (spec < 1 - inc_ugriz_19)) | ((imag__ugriz > 19) & (imag__ugriz <= 19.5) & (spec < 1 - inc_ugriz_195)) | ((imag__ugriz > 19.5) & (imag__ugriz <= 20) & (spec < 1 - inc_ugriz_20)) | ((imag__ugriz > 20) & (imag__ugriz <= 20.5) & (spec < 1 - inc_ugriz_205)) | ((imag__ugriz > 20.5) & (imag__ugriz <= 21) & (spec < 1 - inc_ugriz_21)) | ((imag__ugriz > 21) & (imag__ugriz <= 21.5) & (spec < 1 - inc_ugriz_215)) | ((imag__ugriz > 21.5) & (imag__ugriz <= 22) & (spec < 1 - inc_ugriz_22)) | ((imag__ugriz > 22) & (imag__ugriz <= 22.5) & (spec < 1 - inc_ugriz_225)) | ((imag__ugriz > 22.5) & (imag__ugriz <= 23) & (spec < 1 - inc_ugriz_23)) | ((imag__ugriz > 23) & (imag__ugriz <= 23.5) & (spec < 1 - inc_ugriz_235)) | ((imag__ugriz > 23.5) & (imag__ugriz <= 24) & (spec < 1 - inc_ugriz_24)))]
            imag__ugriz = imag__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf)) & (((imag__ugriz > brightmag) & (imag__ugriz <= 18) & (spec < 1 - inc_ugriz_18)) | ((imag__ugriz > 18) & (imag__ugriz <= 18.5) & (spec < 1 - inc_ugriz_185)) | ((imag__ugriz > 18.5) & (imag__ugriz <= 19) & (spec < 1 - inc_ugriz_19)) | ((imag__ugriz > 19) & (imag__ugriz <= 19.5) & (spec < 1 - inc_ugriz_195)) | ((imag__ugriz > 19.5) & (imag__ugriz <= 20) & (spec < 1 - inc_ugriz_20)) | ((imag__ugriz > 20) & (imag__ugriz <= 20.5) & (spec < 1 - inc_ugriz_205)) | ((imag__ugriz > 20.5) & (imag__ugriz <= 21) & (spec < 1 - inc_ugriz_21)) | ((imag__ugriz > 21) & (imag__ugriz <= 21.5) & (spec < 1 - inc_ugriz_215)) | ((imag__ugriz > 21.5) & (imag__ugriz <= 22) & (spec < 1 - inc_ugriz_22)) | ((imag__ugriz > 22) & (imag__ugriz <= 22.5) & (spec < 1 - inc_ugriz_225)) | ((imag__ugriz > 22.5) & (imag__ugriz <= 23) & (spec < 1 - inc_ugriz_23)) | ((imag__ugriz > 23) & (imag__ugriz <= 23.5) & (spec < 1 - inc_ugriz_235)) | ((imag__ugriz > 23.5) & (imag__ugriz <= 24) & (spec < 1 - inc_ugriz_24)))]
        
            if limmag == 23:
                z_ugriz = np.concatenate((z_ugriz,z__ugriz,cont_z_18,cont_z_185,cont_z_19,cont_z_195,cont_z_20,cont_z_205,cont_z_21,cont_z_215,cont_z_22,cont_z_225,cont_z_23))
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23))
                imag_ugriz = np.concatenate((imag_ugriz,imag__ugriz,cont_imag_18,cont_imag_185,cont_imag_19,cont_imag_195,cont_imag_20,cont_imag_205,cont_imag_21,cont_imag_215,cont_imag_22,cont_imag_225,cont_imag_23))
            if limmag == 24:
                z_ugriz = np.concatenate((z_ugriz,z__ugriz,cont_z_18,cont_z_185,cont_z_19,cont_z_195,cont_z_20,cont_z_205,cont_z_21,cont_z_215,cont_z_22,cont_z_225,cont_z_23,cont_z_235,cont_z_24))
                posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23,cont_posx_235,cont_posx_24))
                posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23,cont_posy_235,cont_posy_24))
                imag_ugriz = np.concatenate((imag_ugriz,imag__ugriz,cont_imag_18,cont_imag_185,cont_imag_19,cont_imag_195,cont_imag_20,cont_imag_205,cont_imag_21,cont_imag_215,cont_imag_22,cont_imag_225,cont_imag_23,cont_imag_235,cont_imag_24))
            #if len(imag_ugriz) > 0: print np.max(imag_ugriz)

        else:
            z___ugrizJHK = z__ugrizJHK
            z__ugrizJHK = z___ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s) & ((z___ugrizJHK >= zsup) | (z___ugrizJHK <= zinf))]
            posx__ugrizJHK = posx__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s)]
            posy__ugrizJHK = posy__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s)]
            imag__ugrizJHK = imag__ugrizJHK[(imag__ugrizJHK > brightmag) & (imag__ugrizJHK <= limmag) & (z___ugrizJHK <= z_s)]
            z_ugrizJHK = np.concatenate((z_ugrizJHK,z__ugrizJHK))
            posx_ugrizJHK = np.concatenate((posx_ugrizJHK,posx__ugrizJHK))
            posy_ugrizJHK = np.concatenate((posy_ugrizJHK,posy__ugrizJHK))
            imag_ugrizJHK = np.concatenate((imag_ugrizJHK,imag__ugrizJHK))

            z___ugriz = z__ugriz
            z__ugriz = z___ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf))]
            posx__ugriz = posx__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf))]
            posy__ugriz = posy__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf))]
            imag__ugriz = imag__ugriz[(imag__ugriz > brightmag) & (imag__ugriz <= limmag) & (z___ugriz <= z_s) & ((z___ugriz >= zsup) | (z___ugriz <= zinf))]
            z_ugriz = np.concatenate((z_ugriz,z__ugriz))
            posx_ugriz = np.concatenate((posx_ugriz,posx__ugriz))
            posy_ugriz = np.concatenate((posy_ugriz,posy__ugriz))
            imag_ugriz = np.concatenate((imag_ugriz,imag__ugriz))

cat_ugrizJHK = np.c_[z_ugrizJHK,posx_ugrizJHK,posy_ugrizJHK,imag_ugrizJHK]
del z___ugrizJHK
del z_ugrizJHK
del posx_ugrizJHK
del posy_ugrizJHK
del imag_ugrizJHK
del z__ugrizJHK
del posx__ugrizJHK
del posy__ugrizJHK
del imag__ugrizJHK
cat_ugriz = np.c_[z_ugriz,posx_ugriz,posy_ugriz,imag_ugriz]
del z___ugriz
del z_ugriz
del posx_ugriz
del posy_ugriz
del imag_ugriz
del z__ugriz
del posx__ugriz
del posy__ugriz
del imag__ugriz

index_z = 0
index_posx = 1
index_posy = 2
index_mstar = 3
index_imag = 4
index_Mhalo = 5
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
bands = "griK"
weightedcounts(cat,spacing,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,bands)

cat = cat_ugriz
bands = "ugriz"
weightedcounts(cat,spacing,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,bands)

print "Computed weights in ", time.time() - start_weights, "seconds"

print(" Field done in --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
