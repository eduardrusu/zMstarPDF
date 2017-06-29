# This code uses the Millenium Sumilation convergence and shear maps as well as the associated SA catalogue of  galaxies, in order to compute the weighted counts for fields centered around each kappa and gamma point. This is done for a variety of limiting magnitudes, aperture radii, and weights. This version uses the gamma & kappa from MS plane 35, and only considers i<23 galaxies (designed for WFI2033).

# run from the working directory (with all the files) for each of the 64 simulation fields, with the following arguments: lens name, field name, outer mask radius, type, inner mask radius; e.g.: python kappamed_insertstarsnobeta_griK_i23.py WFI2033 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_35_f 120 measured 5

import numpy as np
import scipy
import sys
import os
from os import system
from scipy import special
import pandas as pd
import time
import distances
from scipy.interpolate import griddata

############################
# function definitions
############################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...") # use for debugging

def readbinary(replacestr,plane,radiusstr):
    replace = plane + replacestr
    os.system("sed \"11s/.*/  const char kappa_file_name[]   = \\\"\%s\\\";/\" readKappaBinary.c > readKappaBinary_%smed%s.c_" % (rootkappagamma + replace,plane,radiusstr))
    os.system("sed \"35s/.*/  fpt = fopen (\\\"kappa_values_%smed%s.dat\\\", \\\"w\\\");/\"  readKappaBinary_%smed%s.c_ > readKappaBinary_%smed%s.c" % (plane,radiusstr,plane,radiusstr,plane,radiusstr))
    os.system("rm readKappaBinary_%smed%s.c_" % (plane,radiusstr))
    os.system("gcc readKappaBinary_%smed%s.c -o compiled_%smed%s.out" % (plane,radiusstr,plane,radiusstr))
    os.system("./compiled_%smed%s.out" % (plane,radiusstr))
    os.system("rm readKappaBinary_%smed%s.c" % (plane,radiusstr))
    os.system("rm compiled_%smed%s.out" % (plane,radiusstr))

def contaminants(count,cont_ugr,posxmin,posxmax,posymin,posymax,star_imag,star_z,star_mstar):
    cont = np.random.random_integers(0,499,int(count * cont_ugr)) # randomly select from the star catalogues which contain 500 stars each
    cont_posx = np.random.uniform(posxmin,posxmax,int(count * cont_ugr))
    cont_posy = np.random.uniform(posymin,posymax,int(count * cont_ugr))
    cont_imag = star_imag[cont]
    cont_z = star_z[cont]
    cont_mstar = star_mstar[cont]
    return cont, cont_posx, cont_posy, cont_imag, cont_z, cont_mstar

def weightedcounts(cat,spacing,radius,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,type,lens,plane,bands,inner):
    initialized = 0
    for i in range(spacing):
        for j in range(spacing):
            print "(",i+1,",",j+1,") / (",spacing,",",spacing, ") for radius", radius
            grid_x, grid_y = np.mgrid[lim1D + i:4096 - lim1D - (4096 - 2 * lim1D) % spacing - spacing + i:complex(0,cells_on_a_side), lim1D + j:4096 - lim1D - (4096 - 2 * lim1D) % spacing - spacing + j:complex(0,cells_on_a_side)] # the grid containing the kappa pixel at the center of cells
            cellx = grid_x.flatten()
            celly = grid_y.flatten()
            cellxy = np.array([cellx,celly])
            posxy = np.array([-0.5 * L_field  + (1 + cellxy[0] + 0.5) * L_pix, -0.5 * L_field  + (1 + cellxy[1] + 0.5) * L_pix])
            cellkappagamma = np.c_[cells,kappagamma[:,][(cellxy[0] * 4096 + cellxy[1]).astype(int)]]
            index = griddata(posxy.T, cells, (cat[:,index_posx], cat[:,index_posy]), method='nearest')
            sep = np.sqrt((posxy[0][index.astype(int)]-cat[:,index_posx])**2 + (posxy[1][index.astype(int)]-cat[:,index_posy])**2)*1/degree*3600
            cat_msk = np.c_[cat,index,sep]
            cat_msk = cat_msk[cat_msk[:,index_sep] <= radius] # mask objects at distance larger than the aperture from the center
            cat_msk = cat_msk[cat_msk[:,index_sep] > innermsk] # inner mask
            cat_msk[:,index_sep][cat_msk[:,index_sep] < 10] = 10
            #cat = np.c_[z,posx,posy,mstar,imag...]
            
            #cat_msk = cat_msk[cat_msk[:,index_imag] <= 23]
            p_zweight = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'zweight':1.0 * (z_s * cat_msk[:,index_z]) - cat_msk[:,index_z]**2})
            p_zoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'zoverr':1.0 * ((z_s * cat_msk[:,index_z]) - cat_msk[:,index_z]**2) / cat_msk[:,index_sep]})
            p_oneoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'oneoverr':1.0 / cat_msk[:,index_sep]})
            p_mass = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'mass':1.0 * cat_msk[:,index_mstar]})
            p_mass2 = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'mass2':1.0 * (cat_msk[:,index_mstar]) ** 2})
            p_mass3 = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'mass3':1.0 * (cat_msk[:,index_mstar]) ** 3})
            p_massoverr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'massoverr':1.0 * cat_msk[:,index_mstar] / cat_msk[:,index_sep]})
            p_mass2overr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'mass2overr':1.0 * (cat_msk[:,index_mstar] ** 2) / cat_msk[:,index_sep]})
            p_mass3overr = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'mass3overr':1.0 * (cat_msk[:,index_mstar] ** 3) / cat_msk[:,index_sep]})
            p_flexion = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'flexion':cat_msk[:,index_mstar]  / (cat_msk[:,index_sep] ** 3)})
            p_tidal = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'tidal':cat_msk[:,index_mstar] / (cat_msk[:,index_sep] ** 2)})
            p_SIS = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'SIS':np.sqrt(cat_msk[:,index_mstar]) / cat_msk[:,index_sep]})
            p_SIShalo = pd.DataFrame({'cell':cat_msk[:,index_index].astype(int),'SIShalo':np.sqrt(cat_msk[:,index_Mhalo]) / cat_msk[:,index_sep]})
            
            w_gal_23 = np.bincount(cat_msk[:,index_index].astype(int))
            if (p_zweight.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_zoverr.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_oneoverr.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_mass.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_mass2.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_mass3.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_massoverr.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_mass2overr.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_mass3overr.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_flexion.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_tidal.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_SIS.groupby(['cell']).median().values[:,0].size == w_gal_23.size) and (p_SIShalo.groupby(['cell']).median().values[:,0].size == w_gal_23.size): # it seems for 45 arcsec aperture only, in some cases the two vary by 1 (e.g. 24335 and 24336)
                w_gal_23 = np.bincount(cat_msk[:,index_index].astype(int))
                w_zweight_23 = p_zweight.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass_23 = p_mass.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass2_23 = p_mass2.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass3_23 = p_mass3.groupby(['cell']).median().values[:,0] * w_gal_23
                w_oneoverr_23 = p_oneoverr.groupby(['cell']).median().values[:,0] * w_gal_23
                w_zoverr_23 = p_zoverr.groupby(['cell']).median().values[:,0] * w_gal_23
                w_massoverr_23 = p_massoverr.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass2overr_23 = p_mass2overr.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass3overr_23 = p_mass3overr.groupby(['cell']).median().values[:,0] * w_gal_23
                w_flexion_23 = p_flexion.groupby(['cell']).median().values[:,0] * w_gal_23
                w_tidal_23 = p_tidal.groupby(['cell']).median().values[:,0] * w_gal_23
                w_SIS_23 = p_SIS.groupby(['cell']).median().values[:,0] * w_gal_23
                w_SIShalo_23 = p_SIShalo.groupby(['cell']).median().values[:,0] * w_gal_23
                w_mass2rms_23 = np.sqrt(w_mass2_23)
                w_mass3rms_23 = scipy.special.cbrt(w_mass3_23)
                w_mass2overrms_23 = np.sqrt(w_mass2overr_23)
                w_mass3overrms_23 = scipy.special.cbrt(w_mass3overr_23)
            
                cellkappagamma = np.c_[cellkappagamma,w_gal_23,w_zweight_23,w_mass_23,w_mass2_23,w_mass3_23,w_oneoverr_23,w_zoverr_23,w_massoverr_23,w_mass2overr_23,w_mass3overr_23,w_mass2rms_23,w_mass3rms_23,w_mass2overrms_23,w_mass3overrms_23,w_flexion_23,w_tidal_23,w_SIS_23,w_SIShalo_23]
                if initialized == 0:
                    os.system('rm nobeta%s%smedinject_%s_%s_%s_%s_%sarcsecinnermsk.cat' % (pln,type,bands,lens,plane[0:13],radius,inner))
                    f=open('nobeta%s%smedinject_%s_%s_%s_%s_%sarcsecinnermsk.cat' % (pln,type,bands,lens,plane[0:13],radius,inner),'ab') # open in binary
                    output = "# ID kappa gamma1 gamma2 w_gal_23 w_zweight_23 w_mass_23 w_mass2_23 w_mass3_23 w_oneoverr_23 w_zoverr_23 w_massoverr_23 w_mass2overr_23 w_mass3overr_23 w_mass2rms_23 w_mass3rms_23 w_mass2overrms_23 w_mass3overrms_23 w_flexion_23 w_tidal_23 w_SIS_23 w_SIShalo_23 \n"
                    f.write(output) # needs to be done inside the loop because otherwise it crashes
                    initialized = 1
                output = ""
                for k in range(cellkappagamma.shape[0]):
                    output = output + str('%.1d' % cellkappagamma[k][1])+ "\t" + str('%.5e' % cellkappagamma[k][2])+ "\t" + str('%.5e' % cellkappagamma[k][3])+ "\t" + str('%.5e' % cellkappagamma[k][4])+ "\t" + str('%.1d' % cellkappagamma[k][5])+ "\t" + str('%.5e' % cellkappagamma[k][6])+ "\t" + str('%.4e' % cellkappagamma[k][7])+ "\t" + str('%.4e' % cellkappagamma[k][8])+ "\t" + str('%.4e' % cellkappagamma[k][9])+ "\t" + str('%.4e' % cellkappagamma[k][10])+ "\t" + str('%.4e' % cellkappagamma[k][11])+ "\t" + str('%.4e' % cellkappagamma[k][12])+ "\t" + str('%.4e' % cellkappagamma[k][13])+ "\t" + str('%.4e' % cellkappagamma[k][14])+ "\t" + str('%.4e' % cellkappagamma[k][15])+ "\t" + str('%.4e' % cellkappagamma[k][16])+ "\t" + str('%.4e' % cellkappagamma[k][17])+ "\t" + str('%.4e' % cellkappagamma[k][18])+ "\t" + str('%.4e' % cellkappagamma[k][19])+ "\t" + str('%.4e' % cellkappagamma[k][20])+ "\t" + str('%.4e' % cellkappagamma[k][21])+ "\t" + str('%.4e' % cellkappagamma[k][22])+ "\n"
                f.write(output)
    f.close()

############################
# lens information
############################

start_time = time.time()

rootkappagamma = "\/mfst01a\/rusucs\/results_kappagamma\/mstarsim\/" # this is where I orininally saved the files
#rootwghtratios = "/wam03b/rusucs/MSwghtratios/" # save to the partition with more available space
rootgals = "/mfst01a/rusucs/WFI2033/MSgals/"
#rootkappagamma = "/mfst01a/rusucs/results_kappagamma/mstarsim/"
lens = str(sys.argv[1])
plane = str(sys.argv[2])
radiusstr = str(sys.argv[3])
type = str(sys.argv[4]) # computed or measured
innermsk = int(str(sys.argv[5])) # inner mask in arcsec

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
    if (radiusstr == "120"):
        radius = 120
        fracspec20 = 1
        fracspec21 = 0.87
        fracspec22 = 0.44
        fracspec23 = 0.08
    if (radiusstr == "60"):
        radius = 60
        fracspec20 = 1 # gal+stars
        fracspec21 = 0.8
        fracspec22 = 0.78
        fracspec23 = 0.15
    if (radiusstr == "90"):
        radius = 90
        fracspec20 = 1
        fracspec21 = 0.85
        fracspec22 = 0.5
        fracspec23 = 0.1
if lens == "WFI2033":
    z_s = 1.66
    z_l = 0.66
    brightmag = 16.90
    #innermsk = 8 # inner mask in arcsec
    pln = 35
    if (radiusstr == "45"):
        hstcoverage = 1
        radius = 45
        fracspec20 = 1 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.73
        fracspec23 = 0.15
    if (radiusstr == "120"):
        hstcoverage = 1 * 0.47
        radius = 120
        fracspec20 = 0.69
        fracspec21 = 0.81
        fracspec22 = 0.52
        fracspec23 = 0.07
    if (radiusstr == "60"):
        hstcoverage = 1
        radius = 60
        fracspec20 = 1 # gal+stars
        fracspec21 = 1
        fracspec22 = 0.70
        fracspec23 = 0.11
    if (radiusstr == "90"):
        hstcoverage = 1 * 0.47 * (90.0**2)/(120.0**2)
        radius = 90
        fracspec20 = 0.71
        fracspec21 = 0.88
        fracspec22 = 0.6
        fracspec23 = 0.06
if lens == "HE1104":
    z_s = 2.32
#pln = 30 & 31
if lens == "RX1131":
    z_s = 0.66
#pln = 45 & 46
if lens == "J1206":
    z_s = 1.79
    pln = 34

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
nospec_235 = 1
nospec_24 = 1

cont_griK_18 = cont_h12_18 * nospec_18 * (1 - hstcoverage)
cont_griK_185 = cont_h12_185 * nospec_185 * (1 - hstcoverage)
cont_griK_19 = cont_h12_19 * nospec_19 * (1 - hstcoverage)
cont_griK_195 = cont_h12_195 * nospec_195 * (1 - hstcoverage)
cont_griK_20 = cont_h12_20 * nospec_20 * (1 - hstcoverage)
cont_griK_205 = cont_h12_205 * nospec_205 * (1 - hstcoverage)
cont_griK_21 = cont_h12_21 * nospec_21 * (1 - hstcoverage)
cont_griK_215 = cont_h12_215 * nospec_215 * (1 - hstcoverage)
cont_griK_22 = cont_h12_22 * nospec_22 * (1 - hstcoverage)
cont_griK_225 = cont_h12_225 * nospec_225 * (1 - hstcoverage)
cont_griK_23 = cont_h12_23 * nospec_23 * (1 - hstcoverage)
cont_griK_235 = cont_h12_235 * nospec_235
cont_griK_24 = cont_h12_24 * nospec_24

inc_griK_18 = inc_h12_18 * nospec_18 * (1 - hstcoverage)
inc_griK_185 = inc_h12_185 * nospec_185 * (1 - hstcoverage)
inc_griK_19 = inc_h12_19 * nospec_19 * (1 - hstcoverage)
inc_griK_195 = inc_h12_195 * nospec_195 * (1 - hstcoverage)
inc_griK_20 = inc_h12_20 * nospec_20 * (1 - hstcoverage)
inc_griK_205 = inc_h12_205 * nospec_205 * (1 - hstcoverage)
inc_griK_21 = inc_h12_21 * nospec_21 * (1 - hstcoverage)
inc_griK_215 = inc_h12_215 * nospec_215 * (1 - hstcoverage)
inc_griK_22 = inc_h12_22 * nospec_22 * (1 - hstcoverage)
inc_griK_225 = inc_h12_225 * nospec_225 * (1 - hstcoverage)
inc_griK_23 = inc_h12_23 * nospec_23 * (1 - hstcoverage)
inc_griK_235 = inc_h12_235 * nospec_235
inc_griK_24 = inc_h12_24 * nospec_24

# read the stellar contaminants I wil insert

star_imag_18, star_z_18, star_mstar_18 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star018zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_185, star_z_185, star_mstar_185 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star18185zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_19, star_z_19, star_mstar_19 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star18519zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_195, star_z_195, star_mstar_195 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star19195zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_20, star_z_20, star_mstar_20 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star19520zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_205, star_z_205, star_mstar_205 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star20205zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_21, star_z_21, star_mstar_21 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star20521zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_215, star_z_215, star_mstar_215 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star21215zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_22, star_z_22, star_mstar_22 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star21522zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_225, star_z_225, star_mstar_225 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star22225zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_23, star_z_23, star_mstar_23 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star22523zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_235, star_z_235, star_mstar_235 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star23235zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)
star_imag_24, star_z_24, star_mstar_24 = np.loadtxt("/mfst01a/rusucs/results_kappagamma/mstarsim/star23524zcut_catpdzmstar_magrepaired.cat",usecols = [0,2,3], unpack=True)

############################
# read the binary kappa file
############################

start_readkappa = time.time()

if str(pln) in plane:
    readbinary(".kappa",plane,radiusstr)
    pos1D,kappa = np.loadtxt("kappa_values_%smed%s.dat" % (plane,radiusstr), unpack=True)
    readbinary(".gamma_1",plane,radiusstr)
    gamma1 = np.loadtxt("kappa_values_%smed%s.dat" % (plane,radiusstr), usecols = [1], unpack=True)
    readbinary(".gamma_2",plane,radiusstr)
    gamma2 = np.loadtxt("kappa_values_%smed%s.dat" % (plane,radiusstr), usecols = [1], unpack=True)
    kappagamma = np.c_[pos1D,kappa,gamma1,gamma2]
    os.system("rm kappa_values_%smed%s.dat" % (plane,radiusstr))
    print "Read kappa and shear in ", time.time() - start_readkappa, "seconds"
else: sys.exit('Wrong MS plane for this lens!!!')

############################
# prepare the galaxy catalogues
############################

start_readcat = time.time()

z_griK = np.array([])
posx_griK = np.array([])
posy_griK = np.array([])
mstar_griK = np.array([])
imag_griK = np.array([])
Mhalo_griK = np.array([])

root = plane[0:13]

for i in range(4):
    for j in range(4):
        file_griK = '%s%s_%d_%d_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griK_%s.images_forNAOJ.txt' % (rootgals,root,i,j,lens)
        if "measured" in type:
            posx__griK, posy__griK, imag__griK, z__griK, zspec__griK, mstar__griK, mstarspec__griK, mhalo__griK, mhalospec__griK = np.loadtxt(file_griK, usecols = (2,3,7,8,1,10,9,12,11), unpack=True)
        else: sys.exit('Only \"measured\" is accepted as the fourth argument of the code. Execution stopped.')

        if "measured" in type:
            spec = np.random.uniform(0,1,len(z__griK)) # is used to indicate the objects which have "spectroscopic" redshifts
            z__griK[(spec < fracspec20) & (imag__griK <= 20)] = zspec__griK[(spec < fracspec20) & (imag__griK <= 20)] # use "spectroscopic" redshifts for the randomy selected objects
            z__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)] = zspec__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)]
            z__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)] = zspec__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)]
            z__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)] = zspec__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)]
            mstar__griK[(spec < fracspec20) & (imag__griK <= 20)] = mstarspec__griK[(spec < fracspec20) & (imag__griK <= 20)] # use the corresponding stellar masses for the "spectroscopic" redshifts objects
            mstar__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)] = mstarspec__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)]
            mstar__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)] = mstarspec__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)]
            mstar__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)] = mstarspec__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)]
            mhalo__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)] = mhalospec__griK[(spec < fracspec21) & (imag__griK > 20) & (imag__griK <= 21)]
            mhalo__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)] = mhalospec__griK[(spec < fracspec22) & (imag__griK > 21) & (imag__griK <= 22)]
            mhalo__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)] = mhalospec__griK[(spec < fracspec23) & (imag__griK > 22) & (imag__griK <= 23)]

            # count the galaxies in order to implement the fraction of stars
            count18 = z__griK[(imag__griK > brightmag) & (imag__griK <= 18) & (z__griK <= z_s)].size
            count185 = z__griK[(imag__griK > 18) & (imag__griK <= 18.5) & (z__griK <= z_s)].size
            count19 = z__griK[(imag__griK > 18.5) & (imag__griK <= 19) & (z__griK <= z_s)].size
            count195 = z__griK[(imag__griK > 19) & (imag__griK <= 19.5) & (z__griK <= z_s)].size
            count20 = z__griK[(imag__griK > 19.5) & (imag__griK <= 20) & (z__griK <= z_s)].size
            count205 = z__griK[(imag__griK > 20) & (imag__griK <= 20.5) & (z__griK <= z_s)].size
            count21 = z__griK[(imag__griK > 20.5) & (imag__griK <= 21) & (z__griK <= z_s)].size
            count215 = z__griK[(imag__griK > 21) & (imag__griK <= 21.5) & (z__griK <= z_s)].size
            count22 = z__griK[(imag__griK > 21.5) & (imag__griK <= 22) & (z__griK <= z_s)].size
            count225 = z__griK[(imag__griK > 22) & (imag__griK <= 22.5) & (z__griK <= z_s)].size
            count23 = z__griK[(imag__griK > 22.5) & (imag__griK <= 23) & (z__griK <= z_s)].size
            #count235 = z__griK[(imag__griK > 23) & (imag__griK <= 23.5) & (z__griK <= z_s)].size
            #count24 = z__griK[(imag__griK > 23.5) & (imag__griK <= 24) & (z__griK <= z_s)].size
            #print count18,count185,count19,count195,count20,count205,count21,count215,count22,count225,count23,count235,count24
        
            posxmin = np.min(posx__griK) # used to insert random stellar contaminants
            posxmax = np.max(posx__griK)
            posymin = np.min(posy__griK)
            posymax = np.max(posy__griK)

            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_imag_18,cont_z_18,cont_mstar_18 = contaminants(count18,cont_griK_18,posxmin,posxmax,posymin,posymax,star_imag_18,star_z_18,star_mstar_18)
            cont_185,cont_posx_185,cont_posy_185,cont_imag_185,cont_z_185,cont_mstar_185 = contaminants(count185,cont_griK_185,posxmin,posxmax,posymin,posymax,star_imag_185,star_z_185,star_mstar_185)
            cont_19,cont_posx_19,cont_posy_19,cont_imag_19,cont_z_19,cont_mstar_19 = contaminants(count19,cont_griK_19,posxmin,posxmax,posymin,posymax,star_imag_19,star_z_19,star_mstar_19)
            cont_195,cont_posx_195,cont_posy_195,cont_imag_195,cont_z_195,cont_mstar_195 = contaminants(count195,cont_griK_195,posxmin,posxmax,posymin,posymax,star_imag_195,star_z_195,star_mstar_195)
            cont_20,cont_posx_20,cont_posy_20,cont_imag_20,cont_z_20,cont_mstar_20 = contaminants(count20,cont_griK_20,posxmin,posxmax,posymin,posymax,star_imag_20,star_z_20,star_mstar_20)
            cont_205,cont_posx_205,cont_posy_205,cont_imag_205,cont_z_205,cont_mstar_205 = contaminants(count205,cont_griK_205,posxmin,posxmax,posymin,posymax,star_imag_205,star_z_205,star_mstar_205)
            cont_21,cont_posx_21,cont_posy_21,cont_imag_21,cont_z_21,cont_mstar_21 = contaminants(count21,cont_griK_21,posxmin,posxmax,posymin,posymax,star_imag_21,star_z_21,star_mstar_21)
            cont_215,cont_posx_215,cont_posy_215,cont_imag_215,cont_z_215,cont_mstar_215 = contaminants(count215,cont_griK_215,posxmin,posxmax,posymin,posymax,star_imag_215,star_z_215,star_mstar_215)
            cont_22,cont_posx_22,cont_posy_22,cont_imag_22,cont_z_22,cont_mstar_22 = contaminants(count22,cont_griK_22,posxmin,posxmax,posymin,posymax,star_imag_22,star_z_22,star_mstar_22)
            cont_225,cont_posx_225,cont_posy_225,cont_imag_225,cont_z_225,cont_mstar_225 = contaminants(count225,cont_griK_225,posxmin,posxmax,posymin,posymax,star_imag_225,star_z_225,star_mstar_225)
            cont_23,cont_posx_23,cont_posy_23,cont_imag_23,cont_z_23,cont_mstar_23 = contaminants(count23,cont_griK_23,posxmin,posxmax,posymin,posymax,star_imag_23,star_z_23,star_mstar_23)
            #cont_235,cont_posx_235,cont_posy_235,cont_imag_235,cont_z_235,cont_mstar_235 = contaminants(count235,cont_griK_235,posxmin,posxmax,posymin,posymax,star_imag_235,star_z_235,star_mstar_235)
            #cont_24,cont_posx_24,cont_posy_24,cont_imag_24,cont_z_24,cont_mstar_24 = contaminants(count24,cont_griK_24,posxmin,posxmax,posymin,posymax,star_imag_24,star_z_24,star_mstar_24)

            #print cont_18.size,cont_185.size,cont_19.size,cont_195.size,cont_20.size,cont_205.size,cont_21.size,cont_215.size,cont_22.size,cont_225.size,cont_23.size,cont_235.size,cont_24.size
        
            # masking the fraction of galaxies expected to not be detected as galaxies
            z___griK = z__griK
            z__griK = z___griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]
            posx__griK = posx__griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]
            posy__griK = posy__griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]
            mstar__griK = mstar__griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]
            mhalo__griK = mhalo__griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]
            imag__griK = imag__griK[(imag__griK > brightmag) & (imag__griK <= 23) & (z___griK <= z_s) & (((imag__griK > brightmag) & (imag__griK <= 18) & (spec < 1 - inc_griK_18)) | ((imag__griK > 18) & (imag__griK <= 18.5) & (spec < 1 - inc_griK_185)) | ((imag__griK > 18.5) & (imag__griK <= 19) & (spec < 1 - inc_griK_19)) | ((imag__griK > 19) & (imag__griK <= 19.5) & (spec < 1 - inc_griK_195)) | ((imag__griK > 19.5) & (imag__griK <= 20) & (spec < 1 - inc_griK_20)) | ((imag__griK > 20) & (imag__griK <= 20.5) & (spec < 1 - inc_griK_205)) | ((imag__griK > 20.5) & (imag__griK <= 21) & (spec < 1 - inc_griK_21)) | ((imag__griK > 21) & (imag__griK <= 21.5) & (spec < 1 - inc_griK_215)) | ((imag__griK > 21.5) & (imag__griK <= 22) & (spec < 1 - inc_griK_22)) | ((imag__griK > 22) & (imag__griK <= 22.5) & (spec < 1 - inc_griK_225)) | ((imag__griK > 22.5) & (imag__griK <= 23) & (spec < 1 - inc_griK_23)))]

            # inserting the stellar contaminants
            z_griK = np.concatenate((z_griK,z__griK,cont_z_18,cont_z_185,cont_z_19,cont_z_195,cont_z_20,cont_z_205,cont_z_21,cont_z_215,cont_z_22,cont_z_225,cont_z_23))
            posx_griK = np.concatenate((posx_griK,posx__griK,cont_posx_18,cont_posx_185,cont_posx_19,cont_posx_195,cont_posx_20,cont_posx_205,cont_posx_21,cont_posx_215,cont_posx_22,cont_posx_225,cont_posx_23))
            posy_griK = np.concatenate((posy_griK,posy__griK,cont_posy_18,cont_posy_185,cont_posy_19,cont_posy_195,cont_posy_20,cont_posy_205,cont_posy_21,cont_posy_215,cont_posy_22,cont_posy_225,cont_posy_23))
            mstar_griK = np.concatenate((mstar_griK,mstar__griK,cont_mstar_18,cont_mstar_185,cont_mstar_19,cont_mstar_195,cont_mstar_20,cont_mstar_205,cont_mstar_21,cont_mstar_215,cont_mstar_22,cont_mstar_225,cont_mstar_23))
            Mhalo_griK = np.concatenate((Mhalo_griK,mhalo__griK,cont_mstar_18,cont_mstar_185,cont_mstar_19,cont_mstar_195,cont_mstar_20,cont_mstar_205,cont_mstar_21,cont_mstar_215,cont_mstar_22,cont_mstar_225,cont_mstar_23))
            imag_griK = np.concatenate((imag_griK,imag__griK,cont_imag_18,cont_imag_185,cont_imag_19,cont_imag_195,cont_imag_20,cont_imag_205,cont_imag_21,cont_imag_215,cont_imag_22,cont_imag_225,cont_imag_23))

            # repeat for ugriz catalogues
        
            # generate the stellar contaminants
            cont_18,cont_posx_18,cont_posy_18,cont_imag_18,cont_z_18,cont_mstar_18 = contaminants(count18,cont_griK_18,posxmin,posxmax,posymin,posymax,star_imag_18,star_z_18,star_mstar_18)
            cont_185,cont_posx_185,cont_posy_185,cont_imag_185,cont_z_185,cont_mstar_185 = contaminants(count185,cont_griK_185,posxmin,posxmax,posymin,posymax,star_imag_185,star_z_185,star_mstar_185)
            cont_19,cont_posx_19,cont_posy_19,cont_imag_19,cont_z_19,cont_mstar_19 = contaminants(count19,cont_griK_19,posxmin,posxmax,posymin,posymax,star_imag_19,star_z_19,star_mstar_19)
            cont_195,cont_posx_195,cont_posy_195,cont_imag_195,cont_z_195,cont_mstar_195 = contaminants(count195,cont_griK_195,posxmin,posxmax,posymin,posymax,star_imag_195,star_z_195,star_mstar_195)
            cont_20,cont_posx_20,cont_posy_20,cont_imag_20,cont_z_20,cont_mstar_20 = contaminants(count20,cont_griK_20,posxmin,posxmax,posymin,posymax,star_imag_20,star_z_20,star_mstar_20)
            cont_205,cont_posx_205,cont_posy_205,cont_imag_205,cont_z_205,cont_mstar_205 = contaminants(count205,cont_griK_205,posxmin,posxmax,posymin,posymax,star_imag_205,star_z_205,star_mstar_205)
            cont_21,cont_posx_21,cont_posy_21,cont_imag_21,cont_z_21,cont_mstar_21 = contaminants(count21,cont_griK_21,posxmin,posxmax,posymin,posymax,star_imag_21,star_z_21,star_mstar_21)
            cont_215,cont_posx_215,cont_posy_215,cont_imag_215,cont_z_215,cont_mstar_215 = contaminants(count215,cont_griK_215,posxmin,posxmax,posymin,posymax,star_imag_215,star_z_215,star_mstar_215)
            cont_22,cont_posx_22,cont_posy_22,cont_imag_22,cont_z_22,cont_mstar_22 = contaminants(count22,cont_griK_22,posxmin,posxmax,posymin,posymax,star_imag_22,star_z_22,star_mstar_22)
            cont_225,cont_posx_225,cont_posy_225,cont_imag_225,cont_z_225,cont_mstar_225 = contaminants(count225,cont_griK_225,posxmin,posxmax,posymin,posymax,star_imag_225,star_z_225,star_mstar_225)
            cont_23,cont_posx_23,cont_posy_23,cont_imag_23,cont_z_23,cont_mstar_23 = contaminants(count23,cont_griK_23,posxmin,posxmax,posymin,posymax,star_imag_23,star_z_23,star_mstar_23)
            #cont_235,cont_posx_235,cont_posy_235,cont_imag_235,cont_z_235,cont_mstar_235 = contaminants(count235,cont_griK_235,posxmin,posxmax,posymin,posymax,star_imag_235,star_z_235,star_mstar_235)
            #cont_24,cont_posx_24,cont_posy_24,cont_imag_24,cont_z_24,cont_mstar_24 = contaminants(count24,cont_griK_24,posxmin,posxmax,posymin,posymax,star_imag_24,star_z_24,star_mstar_24)

if "measured" in type:
    mstar_griK = 10 ** mstar_griK
    Mhalo_griK = 10 ** Mhalo_griK

cat_griK = np.c_[z_griK,posx_griK,posy_griK,mstar_griK,imag_griK,Mhalo_griK]
del z___griK
del z_griK
del posx_griK
del posy_griK
del mstar_griK
del Mhalo_griK
del imag_griK
del z__griK
del posx__griK
del posy__griK
del mstar__griK
del imag__griK

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

cat = cat_griK
bands = "griK"
weightedcounts(cat,spacing,radius,lim1D,cells_on_a_side,L_field,L_pix,cells,kappagamma,pln,type,lens,plane,bands,innermsk)

print "Computed weights in ", time.time() - start_weights, "seconds"

print(" Field done in --- %s seconds ---" % (time.time() - start_time))

print 'Done!'
