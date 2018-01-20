# This code checks, for a given file in the simulation catalogue, if photoz has been computed successfuly for each of the i<24 galaxies in the catalogue
# run from the /code folder as: python photozMilleniumcheckpresence.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_8_7_7_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt /Volumes/G-RAIDStudio/simulations/lensing_simulations/Guo_galaxies/GGL_los_8_7_7_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.pdz nr
# DO NOT RUN IT IN THE SAME TIME WITH photozMilleniumcheckpresence.py because they use the same count.lst files
# where nr is just a number, wchich should be different for each concurrent execution

import sys
import os
from os import system
import time

start_timefield = time.time()

os.system("wc %s > count%s.lst" % (str(sys.argv[1]),str(sys.argv[3]))) # need to know how many lines in file 1
with open('count%s.lst' % str(sys.argv[3])) as count:
    for line in count:
        size=line.split()[0]
pos=0 # position in file 2
posremember=0
#print int(size)
ord=0 # position in file 1
with open(str(sys.argv[1])) as search_for:
    for gal_for in search_for:
        ord=ord+1
        if gal_for!="\n":
            if gal_for.split()[0]!="GalID" and (float(gal_for.split()[15]) <= 24):
                gal_for_ID = gal_for.split()[0]
                present="false"
                with open(str(sys.argv[2])) as search_in:
                    for i in xrange(pos):
                        search_in.next()
                    for gal_in in search_in:
                      pos=pos+1
                      if gal_in!="\n":
                          if gal_in.split()[0]==gal_for_ID:
                              present="true"
                              posremember=pos
                              break
                        #print pos
                if (pos==int(size)) and (present=="false"): #need the 2nd condition, or the last element is considered not found
                    pos=posremember
                    print gal_for_ID, "not found!"
        if ord in [1000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000, 220000, 230000, 240000, 250000, 260000, 270000, 280000, 290000, 300000]:
            print ord, "objects..."
os.system("rm count%s.lst" % str(sys.argv[3]))

print("Total time for field: --- %s seconds ---" % (time.time() - start_timefield))

print 'Done!'
