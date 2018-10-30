#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb1.out
#PBS -e Logb1.err
#PBS -N 1
#PBS -l mem=20gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_flexion
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_tidal
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_SIS
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_SIShalo
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_z
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_zoverr
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_massoverr
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overr
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overr

python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr 45_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removelensgrouphandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 120_gal 120_oneoverr
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 45_gal 45_oneoverr
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 120_gal 120_oneoverr 120_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 45_gal 45_oneoverr 45_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr 45_gamma
python kappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 23 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr
