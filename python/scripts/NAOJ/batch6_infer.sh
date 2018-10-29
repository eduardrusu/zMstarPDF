#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb7.out
#PBS -e Logb7.err
#PBS -N 7
#PBS -l mem=20gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_z
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_zoverr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_massoverr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2overr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3overr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2rms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3rms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2overrrms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3overrrms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_flexion
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_tidal
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_SIS
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_SIShalo
