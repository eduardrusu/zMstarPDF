#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb6.out
#PBS -e Logb6.err
#PBS -N 6
#PBS -l mem=20gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_z 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass2 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass3 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_zoverr 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_massoverr 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass2overr 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass3overr 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass2rms 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass3rms 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass2overrrms 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_mass3overrrms 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_flexion 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_tidal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_SIS 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_SIShalo 120_gamma
