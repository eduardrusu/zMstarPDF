#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb5.out
#PBS -e Logb5.err
#PBS -N 5
#PBS -l mem=20gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overrrms 45_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overrrms 45_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_flexion 45_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_tidal 45_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIS 45_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIShalo 45_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 120_gal
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked chameleon empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked chameleon empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked composite empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked composite empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_gamma
