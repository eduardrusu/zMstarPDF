#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb2.out
#PBS -e Logb2.err
#PBS -N 2
#PBS -l mem=20gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2rms
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3rms
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overrrms
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overrrms
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_flexion
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_tidal
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIS
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIShalo
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked chameleon empty notremovegroups 5 22.5 measured med 120_gal 120_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked chameleon empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked composite empty notremovegroups 5 22.5 measured med 120_gal 120_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked composite empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_gamma
#python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 


python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_z 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_massoverr 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overr 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overr 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2rms 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3rms 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overrrms 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overrrms 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_flexion 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_tidal 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIS 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIShalo 120_gal 120_gamma
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_z
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_zoverr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_massoverr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overr
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2rms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3rms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overrrms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overrrms
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_flexion
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_tidal
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIS
python inferkappa_unbiasedwithshear45and120FITSio.py WFI2033 -1.0 -1.0 removehandpicked fiducial empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIShalo
