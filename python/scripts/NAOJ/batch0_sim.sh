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

python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 45_gal 45_gamma 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 45_gal 45_gamma 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 45_gal 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 45_gal 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_zoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_zoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_zoverr 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_zoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 120_gal 120_gamma 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_z 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_massoverr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2rms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3rms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass2overrrms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_mass3overrrms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_flexion 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_tidal 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIS 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_zoverr 45_SIShalo 120_gal 120_gamma

python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_gamma 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_gamma 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 45_gal 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_oneoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 120_oneoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 45_oneoverr 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_oneoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_oneoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 120_gal 120_gamma 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_z 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_zoverr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_massoverr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overr 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2rms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3rms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass2overrrms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_mass3overrrms 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_flexion 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_tidal 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIS 120_gal 120_gamma
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_oneoverr 45_SIShalo 120_gal 120_gamma

python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 120_gal 120_zoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 45_gal 45_zoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 45_gamma 120_gal 120_z 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 120_gamma 45_gal 45_z 45_SIS

python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_z 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_zoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_z 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_zoverr 45_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_mass3overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal 120_gal 120_oneoverr 120_SIShalo
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_z
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_zoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass2
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass3
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_massoverr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass2overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass3overr
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass2rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass3rms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass2overrrms
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_mass3overrrms 
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_flexion
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_tidal
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_SIS
python inferkappasimbiasactual.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal 45_gal 45_oneoverr 45_SIShalo
