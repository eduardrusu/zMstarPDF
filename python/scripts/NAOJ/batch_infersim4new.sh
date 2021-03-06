#!/bin/sh
##PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Log4.out
#PBS -e Log4.err
#PBS -N 4
#PBS -l mem=30gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python inferkappa_unbiasedwithshear45and120FITSioincrem4.py J1206 -1.0 -1.0 removegrouphandpicked fiducial notempty 5 23 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr
python inferkappa_unbiasedwithshear45and120FITSioincrem4.py J1206 -1.0 -1.0 removegrouphandpicked fiducial notempty 5 23 measured med 120_gal 120_zoverr 45_gal 45_zoverr
python inferkappa_unbiasedwithshear45and120FITSioincrem4.py J1206 -1.0 -1.0 removegrouphandpicked fiducial notempty 5 24 measured med 120_gal 120_oneoverr 45_gal 45_oneoverr
python inferkappa_unbiasedwithshear45and120FITSioincrem4.py J1206 -1.0 -1.0 removegrouphandpicked fiducial notempty 5 24 measured med 120_gal 120_zoverr 45_gal 45_zoverr