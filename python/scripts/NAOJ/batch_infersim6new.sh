#!/bin/sh
##PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Log6.out
#PBS -e Log6.err
#PBS -N 6
#PBS -l mem=30gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python inferkappa_unbiasedwithshearJ1206withHE0435.py J1206 -1.0 -1.0 removegrouphandpicked fiducial 5 45 24 meds gal zoverr