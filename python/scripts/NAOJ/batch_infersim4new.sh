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

python inferkappa_unbiasedwithshearJ1206withHE0435.py J1206 -1.0 -1.0 nohandpicked fiducial 5 45 24 meds gal oneoverr
