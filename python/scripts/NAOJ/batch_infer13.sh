#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Log13.out
#PBS -e Log13.err
#PBS -N 13
#PBS -l mem=16gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python inferkappa_unbiasedwithshearincrement2226.py WFI2033 5 45 23 meds gal gamma oneoverr mass3rms
python inferkappa_unbiasedwithshearincrement2224.py WFI2033 5 120 23 meds gal gamma oneoverr mass3rms
