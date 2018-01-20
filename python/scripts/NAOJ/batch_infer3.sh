#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Log3.out
#PBS -e Log3.err
#PBS -N 3
#PBS -l mem=16gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python inferkappa_unbiasedwithshear.py WFI2033 5 45 23 meds gal gamma oneoverr
python inferkappa_unbiasedwithshear.py WFI2033 5 120 23 meds gal gamma oneoverr

