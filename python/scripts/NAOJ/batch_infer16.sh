#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Log16.out
#PBS -e Log16.err
#PBS -N 16
#PBS -l mem=16gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python inferkappa_unbiasedwithshearincrement2226.py WFI2033 5 45 23 meds gal gamma oneoverr flexion
python inferkappa_unbiasedwithshearincrement2224.py WFI2033 5 120 23 meds gal gamma oneoverr flexion
