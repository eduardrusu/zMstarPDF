#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Log2s.out
#PBS -e Log2s.err
#PBS -N 2s
#PBS -l mem=16gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python inferkappasimbias.py WFI2033 5 45 23 meds gal oneoverr
python inferkappasimbias.py WFI2033 5 120 23 meds gal oneoverr
python inferkappasimbiasphil.py WFI2033 5 45 23 meds gal oneoverr
python inferkappasimbiasphil.py WFI2033 5 120 23 meds gal oneoverr
