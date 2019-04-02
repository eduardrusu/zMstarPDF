#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb3.out
#PBS -e Logb3.err
#PBS -N 3
#PBS -l mem=40gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python ascii_to_fits3.py
