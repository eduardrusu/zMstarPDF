#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1
#PBS -o Logb17.out
#PBS -e Logb17.err
#PBS -N 17
#PBS -l mem=4gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetanomasstable.py J1206 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1

