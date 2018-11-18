#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb1.out
#PBS -e Logb1.err
#PBS -N 1
#PBS -l mem=15gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_1_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_2_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_3_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_4_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_5_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_6_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_7_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5

python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_0_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_1_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_2_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_3_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_4_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_5_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_6_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_1_7_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
