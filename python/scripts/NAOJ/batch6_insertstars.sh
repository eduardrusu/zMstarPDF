#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1
#PBS -o Logb6.out
#PBS -e Logb6.err
#PBS -N 6
#PBS -l mem=11gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_1_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_2_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_3_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_4_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_5_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_6_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_7_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5

python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_0_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_1_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_2_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_3_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_4_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_5_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_6_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
python kappamed_insertstarsnobetasingleband.py PG1115 GGL_los_8_6_7_N_4096_ang_4_rays_to_plane_34_f 23 120 measured 5
