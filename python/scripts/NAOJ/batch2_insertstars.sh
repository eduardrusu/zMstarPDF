#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb2.out
#PBS -e Logb2.err
#PBS -N 2
#PBS -l mem=10gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_0_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_1_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_2_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_3_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_4_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_5_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_6_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_7_N_4096_ang_4_rays_to_plane_35_f 22.5 45 measured 5 -1 -1

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_0_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_1_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_2_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_3_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_4_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_5_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_6_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_2_7_N_4096_ang_4_rays_to_plane_35_f 22.5 120 measured 5 -1 -1
