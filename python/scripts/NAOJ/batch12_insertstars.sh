#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Logb12.out
#PBS -e Logb12.err
#PBS -N 12
#PBS -l mem=4gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_0_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_1_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_2_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_3_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_4_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_5_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_6_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_7_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 5 0.61 0.71

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_0_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_1_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_2_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_3_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_4_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_5_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_6_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_4_7_N_4096_ang_4_rays_to_plane_35_f 23 45 measured 15 -1 -1