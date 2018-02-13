#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1m
#PBS -o Logb6.out
#PBS -e Logb6.err
#PBS -N 6
#PBS -l mem=4gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /mfst01a/rusucs/WFI2033/MSwghtratios/

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_0_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_1_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_2_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_3_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_4_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_5_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_6_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_7_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 5 0.61 0.71

python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_0_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_1_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_2_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_3_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_4_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_5_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_6_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
python kappamed_insertstarsnobeta.py WFI2033 GGL_los_8_6_7_N_4096_ang_4_rays_to_plane_35_f 23 120 measured 15 -1 -1
