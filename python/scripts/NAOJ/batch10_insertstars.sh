#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1
#PBS -o Logb10.out
#PBS -e Logb10.err
#PBS -N 10
#PBS -l mem=10gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_1_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_2_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_3_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_4_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_5_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_6_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_7_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_0_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_1_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_2_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_3_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_4_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_5_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_6_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_2_7_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
