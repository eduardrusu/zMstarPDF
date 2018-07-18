#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1
#PBS -o Logb13.out
#PBS -e Logb13.err
#PBS -N 13
#PBS -l mem=4gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_1_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_2_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_3_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_4_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_5_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_6_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_7_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_0_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_1_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_2_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_3_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_4_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_5_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_6_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_5_7_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
