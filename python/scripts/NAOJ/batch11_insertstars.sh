#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q1
#PBS -o Logb11.out
#PBS -e Logb11.err
#PBS -N 11
#PBS -l mem=10gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_1_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_2_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_3_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_4_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_5_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_6_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_7_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1 -1

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_0_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_1_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_2_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_3_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_4_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_5_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_6_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_3_7_N_4096_ang_4_rays_to_plane_34_f 24 45 measured 5 -1 -1
