#!/bin/sh
#PBS -r y
#PBS -m abe
#PBS -q q4
#PBS -o Logb2_.out
#PBS -e Logb2_.err
#PBS -N 2_
#PBS -l mem=15gb
#PBS -M eduardrusu@yahoo.com

# Go to this job's working directory
cd $PBS_O_HOME

# Run your executable
cd /lfs08/rusucs/code/

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_1_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_2_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_3_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_4_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_5_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_6_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_7_N_4096_ang_4_rays_to_plane_34_f 22.5 120 measured 5 -1 -1
