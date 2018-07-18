. /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
cd /u/flashscratch/c/cerusu

python kappamed_insertstarsnobetanomass.py J1206 GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_34_f 23 45 measured 5 -1.0 -1.0

# run as qsub -l h_data=5G,h_rt=336:00:00,highp -m abe -N [jobname] script.sh
