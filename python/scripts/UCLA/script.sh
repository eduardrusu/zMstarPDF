# !/bin/sh
# -pe dc_* 96
# -l h_data=5G,h_rt=336:00:00,highp
# -l h_data=32G
# -m e
# -N test

#module load python/2.7
python test.py

