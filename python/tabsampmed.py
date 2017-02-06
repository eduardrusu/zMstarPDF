import numpy as np
import os
from os import system
tab_gal=np.loadtxt('GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size120_i24_ratioquick.lst',usecols=[6], unpack=True)
tab_oneoverr=np.loadtxt('GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size120_i24_ratioquick.lst',usecols=[8], unpack=True)
tab_massoverr=np.loadtxt('GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size120_i24_ratioquick.lst',usecols=[24], unpack=True)
tab_zoverr=np.loadtxt('GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size120_i24_ratioquick.lst',usecols=[22], unpack=True)
tab_mass3=np.loadtxt('GGL_los_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63_pdzmstar_noJHKs_B1608_tab_size120_i24_ratioquick.lst',usecols=[18], unpack=True)
os.system("rm tab_24.lst")
v = np.zeros(100)
w = np.zeros(100)
x = np.zeros(100)
y = np.zeros(100)
z = np.zeros(100)
for i in range(len(tab_massoverr)/100):
    for j in range(100):
        v[j]=tab_gal[i*100+j]
        w[j]=tab_oneoverr[i*100+j]
        x[j]=tab_massoverr[i*100+j]
        y[j]=tab_zoverr[i*100+j]
        z[j]=tab_mass3[i*100+j]
    f=open('tab_24.lst','a')
    f.write('%s %s %s %s %s \n' %(np.median(v),np.median(w),np.median(x),np.median(y),np.median(z)))
    f.close()
    if i%1000==0:
        print i