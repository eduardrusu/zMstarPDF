# rename all folders given a rule

import os

#rootin = '/lfs08/rusucs/WFI2033/MSwghtratios/'
#rootout = '/lfs08/rusucs/WFI2033/MSwghtratios/'
#init = [x for x in os.listdir(rootin) if ('.lst' in x) and ('sampled.lst' not in x) and ('count.lst' not in x) and ('meds' in x) and ('orig' not in x)]
#init = [x for x in os.listdir(rootin) if ('arcsecinnermsk.cat' in x)]
#for i in range(len(init)):
#    indiv = [x.strip() for x in init[i].split("_")]
#    out = indiv[0]+'_'+indiv[1]+'_'+indiv[2]+'_'+indiv[3]+'_'+indiv[4]+'_'+indiv[5]+'_'+indiv[6]+'_'+indiv[7]+'_23_'+indiv[8]+'_5arcsecinner_gap_-1.0_-1.0.cat'
#    os.system('mv %s%s %s%s' %(rootin,init[i],rootout,out))

rootin = '/lfs08/rusucs/0408/MSgals/'
rootout = '/lfs08/rusucs/WFI2033/MSwghtratios/'
init = [x for x in os.listdir(rootin) if ('arcsecinnermsk.cat' in x)]
for i in range(len(init)):
    indiv = [x.strip() for x in init[i].split("_")]
    out = indiv[0]+'_'+indiv[1]+'_'+indiv[2]+'_'+indiv[3]+'_'+indiv[4]+'_'+indiv[5]+'_'+indiv[6]+'_'+indiv[7]+'_23_'+indiv[8]+'_5arcsecinner_gap_-1.0_-1.0.cat'
    os.system('mv %s%s %s%s' %(rootin,init[i],rootout,out))

#GGL_los_8_7_7_3_3_N_4096_ang_4_SA_galaxies_on_plane_27_to_63_griz_0408.images_forNAOJ.txt
