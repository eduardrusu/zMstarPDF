# rename all folders given a rule

rootin = '/lfs08/rusucs/WFI2033/MSwghtratios/'
rootout = '/lfs08/rusucs/WFI2033/MSwghtratios/'
import os

#init = [x for x in os.listdir(rootin) if ('.lst' in x) and ('sampled.lst' not in x) and ('count.lst' not in x) and ('meds' in x) and ('orig' not in x)]
init = [x for x in os.listdir(rootin) if ('arcsecinnermsk.cat' in x)]
for i in range(len(init)):
    indiv = [x.strip() for x in init[i].split("_")]
    out = indiv[0]+'_'+indiv[1]+'_'+indiv[2]+'_'+indiv[3]+'_'+indiv[4]+'_'+indiv[5]+'_'+indiv[6]+'_'+indiv[7]+'_23_'+indiv[8]+'_5arcsecinner_gap_-1.0_-1.0.cat'
    os.system('mv %s%s %s%s' %(rootin,init[i],rootout,out))

#nobeta35measuredmedinject_ugriz_WFI2033_GGL_los_8_0_0_45_5arcsecinnermsk.cat
#nobeta35measuredmedinject_ugriz_WFI2033_GGL_los_8_0_0_23_45_5arcsecinner_gap_0.61_0.71.cat
