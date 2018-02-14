# rename all folders given a rule

rootin = '/Volumes/LaCieSubaru/weightedcounts/WFI2033/weightedcountsWFI2033_Davis/'
rootout = '/Volumes/LaCieSubaru/weightedcounts/WFI2033/'
import os

#init = [x for x in os.listdir(rootin) if ('.lst' in x) and ('sampled.lst' not in x) and ('count.lst' not in x) and ('meds' in x) and ('orig' not in x)]
init = [x for x in os.listdir(rootin) if ('.lst' in x) and ('count.lst' in x)]
for i in range(len(init)):
    indiv = [x.strip() for x in init[i].split("_")]
    #out = indiv[0]+'_'+indiv[1]+'_'+indiv[2]+'_msk'+indiv[3]+'arcsecrad5arcsecgap_WFI2033_'+indiv[5]+'_'+indiv[6]+'_meds_zgap-1.0_-1.0_'+indiv[9]
    out = indiv[0]+'_'+indiv[1]+'_msk'+indiv[2]+'arcsecrad5arcsecgap_WFI2033_'+indiv[4]+'_'+indiv[5]+'_meds_zgap-1.0_-1.0_count.lst'
    os.system('mv %s%s %s%s' %(rootin,init[i],rootout,out))

#W1p2p3_24galphotmstar_120_WFI2033_detir_noIRAC_meds_5arcsec_count
#W2p3p1_24galphotmstar_msk120arcsecrad5arcsecgap_WFI2033_deti_noIRAC_meds_zgap-1.0_-1.0_count
