# script the execution of inferkappasimbias.py to vary with a specific parameter

import sys
import os
#import numpy as np

file = "/lfs08/rusucs/code/inferkappasimbias.py"
for i in range(10): # 8
    constr_gal_meds45 = 0.6 + i * 0.2
    with open(file, 'r') as f:
        code = f.readlines()
    code[138 - 1] = "constr_gal_meds45 = %s \n" % constr_gal_meds45
    code[139 - 1] = "constrwidth_gal_meds_inf45 = %s \n" % (constr_gal_meds45 - 0.1)
    code[140 - 1] = "constrwidth_gal_meds_sup45 = %s \n" % (constr_gal_meds45 + 0.1)
    code[142 - 1] = "constr_gal_meds120 = %s \n" % constr_gal_meds45
    code[143 - 1] = "constrwidth_gal_meds_inf120 = %s \n" % (constr_gal_meds45 - 0.1)
    code[144 - 1] = "constrwidth_gal_meds_sup120 = %s \n" % (constr_gal_meds45 + 0.1)
    with open(file, 'w') as f:
        f.writelines(code)
        f.close()
    os.system("python inferkappasimbias.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 45_gal")
    #os.system("python inferkappasimbias.py WFI2033 -1.0 -1.0 empty notremovegroups 5 22.5 measured med 120_gal")
    #os.system("python inferkappasimbias.py WFI2033 -1.0 -1.0 empty notremovegroups 5 23.5 measured med 45_gal")
    #os.system("python inferkappasimbias.py WFI2033 -1.0 -1.0 empty notremovegroups 5 23.5 measured med 120_gal")
