# Use to combine a  certain kappa distribution in case it could not be computed at once from the 64 MS fields. Normalizes before summing. 

import numpy as np

root = "/Users/cerusu/Dropbox/"

data1 = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap0.61_0.71_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_23_meds_increments2_2_2_2_part1.cat" % root, usecols=[0], unpack=True)
data2 = np.loadtxt("%skappahist_WFI2033_5innermask_nobeta_zgap0.61_0.71_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_23_meds_increments2_2_2_2_part2.cat" % root, usecols=[0], unpack=True)

data1 = data1 / np.sum(data1)
data2 = data2 / np.sum(data2)

np.savetxt("%skappahist_WFI2033_5innermask_nobeta_zgap0.61_0.71_fiducial_120_gal_120_gamma_120_oneoverr_45_gal_23_meds_increments2_2_2_2.cat" % root, data1 + data2)
