# Use to plot the results of mcmc_einstmagniftime.py. After running mcmc_einstmagniftime.py I need to save the output of its terminal into terminal.dat, then I can run this code

import sys
import os
import numpy as np
import corner

img = 4
fileorig = "pointSIEgamma.input"
fileoutorig = "pointSIEgamma_einstmagniftime_out_.dat"
file = "terminal.dat"
os.system("grep \"M_Sun/h\" %s > %s" % (file,fileoutorig[:-4]+"einst.dat"))
os.system("paste %s %s > %s" % (fileoutorig,fileoutorig[:-4]+"einst.dat",fileoutorig[:-5]+".dat"))
os.system("rm %s %s %s" % (fileoutorig,fileoutorig[:-4]+"einst.dat",file))

mcmc = np.loadtxt(fileoutorig[:-5]+".dat", usecols=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,26,0,1,2])
mcmc = mcmc[mcmc[:,17]==img]
mcmc = np.delete(mcmc,17,1)
mcmcfinal = np.copy(mcmc)

# make sure the order of images is correct:

for i in range(np.shape(mcmc)[0]):
    if i > 0:
        dist1_1 = np.sqrt((mcmc[i][0]-mcmc[0][0])**2 + (mcmc[i][1]-mcmc[0][1])**2)
        dist2_1 = np.sqrt((mcmc[i][4]-mcmc[0][0])**2 + (mcmc[i][5]-mcmc[0][1])**2)
        dist3_1 = np.sqrt((mcmc[i][8]-mcmc[0][0])**2 + (mcmc[i][9]-mcmc[0][1])**2)
        dist4_1 = np.sqrt((mcmc[i][12]-mcmc[0][0])**2 + (mcmc[i][13]-mcmc[0][1])**2)
        if np.min([dist1_1,dist2_1,dist3_1,dist4_1]) == dist1_1:
            pass
        if np.min([dist1_1,dist2_1,dist3_1,dist4_1]) == dist2_1:
            mcmcfinal[i][0] = mcmcfinal[i][4]; mcmcfinal[i][1] = mcmcfinal[i][5]; mcmcfinal[i][2] = mcmcfinal[i][6]; mcmcfinal[i][3] = mcmcfinal[i][7]
            mcmcfinal[i][4] = mcmc[i][0]; mcmcfinal[i][5] = mcmc[i][1]; mcmcfinal[i][6] = mcmc[i][2]; mcmcfinal[i][7] = mcmc[i][3]
        if np.min([dist1_1,dist2_1,dist3_1,dist4_1]) == dist3_1:
            mcmcfinal[i][0] = mcmcfinal[i][8]; mcmcfinal[i][1] = mcmcfinal[i][9]; mcmcfinal[i][2] = mcmcfinal[i][10]; mcmcfinal[i][3] = mcmcfinal[i][11]
            mcmcfinal[i][8] = mcmc[i][0]; mcmcfinal[i][9] = mcmc[i][1]; mcmcfinal[i][10] = mcmc[i][2]; mcmcfinal[i][11] = mcmc[i][3]
        if np.min([dist1_1,dist2_1,dist3_1,dist4_1]) == dist4_1:
            mcmcfinal[i][0] = mcmcfinal[i][12]; mcmcfinal[i][1] = mcmcfinal[i][13]; mcmcfinal[i][2] = mcmcfinal[i][14]; mcmcfinal[i][3] = mcmcfinal[i][15]
            mcmcfinal[i][12] = mcmc[i][0]; mcmcfinal[i][13] = mcmc[i][1]; mcmcfinal[i][14] = mcmc[i][2]; mcmcfinal[i][15] = mcmc[i][3]
mcmc = np.copy(mcmcfinal)
for i in range(np.shape(mcmc)[0]):
    if i > 0:
        dist2_2 = np.sqrt((mcmc[i][4]-mcmc[0][4])**2 + (mcmc[i][5]-mcmc[0][5])**2)
        dist3_2 = np.sqrt((mcmc[i][8]-mcmc[0][4])**2 + (mcmc[i][9]-mcmc[0][5])**2)
        dist4_2 = np.sqrt((mcmc[i][12]-mcmc[0][4])**2 + (mcmc[i][13]-mcmc[0][5])**2)
        if np.min([dist2_2,dist3_2,dist4_2]) == dist2_2:
            pass
        if np.min([dist2_2,dist3_2,dist4_2]) == dist3_2:
            mcmcfinal[i][4] = mcmcfinal[i][8]; mcmcfinal[i][5] = mcmcfinal[i][9]; mcmcfinal[i][6] = mcmcfinal[i][10]; mcmcfinal[i][7] = mcmcfinal[i][11]
            mcmcfinal[i][8] = mcmc[i][4]; mcmcfinal[i][9] = mcmc[i][5]; mcmcfinal[i][10] = mcmc[i][6]; mcmcfinal[i][11] = mcmc[i][7]
        if np.min([dist2_2,dist3_2,dist4_2]) == dist4_2:
            mcmcfinal[i][4] = mcmcfinal[i][12]; mcmcfinal[i][5] = mcmcfinal[i][13]; mcmcfinal[i][6] = mcmcfinal[i][14]; mcmcfinal[i][7] = mcmcfinal[i][15]
            mcmcfinal[i][12] = mcmc[i][4]; mcmcfinal[i][13] = mcmc[i][5]; mcmcfinal[i][14] = mcmc[i][6]; mcmcfinal[i][15] = mcmc[i][7]
mcmc = np.copy(mcmcfinal)
for i in range(np.shape(mcmc)[0]):
    if i > 0:
        dist3_3 = np.sqrt((mcmc[i][8]-mcmc[0][8])**2 + (mcmc[i][9]-mcmc[0][9])**2)
        dist4_3 = np.sqrt((mcmc[i][12]-mcmc[0][8])**2 + (mcmc[i][13]-mcmc[0][9])**2)
        if np.min([dist3_3,dist4_3]) == dist3_3:
            pass
        if np.min([dist3_3,dist4_3]) == dist4_3:
            mcmcfinal[i][12] = mcmcfinal[i][8]; mcmcfinal[i][13] = mcmcfinal[i][9]; mcmcfinal[i][14] = mcmcfinal[i][10]; mcmcfinal[i][15] = mcmcfinal[i][11]
            mcmcfinal[i][8] = mcmc[i][12]; mcmcfinal[i][9] = mcmc[i][13]; mcmcfinal[i][10] = mcmc[i][14]; mcmcfinal[i][11] = mcmc[i][15]

np.savetxt('%s' % (fileoutorig[:-4] + "mcmc.dat"),mcmcfinal,fmt="%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f")
os.system("rm %s" % (fileoutorig[:-5]+".dat"))
mcmc = np.loadtxt(fileoutorig[:-4] + "mcmc.dat",usecols=[2,3,6,7,10,11,14,15,16,17,18])

# now remove the column with no dynamic range (time delay = 0 for the reference image)
mcmcfinal = np.delete(mcmc,1,1)

figure = corner.corner(mcmcfinal, labels=np.linspace(1,np.shape(mcmcfinal)[0],np.shape(mcmcfinal)[0]).astype(int).tolist(),quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig(fileoutorig[:-4] + "mcmc.png", dpi=100)

