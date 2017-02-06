# run as: python kappacoord.py /Volumes/G-RAIDStudio/simulations/lensing_simulations/GGL_los_8_7_6_N_4096_ang_4_rays_to_plane_37_f.dat

import numpy as np
import sys
import os
from os import system

degree=np.pi/180
L_field=4.0*degree
N_pix_per_dim = 4096
L_pix = L_field / N_pix_per_dim

i=0
os.system("rm %s_pos.lst" % str(sys.argv[1])[0:len(str(sys.argv[1]))-4])
with open(str(sys.argv[1])) as f:
    for line in f:
        if (line!="\n"):
            if i%1000==0:
                print i
            kappa=float(line.split()[1])
            x=1+i/4096
            y=1+i%4096
            posx = -0.5 * L_field  + (x + 0.5) * L_pix
            posy = -0.5 * L_field  + (y + 0.5) * L_pix
            g=open('%s_pos.lst' %str(sys.argv[1])[0:len(str(sys.argv[1]))-4],'a')
            g.write('%s %s %s %s %s \n' % (posx, posy, x, y, kappa))
            g.close
            i=i+1

print "Done!"
