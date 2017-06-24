# plots observed galaxy colors versus theoretical tracks from BPZ
# this file consists of multiple codes added one ofter another for each filter combination

import numpy as np
import pylab as plt
u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz.cat",usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],unpack=True)
u = u + 1.42
g = g + 0.01
r = r + 0.03
i = i + 0.00
z = z - 0.02
Y = Y - 0.09
J = J - 0.075
H = H - 0.02
K = K - 0.07
err_x = np.sqrt(J_err**2 + H_err**2)
err_y = np.sqrt(H_err**2 + K_err**2)
plt.clf()
plt.errorbar(J-H,H-K,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
sbc_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
sbc_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
el_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
el_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
el_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
scd_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
scd_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
scd_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
sb3_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
sb3_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
sb3_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
sb2_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
sb2_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
sb2_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
im_J=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.J_HAWKI.AB",usecols=[1],unpack=True)
im_H=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.H_HAWKI.AB",usecols=[1],unpack=True)
im_Ks=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.Ks_HAWKI.AB",usecols=[1],unpack=True)
sbc_J=sbc_J[:100] # in order to constraint redshift between 0 and 1
sbc_H=sbc_H[:100]
sbc_Ks=sbc_Ks[:100]
el_J=el_J[:100]
el_H=el_H[:100]
el_Ks=el_Ks[:100]
scd_J=scd_J[:100]
scd_H=scd_H[:100]
scd_Ks=scd_Ks[:100]
sb3_J=sb3_J[:100]
sb3_H=sb3_H[:100]
sb3_Ks=sb3_Ks[:100]
sb2_J=sb2_J[:100]
sb2_H=sb2_H[:100]
sb2_Ks=sb2_Ks[:100]
im_J=im_J[:100]
im_H=im_H[:100]
im_Ks=im_Ks[:100]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR JHK ONLY I NEED TO DO CONVERSION FROM AB (TEMPLATES) TO VEGA (OBSERVED)
# MOIRCS J: AB-Vega= 0.91
# MOIRCS H: AB-Vega= 1.35
# MOIRCS Ks: AB-Vega= 1.83
plt.plot(-2.5*np.log10(sbc_J/sbc_H)-0.91+1.35,-2.5*np.log10(sbc_H/sbc_Ks)-1.35+1.83,color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_J/el_H)-0.91+1.35,-2.5*np.log10(el_H/el_Ks)-1.35+1.83,color='red',label='El')
plt.plot(-2.5*np.log10(scd_J/scd_H)-0.91+1.35,-2.5*np.log10(scd_H/scd_Ks)-1.35+1.83,color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_J/sb3_H)-0.91+1.35,-2.5*np.log10(sb3_H/sb3_Ks)-1.35+1.83,color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_J/sb2_H)-0.91+1.35,-2.5*np.log10(sb2_H/sb2_Ks)-1.35+1.83,color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_J/im_H)-0.91+1.35,-2.5*np.log10(im_H/im_Ks)-1.35+1.83,color='gray',label='Im')
plt.xlabel("J-H")
plt.ylabel("H-Ks")
plt.xlim((-1,2))
plt.ylim((-1,2))
plt.legend()
plt.show()

##################################

import pylab as plt
u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz.cat",usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],unpack=True)
u = u + 1.42
g = g + 0.01
r = r + 0.03
i = i + 0.00
z = z - 0.02
Y = Y - 0.09
J = J - 0.075
H = H - 0.02
K = K - 0.07
err_x = np.sqrt(u_err**2 + g_err**2)
err_y = np.sqrt(g_err**2 + r_err**2)
plt.clf()
plt.errorbar(u-g,g-r,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.u_DEC.AB",usecols=[1],unpack=True)
sbc_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sbc_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.r_DEC.AB",usecols=[1],unpack=True)
el_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.u_DEC.AB",usecols=[1],unpack=True)
el_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.g_DEC.AB",usecols=[1],unpack=True)
el_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.r_DEC.AB",usecols=[1],unpack=True)
scd_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.u_DEC.AB",usecols=[1],unpack=True)
scd_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.g_DEC.AB",usecols=[1],unpack=True)
scd_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb3_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.u_DEC.AB",usecols=[1],unpack=True)
sb3_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb3_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb2_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.u_DEC.AB",usecols=[1],unpack=True)
sb2_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb2_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.r_DEC.AB",usecols=[1],unpack=True)
im_u=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.u_DEC.AB",usecols=[1],unpack=True)
im_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.g_DEC.AB",usecols=[1],unpack=True)
im_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sbc_u=sbc_u[:100] # in order to constraint redshift between 0 and 1
sbc_g=sbc_g[:100]
sbc_r=sbc_r[:100]
el_u=el_u[:100]
el_g=el_g[:100]
el_r=el_r[:100]
scd_u=scd_u[:100]
scd_g=scd_g[:100]
scd_r=scd_r[:100]
sb3_u=sb3_u[:100]
sb3_g=sb3_g[:100]
sb3_r=sb3_r[:100]
sb2_u=sb2_u[:100]
sb2_g=sb2_g[:100]
sb2_r=sb2_r[:100]
im_u=im_u[:100]
im_g=im_g[:100]
im_r=im_r[:100]
plt.plot(-2.5*np.log10(sbc_u/sbc_g),-2.5*np.log10(sbc_g/sbc_r),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_u/el_g),-2.5*np.log10(el_g/el_r),color='red',label='El')
plt.plot(-2.5*np.log10(scd_u/scd_g),-2.5*np.log10(scd_g/scd_r),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_u/sb3_g),-2.5*np.log10(sb3_g/sb3_r),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_u/sb2_g),-2.5*np.log10(sb2_g/sb2_r),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_u/im_g),-2.5*np.log10(im_g/im_r),color='gray',label='Im')
plt.xlabel("u-g")
plt.ylabel("g-r")
plt.xlim((-0.5,2))
plt.ylim((-0.5,2.5))
plt.legend()
plt.show()

##################################

import pylab as plt
u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz.cat",usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],unpack=True)
u = u + 1.42
g = g + 0.01
r = r + 0.03
i = i + 0.00
z = z - 0.02
Y = Y - 0.09
J = J - 0.075
H = H - 0.02
K = K - 0.07
err_x = np.sqrt(g_err**2 + r_err**2)
err_y = np.sqrt(r_err**2 + i_err**2)
plt.clf()
plt.errorbar(g-r,r-i,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sbc_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sbc_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.i_DEC.AB",usecols=[1],unpack=True)
el_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.g_DEC.AB",usecols=[1],unpack=True)
el_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.r_DEC.AB",usecols=[1],unpack=True)
el_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.i_DEC.AB",usecols=[1],unpack=True)
scd_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.g_DEC.AB",usecols=[1],unpack=True)
scd_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.r_DEC.AB",usecols=[1],unpack=True)
scd_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb3_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb3_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb3_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb2_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb2_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb2_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.i_DEC.AB",usecols=[1],unpack=True)
im_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.g_DEC.AB",usecols=[1],unpack=True)
im_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.r_DEC.AB",usecols=[1],unpack=True)
im_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sbc_g=sbc_g[:100] # in order to constraint redshift between 0 and 1
sbc_r=sbc_r[:100]
sbc_i=sbc_i[:100]
el_g=el_g[:100]
el_r=el_r[:100]
el_i=el_i[:100]
scd_g=scd_g[:100]
scd_r=scd_r[:100]
scd_i=scd_i[:100]
sb3_g=sb3_g[:100]
sb3_r=sb3_r[:100]
sb3_i=sb3_i[:100]
sb2_g=sb2_g[:100]
sb2_r=sb2_r[:100]
sb2_i=sb2_i[:100]
im_g=im_g[:100]
im_r=im_r[:100]
im_i=im_i[:100]
plt.plot(-2.5*np.log10(sbc_g/sbc_r),-2.5*np.log10(sbc_r/sbc_i),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_g/el_r),-2.5*np.log10(el_r/el_i),color='red',label='El')
plt.plot(-2.5*np.log10(scd_g/scd_r),-2.5*np.log10(scd_r/scd_i),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_g/sb3_r),-2.5*np.log10(sb3_r/sb3_i),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_g/sb2_r),-2.5*np.log10(sb2_r/sb2_i),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_g/im_r),-2.5*np.log10(im_r/im_i),color='gray',label='Im')
plt.xlabel("g-r")
plt.ylabel("r-i")
#plt.xlim((0,2.5))
#plt.ylim((-0.5,1.5))
plt.legend()
plt.show()

##################################

import pylab as plt
u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz.cat",usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],unpack=True)
u = u + 1.42
g = g + 0.01
r = r + 0.03
i = i + 0.00
z = z - 0.02
Y = Y - 0.09
J = J - 0.075
H = H - 0.02
K = K - 0.07
err_x = np.sqrt(r_err**2 + i_err**2)
err_y = np.sqrt(i_err**2 + z_err**2)
plt.clf()
plt.errorbar(r-i,i-z,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sbc_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sbc_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.z_DEC.AB",usecols=[1],unpack=True)
el_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.r_DEC.AB",usecols=[1],unpack=True)
el_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.i_DEC.AB",usecols=[1],unpack=True)
el_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.z_DEC.AB",usecols=[1],unpack=True)
scd_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.r_DEC.AB",usecols=[1],unpack=True)
scd_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.i_DEC.AB",usecols=[1],unpack=True)
scd_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb3_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb3_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb3_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb2_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb2_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb2_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.z_DEC.AB",usecols=[1],unpack=True)
im_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.r_DEC.AB",usecols=[1],unpack=True)
im_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.i_DEC.AB",usecols=[1],unpack=True)
im_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sbc_r=sbc_r[:100] # in order to constraint redshift between 0 and 1
sbc_i=sbc_i[:100]
sbc_z=sbc_z[:100]
el_r=el_r[:100]
el_i=el_i[:100]
el_z=el_z[:100]
scd_r=scd_r[:100]
scd_i=scd_i[:100]
scd_z=scd_z[:100]
sb3_r=sb3_r[:100]
sb3_i=sb3_i[:100]
sb3_z=sb3_z[:100]
sb2_r=sb2_r[:100]
sb2_i=sb2_i[:100]
sb2_z=sb2_z[:100]
im_r=im_r[:100]
im_i=im_i[:100]
im_z=im_z[:100]
plt.plot(-2.5*np.log10(sbc_r/sbc_i),-2.5*np.log10(sbc_i/sbc_z),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_r/el_i),-2.5*np.log10(el_i/el_z),color='red',label='El')
plt.plot(-2.5*np.log10(scd_r/scd_i),-2.5*np.log10(scd_i/scd_z),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_r/sb3_i),-2.5*np.log10(sb3_i/sb3_z),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_r/sb2_i),-2.5*np.log10(sb2_i/sb2_z),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_r/im_i),-2.5*np.log10(im_i/im_z),color='gray',label='Im')
plt.xlabel("r-i")
plt.ylabel("i-z")
#plt.xlim((0,2.5))
#plt.ylim((-0.5,1.5))
plt.legend()
plt.show()

##################################

import pylab as plt
u,u_err,g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err,J,J_err,H,H_err,K,K_err=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/test/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_forbpzspecz.cat",usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],unpack=True)
u = u + 1.42
g = g + 0.01
r = r + 0.03
i = i + 0.00
z = z - 0.02
Y = Y - 0.09
J = J - 0.075
H = H - 0.02
K = K - 0.07
err_x = np.sqrt(i_err**2 + z_err**2)
err_y = np.sqrt(z_err**2 + Y_err**2)
plt.clf()
plt.errorbar(i-z,z-Y,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sbc_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sbc_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
el_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.i_DEC.AB",usecols=[1],unpack=True)
el_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.z_DEC.AB",usecols=[1],unpack=True)
el_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
scd_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.i_DEC.AB",usecols=[1],unpack=True)
scd_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.z_DEC.AB",usecols=[1],unpack=True)
scd_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sb3_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb3_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb3_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sb2_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb2_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb2_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
im_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.i_DEC.AB",usecols=[1],unpack=True)
im_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.z_DEC.AB",usecols=[1],unpack=True)
im_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sbc_i=sbc_i[:100] # in order to constraint redshift between 0 and 1
sbc_z=sbc_z[:100]
sbc_Y=sbc_Y[:100]
el_i=el_i[:100]
el_z=el_z[:100]
el_Y=el_Y[:100]
scd_i=scd_i[:100]
scd_z=scd_z[:100]
scd_Y=scd_Y[:100]
sb3_i=sb3_i[:100]
sb3_z=sb3_z[:100]
sb3_Y=sb3_Y[:100]
sb2_i=sb2_i[:100]
sb2_z=sb2_z[:100]
sb2_Y=sb2_Y[:100]
im_i=im_i[:100]
im_z=im_z[:100]
im_Y=im_Y[:100]
plt.plot(-2.5*np.log10(sbc_i/sbc_z),-2.5*np.log10(sbc_z/sbc_Y),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_i/el_z),-2.5*np.log10(el_z/el_Y),color='red',label='El')
plt.plot(-2.5*np.log10(scd_i/scd_z),-2.5*np.log10(scd_z/scd_Y),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_i/sb3_z),-2.5*np.log10(sb3_z/sb3_Y),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_i/sb2_z),-2.5*np.log10(sb2_z/sb2_Y),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_i/im_z),-2.5*np.log10(im_z/im_Y),color='gray',label='Im')
plt.xlabel("i-z")
plt.ylabel("z-Y")
plt.xlim((-1,1))
plt.ylim((-0.5,1))
plt.legend()
plt.show()

##################################

import pylab as plt
g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err=np.loadtxt("DESlarge_matchDominique_gals.cat",usecols=[3,4,5,6,7,8,9,10,11,12],unpack=True)
g = g + 0.00
r = r + 0.00
i = i + 0.00
z = z + 0.00
Y = Y + 0.00
err_x = np.sqrt(g_err**2 + r_err**2)
err_y = np.sqrt(r_err**2 + i_err**2)
plt.clf()
plt.errorbar(g-r,r-i,xerr=err_x,yerr=err_y,color='black',linestyle="None")
sbc_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sbc_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sbc_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.i_DEC.AB",usecols=[1],unpack=True)
el_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.g_DEC.AB",usecols=[1],unpack=True)
el_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.r_DEC.AB",usecols=[1],unpack=True)
el_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.i_DEC.AB",usecols=[1],unpack=True)
scd_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.g_DEC.AB",usecols=[1],unpack=True)
scd_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.r_DEC.AB",usecols=[1],unpack=True)
scd_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb3_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb3_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb3_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb2_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.g_DEC.AB",usecols=[1],unpack=True)
sb2_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.r_DEC.AB",usecols=[1],unpack=True)
sb2_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.i_DEC.AB",usecols=[1],unpack=True)
im_g=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.g_DEC.AB",usecols=[1],unpack=True)
im_r=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.r_DEC.AB",usecols=[1],unpack=True)
im_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sbc_g=sbc_g[:100] # in order to constraint redshift between 0 and 1
sbc_r=sbc_r[:100]
sbc_i=sbc_i[:100]
el_g=el_g[:100]
el_r=el_r[:100]
el_i=el_i[:100]
scd_g=scd_g[:100]
scd_r=scd_r[:100]
scd_i=scd_i[:100]
sb3_g=sb3_g[:100]
sb3_r=sb3_r[:100]
sb3_i=sb3_i[:100]
sb2_g=sb2_g[:100]
sb2_r=sb2_r[:100]
sb2_i=sb2_i[:100]
im_g=im_g[:100]
im_r=im_r[:100]
im_i=im_i[:100]
plt.plot(-2.5*np.log10(sbc_g/sbc_r),-2.5*np.log10(sbc_r/sbc_i),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_g/el_r),-2.5*np.log10(el_r/el_i),color='red',label='El')
plt.plot(-2.5*np.log10(scd_g/scd_r),-2.5*np.log10(scd_r/scd_i),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_g/sb3_r),-2.5*np.log10(sb3_r/sb3_i),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_g/sb2_r),-2.5*np.log10(sb2_r/sb2_i),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_g/im_r),-2.5*np.log10(im_r/im_i),color='gray',label='Im')
plt.xlabel("g-r")
plt.ylabel("r-i")
#plt.xlim((0,2.5))
#plt.ylim((-0.5,1.5))
plt.legend()
plt.show()

##################################

import pylab as plt
g,g_err,r,r_err,i,i_err,z,z_err,Y,Y_err=np.loadtxt("DESlarge_matchDominique_gals.cat",usecols=[3,4,5,6,7,8,9,10,11,12],unpack=True)
g = g + 0.00
r = r + 0.00
i = i + 0.00
z = z + 0.00
Y = Y + 0.00
errx_x = np.sqrt(i_err**2 + z_err**2)
errx_y = np.sqrt(z_err**2 + Y_err**2)
plt.clf()
plt.errorbar(i-z,z-Y,xerr=errx_x,yerr=errx_y,color='black',linestyle="None")
sbc_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sbc_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sbc_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Sbc_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
el_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.i_DEC.AB",usecols=[1],unpack=True)
el_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.z_DEC.AB",usecols=[1],unpack=True)
el_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/El_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
scd_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.i_DEC.AB",usecols=[1],unpack=True)
scd_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.z_DEC.AB",usecols=[1],unpack=True)
scd_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Scd_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sb3_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb3_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb3_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB3_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sb2_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.i_DEC.AB",usecols=[1],unpack=True)
sb2_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.z_DEC.AB",usecols=[1],unpack=True)
sb2_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/SB2_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
im_i=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.i_DEC.AB",usecols=[1],unpack=True)
im_z=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.z_DEC.AB",usecols=[1],unpack=True)
im_Y=np.loadtxt("/Users/eduardrusu/software/bpz-1.99.3/AB/Im_B2004a.Y_DEC.AB",usecols=[1],unpack=True)
sbc_i=sbc_i[:100] # in order to constraint redshift between 0 and 1
sbc_z=sbc_z[:100]
sbc_Y=sbc_Y[:100]
el_i=el_i[:100]
el_z=el_z[:100]
el_Y=el_Y[:100]
scd_i=scd_i[:100]
scd_z=scd_z[:100]
scd_Y=scd_Y[:100]
sb3_i=sb3_i[:100]
sb3_z=sb3_z[:100]
sb3_Y=sb3_Y[:100]
sb2_i=sb2_i[:100]
sb2_z=sb2_z[:100]
sb2_Y=sb2_Y[:100]
im_i=im_i[:100]
im_z=im_z[:100]
im_Y=im_Y[:100]
plt.plot(-2.5*np.log10(sbc_i/sbc_z),-2.5*np.log10(sbc_z/sbc_Y),color='blue',label='Sbc')
plt.plot(-2.5*np.log10(el_i/el_z),-2.5*np.log10(el_z/el_Y),color='red',label='El')
plt.plot(-2.5*np.log10(scd_i/scd_z),-2.5*np.log10(scd_z/scd_Y),color='cyan',label='Scd')
plt.plot(-2.5*np.log10(sb3_i/sb3_z),-2.5*np.log10(sb3_z/sb3_Y),color='magenta',label='SB3')
plt.plot(-2.5*np.log10(sb2_i/sb2_z),-2.5*np.log10(sb2_z/sb2_Y),color='brown',label='SB2')
plt.plot(-2.5*np.log10(im_i/im_z),-2.5*np.log10(im_z/im_Y),color='gray',label='Im')
plt.xlabel("i-z")
plt.ylabel("z-Y")
plt.xlim((-1,1.5))
plt.ylim((-2.5,2))
plt.legend()
plt.show()
