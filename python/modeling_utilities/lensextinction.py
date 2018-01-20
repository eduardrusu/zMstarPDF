# For a given pair of images in multiple bands, estimate A_V, R_V and lens redshift by minimizing the quantity in Falco et al. 1999

# Observables: 8 filters for one pair of images = 8; params: z, R_V, A_V = 3 -> dof = 5
# Observables: 8 filters for 6 pairs of images = 48; params: z, R_V, A_V (6), M (6, but only 3 independent) = 11 -> dof = 37
# Observables: 8 filters for 4 pairs of images = 48; params: z, R_V, A_V (4), M (only 3 independent) = 9 -> dof = 23
# Observables: 8 filters for 3 pairs of images = 24; params: z, R_V, A_V (3), M (2) = 7 -> dof = 17

import numpy as np
import extinction

w_g = 4866.
w_r = 6215.
w_i = 7545.
w_z = 8679.
w_y = 9633.
w_Y = 10200.
w_J = 12520.
w_K = 21470.

SIS_ba = 1.45
SIS_ca = 1.54
SIS_da = 0.52
SIS_bc = 1.05
SIS_db = 0.36
SIS_dc = 0.38
SIE_ba = 1.47
SIE_ca = 1.39
SIE_da = 0.58
SIE_bc = 1.05
SIE_db = 0.39
SIE_dc = 0.42

wave = np.array([w_g,w_r,w_i,w_z,w_y,w_Y,w_J,w_K])
#print extinction.ccm89(wave/(1+1), 0.05, 1.1)
#print extinction.ccm89(wave/(1+1), 0.05, -1.1)
#print extinction.odonnell94(wave, 1.0, 3.1)

data = np.loadtxt('fluxratioerr.cat', dtype={'names': ('filter', 'B/A', 'B/Ae', 'C/A', 'C/Ae', 'D/A', 'D/Ae', 'B/C', 'B/Ce', 'D/B', 'D/Be', 'D/C', 'D/Ce'),'formats': ('S2', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
f = [x[7] for x in data]
f_e = [x[8] for x in data]
SIS = SIS_bc
#print 2.5 * np.log10(SIS)
#for i in range(8):
#print -2.5 * np.log10([x[11] for x in data]) + 2.5 * np.log10(SIS_dc)

#z = 0.9
#A_V = 0.05
R_V = 3.1

def chi2(f,f_e,w,SIS,z,A,R):
    return (-2.5 * np.log10(f) + 2.5 * np.log10(SIS) - extinction.ccm89(np.array([w/(1+abs(z))]), A, R))**2/(2.5*f_e/(f*np.log(10)))**2

def chi2sumtwo(arg):
    A_V,R_V = arg
    return chi2(f[0],f_e[0],w_g,SIS,z,A_V,R_V) + chi2(f[1],f_e[1],w_r,SIS,z,A_V,R_V) + chi2(f[2],f_e[2],w_i,SIS,z,A_V,R_V) + chi2(f[3],f_e[3],w_z,SIS,z,A_V,R_V) + chi2(f[4],f_e[4],w_y,SIS,z,A_V,R_V) + chi2(f[5],f_e[5],w_Y,SIS,z,A_V,R_V) + chi2(f[6],f_e[6],w_J,SIS,z,A_V,R_V) + chi2(f[7],f_e[7],w_K,SIS,z,A_V,R_V)

def chi2sumthree(arg):
    z,A_V,R_V = arg
    #print z
    return chi2(f[0],f_e[0],w_g,SIS,z,A_V,R_V) + chi2(f[1],f_e[1],w_r,SIS,z,A_V,R_V) + chi2(f[2],f_e[2],w_i,SIS,z,A_V,R_V) + chi2(f[3],f_e[3],w_z,SIS,z,A_V,R_V) + chi2(f[4],f_e[4],w_y,SIS,z,A_V,R_V) + chi2(f[5],f_e[5],w_Y,SIS,z,A_V,R_V) + chi2(f[6],f_e[6],w_J,SIS,z,A_V,R_V) + chi2(f[7],f_e[7],w_K,SIS,z,A_V,R_V)

def chi2sumfixR(arg):
    z,A_V = arg
    return chi2(f[0],f_e[0],w_g,SIS,z,A_V,R_V) + chi2(f[1],f_e[1],w_r,SIS,z,A_V,R_V) + chi2(f[2],f_e[2],w_i,SIS,z,A_V,R_V) + chi2(f[3],f_e[3],w_z,SIS,z,A_V,R_V) + chi2(f[4],f_e[4],w_y,SIS,z,A_V,R_V) + chi2(f[5],f_e[5],w_Y,SIS,z,A_V,R_V) + chi2(f[6],f_e[6],w_J,SIS,z,A_V,R_V) + chi2(f[7],f_e[7],w_K,SIS,z,A_V,R_V)

def chi2sumfixz(arg):
    A_V,R_V = arg
    return chi2(f[0],f_e[0],w_g,SIS,z,A_V,R_V) + chi2(f[1],f_e[1],w_r,SIS,z,A_V,R_V) + chi2(f[2],f_e[2],w_i,SIS,z,A_V,R_V) + chi2(f[3],f_e[3],w_z,SIS,z,A_V,R_V) + chi2(f[4],f_e[4],w_y,SIS,z,A_V,R_V) + chi2(f[5],f_e[5],w_Y,SIS,z,A_V,R_V) + chi2(f[6],f_e[6],w_J,SIS,z,A_V,R_V) + chi2(f[7],f_e[7],w_K,SIS,z,A_V,R_V)

def chi2sumfixA(arg):
    z,R_V = arg
    return chi2(f[0],f_e[0],w_g,SIS,z,A_V,R_V) + chi2(f[1],f_e[1],w_r,SIS,z,A_V,R_V) + chi2(f[2],f_e[2],w_i,SIS,z,A_V,R_V) + chi2(f[3],f_e[3],w_z,SIS,z,A_V,R_V) + chi2(f[4],f_e[4],w_y,SIS,z,A_V,R_V) + chi2(f[5],f_e[5],w_Y,SIS,z,A_V,R_V) + chi2(f[6],f_e[6],w_J,SIS,z,A_V,R_V) + chi2(f[7],f_e[7],w_K,SIS,z,A_V,R_V)

# start optimization from random places
vary = 100
vary_z = np.random.uniform(0, 2.7, vary)
vary_A_V = np.random.uniform(0, 0.5, vary)
vary_R_V = np.random.uniform(0.5, 10, vary)

chimin = 0
chimin1 = 0
chimin2 = 0
chimin3 = 0
for i in range(vary):
    from scipy.optimize import minimize
    optimal = minimize(chi2sumthree,[vary_z[i],vary_A_V[i],vary_R_V[i]],method='Nelder-Mead',options={'maxiter': 2000})
    if i == 0:
        chimin = chi2sumthree([optimal.x[0],optimal.x[1],optimal.x[2]])[0]
    else:
        if chi2sumthree([optimal.x[0],optimal.x[1],optimal.x[2]])[0] < chimin:
            chimin = chi2sumthree([optimal.x[0],optimal.x[1],optimal.x[2]])[0]
            chimin1 = optimal.x[0]
            chimin2 = optimal.x[1]
            chimin3 = optimal.x[2]
    print i
print chimin1,chimin2,chimin3,chimin

# vary around the solution to get confidence levels
x = np.linspace(0.,1,100)
for i in range(0):
    R_V = x[i]
    optimal = minimize(chi2sumfixR,[chimin1, chimin2],method='Nelder-Mead',options={'maxiter': 2000})
    print R_V,chi2sumfixR([optimal.x[0],optimal.x[1]])

x = np.linspace(0.,0.008,100)
for i in range(0):
    A_V = x[i]
    optimal = minimize(chi2sumfixA,[chimin1, chimin3],method='Nelder-Mead',options={'maxiter': 2000})
    print A_V,chi2sumfixA([optimal.x[0],optimal.x[1]])

x = np.linspace(2.,2.3,100)
for i in range(0):
    z = x[i]
    optimal = minimize(chi2sumfixz,[chimin2, chimin3],method='Nelder-Mead',options={'maxiter': 2000})
    print z,chi2sumfixz([optimal.x[0],optimal.x[1]])

############################
# use all possible image pairs

fba = [x[1] for x in data]
fba_e = [x[2] for x in data]
fca = [x[3] for x in data]
fca_e = [x[4] for x in data]
fda = [x[5] for x in data]
fda_e = [x[6] for x in data]
fbc = [x[7] for x in data]
fbc_e = [x[8] for x in data]
fdb = [x[9] for x in data]
fdb_e = [x[10] for x in data]
fdc = [x[11] for x in data]
fdc_e = [x[12] for x in data]

def chi2sumall(arg):
    z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,A_Vdb,A_Vdc,Mba,Mca,Mda = arg # because the other magnifications can be derived from these
    return chi2(fba[0],fba_e[0],w_g,Mba,z,A_Vba,R_V) + chi2(fba[1],fba_e[1],w_r,Mba,z,A_Vba,R_V) + chi2(fba[2],fba_e[2],w_i,Mba,z,A_Vba,R_V) + chi2(fba[3],fba_e[3],w_z,Mba,z,A_Vba,R_V) + chi2(fba[4],fba_e[4],w_y,Mba,z,A_Vba,R_V) + chi2(fba[5],fba_e[5],w_Y,Mba,z,A_Vba,R_V) + chi2(fba[6],fba_e[6],w_J,Mba,z,A_Vba,R_V) + chi2(fba[7],fba_e[7],w_K,Mba,z,A_Vba,R_V) + chi2(fca[0],fca_e[0],w_g,Mca,z,A_Vca,R_V) + chi2(fca[1],fca_e[1],w_r,Mca,z,A_Vca,R_V) + chi2(fca[2],fca_e[2],w_i,Mca,z,A_Vca,R_V) + chi2(fca[3],fca_e[3],w_z,Mca,z,A_Vca,R_V) + chi2(fca[4],fca_e[4],w_y,Mca,z,A_Vca,R_V) + chi2(fca[5],fca_e[5],w_Y,Mca,z,A_Vca,R_V) + chi2(fca[6],fca_e[6],w_J,Mca,z,A_Vca,R_V) + chi2(fca[7],fca_e[7],w_K,Mca,z,A_Vca,R_V) + chi2(fda[0],fda_e[0],w_g,Mda,z,A_Vda,R_V) + chi2(fda[1],fda_e[1],w_r,Mda,z,A_Vda,R_V) + chi2(fda[2],fda_e[2],w_i,Mda,z,A_Vda,R_V) + chi2(fda[3],fda_e[3],w_z,Mda,z,A_Vda,R_V) + chi2(fda[4],fda_e[4],w_y,Mda,z,A_Vda,R_V) + chi2(fda[5],fda_e[5],w_Y,Mda,z,A_Vda,R_V) + chi2(fda[6],fda_e[6],w_J,Mda,z,A_Vda,R_V) + chi2(fda[7],fda_e[7],w_K,Mda,z,A_Vda,R_V) + chi2(fbc[0],fbc_e[0],w_g,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[1],fbc_e[1],w_r,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[2],fbc_e[2],w_i,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[3],fbc_e[3],w_z,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[4],fbc_e[4],w_y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[5],fbc_e[5],w_Y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[6],fbc_e[6],w_J,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[7],fbc_e[7],w_K,Mba/Mca,z,A_Vbc,R_V) + chi2(fdb[0],fdb_e[0],w_g,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[1],fdb_e[1],w_r,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[2],fdb_e[2],w_i,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[3],fdb_e[3],w_z,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[4],fdb_e[4],w_y,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[5],fdb_e[5],w_Y,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[6],fdb_e[6],w_J,Mda/Mba,z,A_Vdb,R_V) + chi2(fdb[7],fdb_e[7],w_K,Mda/Mba,z,A_Vdb,R_V) + chi2(fdc[0],fdc_e[0],w_g,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[1],fdc_e[1],w_r,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[2],fdc_e[2],w_i,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[3],fdc_e[3],w_z,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[4],fdc_e[4],w_y,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[5],fdc_e[5],w_Y,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[6],fdc_e[6],w_J,Mda/Mca,z,A_Vdc,R_V) + chi2(fdc[7],fdc_e[7],w_K,Mda/Mca,z,A_Vdc,R_V)

def chi2sumnodbdc(arg):
    z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,Mba,Mca,Mda = arg # because the other magnifications can be derived from these
    return chi2(fba[0],fba_e[0],w_g,Mba,z,A_Vba,R_V) + chi2(fba[1],fba_e[1],w_r,Mba,z,A_Vba,R_V) + chi2(fba[2],fba_e[2],w_i,Mba,z,A_Vba,R_V) + chi2(fba[3],fba_e[3],w_z,Mba,z,A_Vba,R_V) + chi2(fba[4],fba_e[4],w_y,Mba,z,A_Vba,R_V) + chi2(fba[5],fba_e[5],w_Y,Mba,z,A_Vba,R_V) + chi2(fba[6],fba_e[6],w_J,Mba,z,A_Vba,R_V) + chi2(fba[7],fba_e[7],w_K,Mba,z,A_Vba,R_V) + chi2(fca[0],fca_e[0],w_g,Mca,z,A_Vca,R_V) + chi2(fca[1],fca_e[1],w_r,Mca,z,A_Vca,R_V) + chi2(fca[2],fca_e[2],w_i,Mca,z,A_Vca,R_V) + chi2(fca[3],fca_e[3],w_z,Mca,z,A_Vca,R_V) + chi2(fca[4],fca_e[4],w_y,Mca,z,A_Vca,R_V) + chi2(fca[5],fca_e[5],w_Y,Mca,z,A_Vca,R_V) + chi2(fca[6],fca_e[6],w_J,Mca,z,A_Vca,R_V) + chi2(fca[7],fca_e[7],w_K,Mca,z,A_Vca,R_V) + chi2(fda[0],fda_e[0],w_g,Mda,z,A_Vda,R_V) + chi2(fda[1],fda_e[1],w_r,Mda,z,A_Vda,R_V) + chi2(fda[2],fda_e[2],w_i,Mda,z,A_Vda,R_V) + chi2(fda[3],fda_e[3],w_z,Mda,z,A_Vda,R_V) + chi2(fda[4],fda_e[4],w_y,Mda,z,A_Vda,R_V) + chi2(fda[5],fda_e[5],w_Y,Mda,z,A_Vda,R_V) + chi2(fda[6],fda_e[6],w_J,Mda,z,A_Vda,R_V) + chi2(fda[7],fda_e[7],w_K,Mda,z,A_Vda,R_V) + chi2(fbc[0],fbc_e[0],w_g,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[1],fbc_e[1],w_r,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[2],fbc_e[2],w_i,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[3],fbc_e[3],w_z,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[4],fbc_e[4],w_y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[5],fbc_e[5],w_Y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[6],fbc_e[6],w_J,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[7],fbc_e[7],w_K,Mba/Mca,z,A_Vbc,R_V)

def chi2sumnodadbdc(arg):
    z,R_V,A_Vba,A_Vca,A_Vbc,Mba,Mca = arg # because the other magnifications can be derived from these
    return chi2(fba[0],fba_e[0],w_g,Mba,z,A_Vba,R_V) + chi2(fba[1],fba_e[1],w_r,Mba,z,A_Vba,R_V) + chi2(fba[2],fba_e[2],w_i,Mba,z,A_Vba,R_V) + chi2(fba[3],fba_e[3],w_z,Mba,z,A_Vba,R_V) + chi2(fba[4],fba_e[4],w_y,Mba,z,A_Vba,R_V) + chi2(fba[5],fba_e[5],w_Y,Mba,z,A_Vba,R_V) + chi2(fba[6],fba_e[6],w_J,Mba,z,A_Vba,R_V) + chi2(fba[7],fba_e[7],w_K,Mba,z,A_Vba,R_V) + chi2(fca[0],fca_e[0],w_g,Mca,z,A_Vca,R_V) + chi2(fca[1],fca_e[1],w_r,Mca,z,A_Vca,R_V) + chi2(fca[2],fca_e[2],w_i,Mca,z,A_Vca,R_V) + chi2(fca[3],fca_e[3],w_z,Mca,z,A_Vca,R_V) + chi2(fca[4],fca_e[4],w_y,Mca,z,A_Vca,R_V) + chi2(fca[5],fca_e[5],w_Y,Mca,z,A_Vca,R_V) + chi2(fca[6],fca_e[6],w_J,Mca,z,A_Vca,R_V) + chi2(fca[7],fca_e[7],w_K,Mca,z,A_Vca,R_V) + chi2(fbc[0],fbc_e[0],w_g,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[1],fbc_e[1],w_r,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[2],fbc_e[2],w_i,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[3],fbc_e[3],w_z,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[4],fbc_e[4],w_y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[5],fbc_e[5],w_Y,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[6],fbc_e[6],w_J,Mba/Mca,z,A_Vbc,R_V) + chi2(fbc[7],fbc_e[7],w_K,Mba/Mca,z,A_Vbc,R_V)


# start optimization from random places
vary = 0
vary_z = np.random.uniform(0, 2.7, vary)
vary_R_V = np.random.uniform(0.5, 10, vary)
vary_A_Vba = np.random.uniform(0, 0.5, vary)
vary_A_Vca = np.random.uniform(0, 0.5, vary)
vary_A_Vda = np.random.uniform(0, 0.5, vary)
vary_A_Vbc = np.random.uniform(0, 0.5, vary)
vary_A_Vdb = np.random.uniform(0, 0.5, vary)
vary_A_Vdc = np.random.uniform(0, 0.5, vary)
vary_Mba = np.random.uniform(1, 2, vary)
vary_Mca = np.random.uniform(1, 2, vary)
vary_Mda = np.random.uniform(0.1, 0.4, vary)

chimin = 0
chimin_z = 0
chimin_R_V = 0
chimin_A_Vba = 0
chimin_A_Vca = 0
chimin_A_Vda = 0
chimin_A_Vbc = 0
chimin_A_Vdb = 0
chimin_A_Vdc = 0
chimin_Mba = 0
chimin_Mca = 0
chimin_Mda = 0
for i in range(vary):
    from scipy.optimize import minimize
    optimal = minimize(chi2sumall,[vary_z[i],vary_R_V[i],vary_A_Vba[i],vary_A_Vca[i],vary_A_Vda[i],vary_A_Vbc[i],vary_A_Vdb[i],vary_A_Vdc[i],vary_Mba[i],vary_Mca[i],vary_Mda[i]],method='Nelder-Mead',options={'maxiter': 2000})
    if i == 0:
        chimin = chi2sumall([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8],optimal.x[9],optimal.x[10]])[0]
    else:
        if chi2sumall([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8],optimal.x[9],optimal.x[10]])[0] < chimin:
            chimin = chi2sumall([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8],optimal.x[9],optimal.x[10]])[0]
            chimin_z = optimal.x[0]
            chimin_R_V = optimal.x[1]
            chimin_A_Vba = optimal.x[2]
            chimin_A_Vca = optimal.x[3]
            chimin_A_Vda = optimal.x[4]
            chimin_A_Vbc = optimal.x[5]
            chimin_A_Vdb = optimal.x[6]
            chimin_A_Vdc = optimal.x[7]
            chimin_Mba = optimal.x[8]
            chimin_Mca = optimal.x[9]
            chimin_Mda = optimal.x[10]
    print i
print chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vda,chimin_A_Vbc,chimin_A_Vdb,chimin_A_Vdc,chimin_Mba,chimin_Mca,chimin_Mda,chimin

vary = 0
vary_z = np.random.uniform(0, 2.7, vary)
vary_R_V = np.random.uniform(0.5, 10, vary)
vary_A_Vba = np.random.uniform(0, 0.5, vary)
vary_A_Vca = np.random.uniform(0, 0.5, vary)
vary_A_Vda = np.random.uniform(0, 0.5, vary)
vary_A_Vbc = np.random.uniform(0, 0.5, vary)
vary_Mba = np.random.uniform(1, 2, vary)
vary_Mca = np.random.uniform(1, 2, vary)
vary_Mda = np.random.uniform(0.1, 0.4, vary)

for i in range(vary):
    from scipy.optimize import minimize
    optimal = minimize(chi2sumnodbdc,[vary_z[i],vary_R_V[i],vary_A_Vba[i],vary_A_Vca[i],vary_A_Vda[i],vary_A_Vbc[i],vary_Mba[i],vary_Mca[i],vary_Mda[i]],method='Nelder-Mead',options={'maxiter': 2000})
    if i == 0:
        chimin = chi2sumnodbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8]])[0]
    else:
        if chi2sumnodbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8]])[0] < chimin:
            chimin = chi2sumnodbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6],optimal.x[7],optimal.x[8]])[0]
            chimin_z = optimal.x[0]
            chimin_R_V = optimal.x[1]
            chimin_A_Vba = optimal.x[2]
            chimin_A_Vca = optimal.x[3]
            chimin_A_Vda = optimal.x[4]
            chimin_A_Vbc = optimal.x[5]
            chimin_Mba = optimal.x[6]
            chimin_Mca = optimal.x[7]
            chimin_Mda = optimal.x[8]
    print i
print chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vda,chimin_A_Vbc,chimin_Mba,chimin_Mca,chimin_Mda,chimin

vary = 0
vary_z = np.random.uniform(0, 2.7, vary)
vary_R_V = np.random.uniform(0.5, 10, vary)
vary_A_Vba = np.random.uniform(0, 0.5, vary)
vary_A_Vca = np.random.uniform(0, 0.5, vary)
vary_A_Vbc = np.random.uniform(0, 0.5, vary)
vary_Mba = np.random.uniform(1, 2, vary)
vary_Mca = np.random.uniform(1, 2, vary)

for i in range(vary):
    from scipy.optimize import minimize
    optimal = minimize(chi2sumnodadbdc,[vary_z[i],vary_R_V[i],vary_A_Vba[i],vary_A_Vca[i],vary_A_Vbc[i],vary_Mba[i],vary_Mca[i]],method='Nelder-Mead',options={'maxiter': 2000})
    if i == 0:
        chimin = chi2sumnodadbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6]])[0]
    else:
        if chi2sumnodadbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6]])[0] < chimin:
            chimin = chi2sumnodadbdc([optimal.x[0],optimal.x[1],optimal.x[2],optimal.x[3],optimal.x[4],optimal.x[5],optimal.x[6]])[0]
            chimin_z = optimal.x[0]
            chimin_R_V = optimal.x[1]
            chimin_A_Vba = optimal.x[2]
            chimin_A_Vca = optimal.x[3]
            chimin_A_Vbc = optimal.x[4]
            chimin_Mba = optimal.x[5]
            chimin_Mca = optimal.x[6]
    print i
print chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vbc,chimin_Mba,chimin_Mca,chimin

################
#MCMC

import emcee
import time

multithread = 1
if multithread == 1:
    def log_prior(arg):
        z,A_V,R_V = arg
        if z < 0 or z > 2.7 or A_V < 0 or R_V <= 0 or R_V >= 10:
            return -np.inf # log(0)
        return 0 #-0.5 * ((R_V - 3.1)**2)/(0.5**2)
    def log_like(arg):
        z,A_V,R_V = arg
        return -0.5 * chi2sumthree(arg)
    def log_posterior(arg):
        if not np.isfinite(log_prior(arg)):
            return -np.inf
        return log_prior(arg) + log_like(arg)

    ndim = 3 # number of parameters in the model
    nwalkers = 20 # number of MCMC walkers
    nsteps = 10000 # number of MCMC steps to take including burn-in
    # initial guesses in a small ball around the best-fit
    #chimin1 = 0.2
    #chimin2 = 0.12
    #chimin3 = 3.1
    starting_guesses = np.c_[abs(np.random.normal(chimin1,0.1*chimin1+1e-05,nwalkers)),abs(np.random.normal(chimin2,0.1*chimin2+1e-05,nwalkers)),abs(np.random.normal(chimin3,0.1*chimin3+1e-05,nwalkers))] # use this in case I want to run emcee instead of parallel tempering
    print("Running MCMC...")
    start_timeemcee8 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior,threads = 7)
    # run while showing the progress
    for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations*
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Time for emcee with 7 threads: --- %s seconds ---" % (time.time() - start_timeemcee8))

    # plot the time laps
    import corner
    import pylab as plt
    from matplotlib.ticker import MaxNLocator
    plt.clf()
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].axhline(chimin1, color="#888888", lw=2)
    axes[0].set_ylabel("1")
    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].axhline(chimin2, color="#888888", lw=2)
    axes[1].set_ylabel("2")
    axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
    axes[2].yaxis.set_major_locator(MaxNLocator(5))
    axes[2].axhline(chimin3, color="#888888", lw=2)
    axes[2].set_ylabel("3")
    axes[2].set_xlabel("step number")
    fig.tight_layout(h_pad=0.0)
    fig.show()

    # print diagnostics and do corner plot
    nburn = nsteps/2 # "burn-in" to stabilize chains
    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
    #alpha_samp = sampler.flatchain.T[0]
    #beta_samp = sampler.flatchain.T[1]
    #sigma_samp = sampler.flatchain.T[2]
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
    print "median, std 1: ", np.median(samples[:,0]), np.std(samples[:,0])
    print "median, std 2: ", np.median(samples[:,1]), np.std(samples[:,1])
    print "median, std 3: ", np.median(samples[:,2]), np.std(samples[:,2])
    print chi2sumthree([np.median(samples[:,0]),np.median(samples[:,1]),np.median(samples[:,2])])[0]
    fig = corner.corner(samples, labels=["1", "2", "3",],truths=[chimin1,chimin2,chimin3])
    fig.show()

pt = 0
if pt == 1:
# parallel tempering; takes longer to run because it runs emcee for each temperature; THIS DOES NOT WORK CURRENTLY BECAUSE EMCEE NO LONGER INCLUDED PTSAMPLER
    start_timept = time.time()
    ntemps = 5
    starting_guesses = np.zeros((ntemps,nwalkers,ndim))
    for i in range(ntemps):
        starting_guesses[i] = np.c_[np.c_[abs(np.random.normal(chimin1,0.1*chimin1+1e-05,nwalkers)),abs(np.random.normal(chimin2,0.1*chimin2+1e-05,nwalkers)),abs(np.random.normal(chimin3,0.1*chimin3+1e-05,nwalkers))]]
    import ptemcee
    sampler = Sampler(ntemps,nwalkers,ndim,log_like,log_prior,threads = 7)
    print("Running MCMC...")
    for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations*
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Time for PT with 7 threads: --- %s seconds ---" % (time.time() - start_timept))
    print "Evidence: ", sampler.thermodynamic_integration_log_evidence()
    # plot the time laps
    plt.clf()
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, :, 0].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].axhline(chimin1, color="#888888", lw=2)
    axes[0].set_ylabel("1")
    axes[1].plot(sampler.chain[:, :, :, 1].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].axhline(chimin2, color="#888888", lw=2)
    axes[1].set_ylabel("2")
    axes[2].plot(sampler.chain[:, :, :, 2].T.reshape((nwalkers*ntemps, nsteps)), color="k", alpha=0.4)
    axes[2].yaxis.set_major_locator(MaxNLocator(5))
    axes[2].axhline(chimin3, color="#888888", lw=2)
    axes[2].set_ylabel("3")
    axes[2].set_xlabel("step number")
    fig.tight_layout(h_pad=0.0)
    fig.show()
    plt.clf()
    samples = sampler.chain[:, :, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
    fig = corner.corner(samples, labels=["1", "2", "3"],truths=[chimin1,chimin2,chimin3])
    fig.show()
    emcee.autocorr.integrated_time(samples,axis=0)
    emcee.autocorr.function(samples,axis=0)

all = 0
if all == 1:
    
#vary_Mba = np.random.uniform(1, 2, vary)
#vary_Mca = np.random.uniform(1, 2, vary)
#vary_Mda = np.random.uniform(0.1, 0.4, vary)

    def log_prior(arg):
        z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,A_Vdb,A_Vdc,Mba,Mca,Mda = arg
        if z < 0 or z > 2.7 or R_V <= 0 or R_V >= 10 or A_Vba < 0 or A_Vca < 0 or A_Vda < 0 or A_Vbc < 0 or A_Vdb < 0 or A_Vdc < 0:
            return -np.inf # log(0)
        return 0 #-0.5 * ((R_V - 3.1)**2)/(0.5**2)
    def log_like(arg):
        z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,A_Vdb,A_Vdc,Mba,Mca,Mda = arg
        return -0.5 * chi2sumall(arg)
    def log_posterior(arg):
        if not np.isfinite(log_prior(arg)):
            return -np.inf
        return log_prior(arg) + log_like(arg)

    ndim = 11 # number of parameters in the model
    nwalkers = 24 # number of MCMC walkers
    nsteps = 100000 # number of MCMC steps to take including burn-in
    # initial guesses in a small ball around the best-fit
    starting_guesses = np.c_[abs(np.random.normal(chimin_z,0.1*abs(chimin_z)+1e-05,nwalkers)),abs(np.random.normal(chimin_R_V,0.1*abs(chimin_R_V)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vba,0.1*abs(chimin_A_Vba)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vca,0.1*abs(chimin_A_Vca)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vda,0.1*abs(chimin_A_Vda)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vbc,0.1*abs(chimin_A_Vbc)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vdb,0.1*abs(chimin_A_Vdb)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vdc,0.1*abs(chimin_A_Vdc)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mba,0.1*abs(chimin_Mba)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mca,0.1*abs(chimin_Mca)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mda,0.1*abs(chimin_Mda)+1e-05,nwalkers))] # use this in case I want to run emcee instead of parallel tempering
    print("Running MCMC...")
    start_timeemcee8 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior,threads = 7)
    # run while showing the progress
    for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations* and I set ndim correctly
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Time for emcee with 7 threads: --- %s seconds ---" % (time.time() - start_timeemcee8))

    # plot the time laps
    import corner
    import pylab as plt
    from matplotlib.ticker import MaxNLocator
    plt.clf()
    fig, axes = plt.subplots(11, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].axhline(chimin_z, color="#888888", lw=2)
    axes[0].set_ylabel("1")
    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].axhline(chimin_R_V, color="#888888", lw=2)
    axes[1].set_ylabel("2")
    axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
    axes[2].yaxis.set_major_locator(MaxNLocator(5))
    axes[2].axhline(chimin_A_Vba, color="#888888", lw=2)
    axes[2].set_ylabel("3")
    axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
    axes[3].yaxis.set_major_locator(MaxNLocator(5))
    axes[3].axhline(chimin_A_Vca, color="#888888", lw=2)
    axes[3].set_ylabel("4")
    axes[4].plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)
    axes[4].yaxis.set_major_locator(MaxNLocator(5))
    axes[4].axhline(chimin_A_Vda, color="#888888", lw=2)
    axes[4].set_ylabel("5")
    axes[5].plot(sampler.chain[:, :, 5].T, color="k", alpha=0.4)
    axes[5].yaxis.set_major_locator(MaxNLocator(5))
    axes[5].axhline(chimin_A_Vbc, color="#888888", lw=2)
    axes[5].set_ylabel("6")
    axes[6].plot(sampler.chain[:, :, 6].T, color="k", alpha=0.4)
    axes[6].yaxis.set_major_locator(MaxNLocator(5))
    axes[6].axhline(chimin_A_Vdb, color="#888888", lw=2)
    axes[6].set_ylabel("7")
    axes[7].plot(sampler.chain[:, :, 7].T, color="k", alpha=0.4)
    axes[7].yaxis.set_major_locator(MaxNLocator(5))
    axes[7].axhline(chimin_A_Vdc, color="#888888", lw=2)
    axes[7].set_ylabel("8")
    axes[8].plot(sampler.chain[:, :, 8].T, color="k", alpha=0.4)
    axes[8].yaxis.set_major_locator(MaxNLocator(5))
    axes[8].axhline(chimin_Mba, color="#888888", lw=2)
    axes[8].set_ylabel("9")
    axes[9].plot(sampler.chain[:, :, 9].T, color="k", alpha=0.4)
    axes[9].yaxis.set_major_locator(MaxNLocator(5))
    axes[9].axhline(chimin_Mca, color="#888888", lw=2)
    axes[9].set_ylabel("10")
    axes[10].plot(sampler.chain[:, :, 10].T, color="k", alpha=0.4)
    axes[10].yaxis.set_major_locator(MaxNLocator(5))
    axes[10].axhline(chimin_Mda, color="#888888", lw=2)
    axes[10].set_ylabel("11")

    axes[2].set_xlabel("step number")
    fig.tight_layout(h_pad=0.0)
    fig.show()
    
    # print diagnostics and do corner plot
    nburn = nsteps/2 # "burn-in" to stabilize chains
    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
    #alpha_samp = sampler.flatchain.T[0]
    #beta_samp = sampler.flatchain.T[1]
    #sigma_samp = sampler.flatchain.T[2]
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
    print "median, std 1: ", np.median(samples[:,0]), np.std(samples[:,0])
    print "median, std 2: ", np.median(samples[:,1]), np.std(samples[:,1])
    print "median, std 3: ", np.median(samples[:,2]), np.std(samples[:,2])
    print "median, std 4: ", np.median(samples[:,3]), np.std(samples[:,3])
    print "median, std 5: ", np.median(samples[:,4]), np.std(samples[:,4])
    print "median, std 6: ", np.median(samples[:,5]), np.std(samples[:,5])
    print "median, std 7: ", np.median(samples[:,6]), np.std(samples[:,6])
    print "median, std 8: ", np.median(samples[:,7]), np.std(samples[:,7])
    print "median, std 9: ", np.median(samples[:,8]), np.std(samples[:,8])
    print "median, std 10: ", np.median(samples[:,9]), np.std(samples[:,9])
    print "median, std 11: ", np.median(samples[:,10]), np.std(samples[:,10])
    print chi2sumall([np.median(samples[:,0]),np.median(samples[:,1]),np.median(samples[:,2]),np.median(samples[:,3]),np.median(samples[:,4]),np.median(samples[:,5]),np.median(samples[:,6]),np.median(samples[:,7]),np.median(samples[:8]),np.median(samples[:,9]),np.median(samples[:,10])])[0]
    fig = corner.corner(samples, labels=["1", "2", "3","4", "5", "6","7", "8", "9", "10", "11"],truths=[chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vda,chimin_A_Vbc,chimin_A_Vdb,chimin_A_Vdc,chimin_Mba,chimin_Mca,chimin_Mda])
    fig.show()

allnodbdc = 0
if allnodbdc == 1:
    
    def log_prior(arg):
        z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,Mba,Mca,Mda = arg
        if z < 0 or z > 2.7 or R_V <= 0 or R_V >= 10 or A_Vba < 0 or A_Vca < 0 or A_Vda < 0 or A_Vbc < 0:
            return -np.inf # log(0)
        return 0 #-0.5 * ((R_V - 3.1)**2)/(0.5**2)
    def log_like(arg):
        z,R_V,A_Vba,A_Vca,A_Vda,A_Vbc,Mba,Mca,Mda = arg
        return -0.5 * chi2sumnodbdc(arg)
    def log_posterior(arg):
        if not np.isfinite(log_prior(arg)):
            return -np.inf
        return log_prior(arg) + log_like(arg)

    ndim = 9 # number of parameters in the model
    nwalkers = 20 # number of MCMC walkers
    nsteps = 100000 # number of MCMC steps to take including burn-in
    # initial guesses in a small ball around the best-fit
    if chimin_R_V > 10: chimin_R_V = 5
    starting_guesses = np.c_[abs(np.random.normal(chimin_z,0.1*abs(chimin_z)+1e-05,nwalkers)),abs(np.random.normal(chimin_R_V,0.1*abs(chimin_R_V)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vba,0.1*abs(chimin_A_Vba)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vca,0.1*abs(chimin_A_Vca)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vda,0.1*abs(chimin_A_Vda)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vbc,0.1*abs(chimin_A_Vbc)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mba,0.1*abs(chimin_Mba)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mca,0.1*abs(chimin_Mca)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mda,0.1*abs(chimin_Mda)+1e-05,nwalkers))] # use this in case I want to run emcee instead of parallel tempering
    print("Running MCMC...")
    start_timeemcee8 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior,threads = 7)
    # run while showing the progress
    for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations* and I set ndim correctly
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Time for emcee with 7 threads: --- %s seconds ---" % (time.time() - start_timeemcee8))

    # plot the time laps
    import corner
    import pylab as plt
    from matplotlib.ticker import MaxNLocator
    plt.clf()
    fig, axes = plt.subplots(9, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].axhline(chimin_z, color="#888888", lw=2)
    axes[0].set_ylabel("1")
    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].axhline(chimin_R_V, color="#888888", lw=2)
    axes[1].set_ylabel("2")
    axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
    axes[2].yaxis.set_major_locator(MaxNLocator(5))
    axes[2].axhline(chimin_A_Vba, color="#888888", lw=2)
    axes[2].set_ylabel("3")
    axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
    axes[3].yaxis.set_major_locator(MaxNLocator(5))
    axes[3].axhline(chimin_A_Vca, color="#888888", lw=2)
    axes[3].set_ylabel("4")
    axes[4].plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)
    axes[4].yaxis.set_major_locator(MaxNLocator(5))
    axes[4].axhline(chimin_A_Vda, color="#888888", lw=2)
    axes[4].set_ylabel("5")
    axes[5].plot(sampler.chain[:, :, 5].T, color="k", alpha=0.4)
    axes[5].yaxis.set_major_locator(MaxNLocator(5))
    axes[5].axhline(chimin_A_Vbc, color="#888888", lw=2)
    axes[5].set_ylabel("6")
    axes[6].plot(sampler.chain[:, :, 6].T, color="k", alpha=0.4)
    axes[6].yaxis.set_major_locator(MaxNLocator(5))
    axes[6].axhline(chimin_A_Vdb, color="#888888", lw=2)
    axes[6].set_ylabel("7")
    axes[7].plot(sampler.chain[:, :, 7].T, color="k", alpha=0.4)
    axes[7].yaxis.set_major_locator(MaxNLocator(5))
    axes[7].axhline(chimin_A_Vdc, color="#888888", lw=2)
    axes[7].set_ylabel("8")
    axes[8].plot(sampler.chain[:, :, 8].T, color="k", alpha=0.4)
    axes[8].yaxis.set_major_locator(MaxNLocator(5))
    axes[8].axhline(chimin_Mba, color="#888888", lw=2)
    axes[8].set_ylabel("9")

    axes[2].set_xlabel("step number")
    fig.tight_layout(h_pad=0.0)
    fig.show()
    
    # print diagnostics and do corner plot
    nburn = nsteps/2 # "burn-in" to stabilize chains
    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
    #alpha_samp = sampler.flatchain.T[0]
    #beta_samp = sampler.flatchain.T[1]
    #sigma_samp = sampler.flatchain.T[2]
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
    print "median, std 1: ", np.median(samples[:,0]), np.std(samples[:,0])
    print "median, std 2: ", np.median(samples[:,1]), np.std(samples[:,1])
    print "median, std 3: ", np.median(samples[:,2]), np.std(samples[:,2])
    print "median, std 4: ", np.median(samples[:,3]), np.std(samples[:,3])
    print "median, std 5: ", np.median(samples[:,4]), np.std(samples[:,4])
    print "median, std 6: ", np.median(samples[:,5]), np.std(samples[:,5])
    print "median, std 7: ", np.median(samples[:,6]), np.std(samples[:,6])
    print "median, std 8: ", np.median(samples[:,7]), np.std(samples[:,7])
    print "median, std 9: ", np.median(samples[:,8]), np.std(samples[:,8])
    print chi2sumnodbdc([np.median(samples[:,0]),np.median(samples[:,1]),np.median(samples[:,2]),np.median(samples[:,3]),np.median(samples[:,4]),np.median(samples[:,5]),np.median(samples[:,6]),np.median(samples[:,7]),np.median(samples[:8])])[0]
    fig = corner.corner(samples, labels=["1", "2", "3", "4", "5", "6", "7", "8", "9"],truths=[chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vda,chimin_A_Vbc,chimin_Mba,chimin_Mca,chimin_Mda])
    fig.show()


allnodadbdc = 0
if allnodadbdc == 1:
    
    def log_prior(arg):
        z,R_V,A_Vba,A_Vca,A_Vbc,Mba,Mca = arg
        if z < 0 or z > 2.7 or R_V <= 0 or R_V >= 10 or A_Vba < 0 or A_Vca < 0 or A_Vbc < 0:
            return -np.inf # log(0)
        return 0 #-0.5 * ((R_V - 3.1)**2)/(0.5**2)
    def log_like(arg):
        z,R_V,A_Vba,A_Vca,A_Vbc,Mba,Mca = arg
        return -0.5 * chi2sumnodadbdc(arg)
    def log_posterior(arg):
        if not np.isfinite(log_prior(arg)):
            return -np.inf
        return log_prior(arg) + log_like(arg)

    ndim = 7 # number of parameters in the model
    nwalkers = 20 # number of MCMC walkers
    nsteps = 100000 # number of MCMC steps to take including burn-in
    # initial guesses in a small ball around the best-fit
    if chimin_z > 2.7 or chimin_z <=0: chimin_z = 0.8
    if chimin_R_V > 10 or chimin_R_V <= 0: chimin_R_V = 3
    if chimin_A_Vba > 1 or chimin_A_Vba <= 0: chimin_A_Vba = 0.05
    if chimin_A_Vca > 1 or chimin_A_Vca <= 0: chimin_A_Vca = 0.05
    if chimin_A_Vbc > 1 or chimin_A_Vbc <= 0: chimin_A_Vbc = 0.05
    if chimin_Mba > 3 or chimin_Mba <= 0: chimin_Mba = 1
    if chimin_Mca > 3 or chimin_Mca <= 0: chimin_Mca = 1
    starting_guesses = np.c_[abs(np.random.normal(chimin_z,0.1*abs(chimin_z)+1e-05,nwalkers)),abs(np.random.normal(chimin_R_V,0.1*abs(chimin_R_V)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vba,0.1*abs(chimin_A_Vba)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vca,0.1*abs(chimin_A_Vca)+1e-05,nwalkers)),abs(np.random.normal(chimin_A_Vbc,0.1*abs(chimin_A_Vbc)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mba,0.1*abs(chimin_Mba)+1e-05,nwalkers)),abs(np.random.normal(chimin_Mca,0.1*abs(chimin_Mca)+1e-05,nwalkers))] # use this in case I want to run emcee instead of parallel tempering
    print("Running MCMC...")
    start_timeemcee8 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers,ndim,log_posterior,threads = 7)
    # run while showing the progress
    for i, result in enumerate(sampler.sample(starting_guesses,iterations=nsteps)): # fails unless I keep the keyword *iterations* and I set ndim correctly
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Time for emcee with 7 threads: --- %s seconds ---" % (time.time() - start_timeemcee8))

    # plot the time laps
    import corner
    import pylab as plt
    from matplotlib.ticker import MaxNLocator
    plt.clf()
    fig, axes = plt.subplots(7, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].axhline(chimin_z, color="#888888", lw=2)
    axes[0].set_ylabel("1")
    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].axhline(chimin_R_V, color="#888888", lw=2)
    axes[1].set_ylabel("2")
    axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
    axes[2].yaxis.set_major_locator(MaxNLocator(5))
    axes[2].axhline(chimin_A_Vba, color="#888888", lw=2)
    axes[2].set_ylabel("3")
    axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
    axes[3].yaxis.set_major_locator(MaxNLocator(5))
    axes[3].axhline(chimin_A_Vca, color="#888888", lw=2)
    axes[3].set_ylabel("4")
    axes[4].plot(sampler.chain[:, :, 4].T, color="k", alpha=0.4)
    axes[4].yaxis.set_major_locator(MaxNLocator(5))
    axes[4].axhline(chimin_A_Vda, color="#888888", lw=2)
    axes[4].set_ylabel("5")
    axes[5].plot(sampler.chain[:, :, 5].T, color="k", alpha=0.4)
    axes[5].yaxis.set_major_locator(MaxNLocator(5))
    axes[5].axhline(chimin_A_Vbc, color="#888888", lw=2)
    axes[5].set_ylabel("6")
    axes[6].plot(sampler.chain[:, :, 6].T, color="k", alpha=0.4)
    axes[6].yaxis.set_major_locator(MaxNLocator(5))
    axes[6].axhline(chimin_A_Vdb, color="#888888", lw=2)
    axes[6].set_ylabel("7")

    axes[2].set_xlabel("step number")
    fig.tight_layout(h_pad=0.0)
    fig.show()
    
    # print diagnostics and do corner plot
    nburn = nsteps/2 # "burn-in" to stabilize chains
    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim)) # combines all walkers, without burn-in
    #alpha_samp = sampler.flatchain.T[0]
    #beta_samp = sampler.flatchain.T[1]
    #sigma_samp = sampler.flatchain.T[2]
    #print("Autocorrelation time:", sampler.get_autocorr_time())
    print "acceptance fraction: ", np.median(sampler.acceptance_fraction)
    print "median, std 1: ", np.median(samples[:,0]), np.std(samples[:,0])
    print "median, std 2: ", np.median(samples[:,1]), np.std(samples[:,1])
    print "median, std 3: ", np.median(samples[:,2]), np.std(samples[:,2])
    print "median, std 4: ", np.median(samples[:,3]), np.std(samples[:,3])
    print "median, std 5: ", np.median(samples[:,4]), np.std(samples[:,4])
    print "median, std 6: ", np.median(samples[:,5]), np.std(samples[:,5])
    print "median, std 7: ", np.median(samples[:,6]), np.std(samples[:,6])
    print chi2sumnodadbdc([np.median(samples[:,0]),np.median(samples[:,1]),np.median(samples[:,2]),np.median(samples[:,3]),np.median(samples[:,4]),np.median(samples[:,5]),np.median(samples[:,6])])[0]
    fig = corner.corner(samples, labels=["1", "2", "3", "4", "5", "6", "7"],truths=[chimin_z,chimin_R_V,chimin_A_Vba,chimin_A_Vca,chimin_A_Vbc,chimin_Mba,chimin_Mca])
    fig.show()



