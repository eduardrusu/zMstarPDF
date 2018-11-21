# Run this code in order to sample possible galaxies part of an incomplete spectroscopic group
# no uncertainties assumed on group centroid because they are small enough to ignore
# if I get error ValueError: Cannot take a larger sample than population when 'replace=False' simply run the code again or increase the faintmagspec or photoztolerance

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
#mode = "poisson"
mode = "mcmc"
samples = 20
photoztolerance = 1.5 # number of sigmas
#zgroup = 0.6588 # WFI2033lens
#zgroup = 0.4956 # WFI2033
zgroup = 0.3097 # PG1115
limmag = 22.5 # WFI2033
limmag = 22.5 # PG1115, slightly fainter than faintest in Momcheva
faintmagspec = 22.5 # WFI2033
faintmagspec = 22.5 # PG1115, slightly fainter than faintest in Momcheva
#center_lensx = '20:33:42.080'; center_lensy = '-47:23:43.00' # WFI2033
center_lensx = '11:18:16.90'; center_lensy = '+07:45:59.00' # PG1115
center_lens = SkyCoord(center_lensx + ' ' + center_lensy, frame='fk5', unit=(u.hourangle, u.deg))
if zgroup == 0.6588:
    center_groupx = '308.43557011'; center_groupy = '-47.37411275'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 60 # converted to arcsec
    virrad = 130 # actually R_200 in arcsec
    virrad_err = 30
if zgroup == 0.4956:
    center_groupx = '308.46337200'; center_groupy = '-47.36336725'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 26
    virrad = 170
    virrad_err = 30
if zgroup == 0.3097: # PG1115
    center_groupx = '169.5681'; center_groupy = '7.7648'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 20
    virrad = 180
    virrad_err = 18 # arbitrary 20%, as average between the ones from Dominique

sep_groupx = center_lens.separation(SkyCoord(center_groupx + ' ' + center_lensy, frame='fk5', unit=(u.deg, u.deg))).arcsec
sep_groupy = center_lens.separation(SkyCoord(center_lensx + ' ' + center_groupy, frame='fk5', unit=(u.hourangle, u.deg))).arcsec
# arcsec # using R_200, which is robust agains Dominique's and is used in the reference papers
# Each of the 2 groups have 3 galaxies outside R_200 so I should not count these
if zgroup == 0.6588: observed_members = 19 # inside virial radius
if zgroup == 0.4956: observed_members = 10
if zgroup == 0.3097: observed_members = 13
#file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat" # WFI2033
if zgroup == 0.3097: file = '/Users/cerusu/Dropbox/Davis_work/code/PG1115/PG1115.cat'
#data = np.loadtxt(file,usecols=[2,3,4,8,28,29,30,40,97])
if zgroup == 0.3097:
    data = np.loadtxt(file,usecols=[2,3,4,8,11])
ra = 0
dec = 1
i = 2
id = 3
z = 4
zinf = 5
zsup = 6
spec = 7
cls = 8
if zgroup == 0.3097:
    ra = 0
    dec = 1
    id = 2
    flux_rad = 3
    i = 4

coord = SkyCoord(ra=data[:,ra]*u.degree, dec=data[:,dec]*u.degree, frame='fk5')
sep = coord.separation(center_lens).arcsec
#all = len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,cls] >= 0)])
all = len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,flux_rad] >= 1.25)])
print "gals: ",all
#specs = len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,spec] > 0)])
specs = 11 # PG1115 mcmc
#specs = 33 # PG1115 poisson
print "specs: ",specs
#observed120_membersID = data[:,id][(data[:,spec] <= zgroup + 0.01) & (data[:,spec] >= zgroup - 0.01) & (sep <= 120) & (data[:,i] <= limmag) & (data[:,cls] >= 0)]
#print 'members',len(observed120_membersID)
#pool = data[:,id][(data[:,spec] == -1) & (sep <= 120) & (data[:,cls] >= 0) & (data[:,z] - photoztolerance * (data[:,z] - data[:,zinf]) <= zgroup) & (data[:,z] + photoztolerance * (data[:,zsup] - data[:,z]) >= zgroup) & (data[:,i] <= faintmagspec)]
#print 'pool',len(pool)

if mode == "mcmc":
    def expected_members():
        # sample uniformly a cube of size 2 x unit radius
        x = np.random.uniform(-1,1,10000)
        y = np.random.uniform(-1,1,10000)
        z = np.random.uniform(-1,1,10000)
        sep_group = np.sqrt(np.random.normal(sep_groupx, err_group, 1)**2 + np.random.normal(sep_groupy, err_group, 1)**2)
        virrad_sample = np.max([10,sep_group,np.abs(np.random.normal(virrad, virrad_err, 1))]) # using lower limit to avoid numerical problems
        cx = 1.0 * sep_group / virrad_sample
        cy = 0
        rad = 120.0 / virrad_sample
        # fraction of volume of the virial sphere which contains the 120"-radius cylinder centered on the lens
        global frac
        frac = 1.0*len(x[(x**2 + y**2 + z**2 <= 1) & ((x-cx)**2 + (y-cy)**2 <=rad**2)]) / len(x[x**2 + y**2 + z**2 <= 1])
        while True:
            #veldisp0.66 = 502+/-83 -> Fig 5 Andreon 2010 [median range: 10-(20)-32; 68% range around central '()' median 13-30, in quadrature 20-12+16] galaxies inside R_200
            #veldisp0.49 = 518+/-99 -> Fig 5 Andreon 2010 [median range: 10-(20)-39; 68% range around central '()' median 13-30, in quadrature 20-12+21] galaxies inside R_200
            #veldisp0.31 = 390+50-60 -> Fig 5 Andreon 2010 [median range: 5-(8)-12; 68% range around central '()' median 6-12, in quadrature 8-4+6] galaxies inside R_200

            if zgroup == 0.6588: med = 20; stdinf = 12; stdsup = 16
            if zgroup == 0.4956: med = 20; stdinf = 12; stdsup = 21
            if zgroup == 0.3097: med = 8; stdinf = 4; stdsup = 6
            if zgroup == 0.3097: med = 13; stdinf = 8.5; stdsup = 8.5
            rand = np.random.uniform(0,1,1)[0]
            if rand <= 0.5: x = med - np.abs(np.random.normal(med, stdinf, 1).astype(int)[0] - med) # based on the velocity dispersion - concentration relation above
            else: x = med + np.abs(np.random.normal(med, stdsup, 1).astype(int)[0] - med)
            #print x,frac * x, observed_members, len(observed120_membersID)
            #if (x > observed_members) and (frac * x > len(observed120_membersID)): break
            if (x > observed_members) and (frac * x > specs): break # PG1115
        return x
if mode == "mcmc":
    def expected_members_noprior():
        # sample uniformly a cube of size 2 x unit radius
        x = np.random.uniform(-1,1,10000)
        y = np.random.uniform(-1,1,10000)
        z = np.random.uniform(-1,1,10000)
        sep_group = np.sqrt(np.random.normal(sep_groupx, err_group, 1)**2 + np.random.normal(sep_groupy, err_group, 1)**2)
        virrad_sample = np.max([10,sep_group,np.abs(np.random.normal(virrad, virrad_err, 1))]) # using lower limit to avoid numerical problems
        cx = 1.0 * sep_group / virrad_sample
        cy = 0
        rad = 120.0 / virrad_sample
        # fraction of volume of the virial sphere which contains the 120"-radius cylinder centered on the lens
        global frac
        frac = 1.0*len(x[(x**2 + y**2 + z**2 <= 1) & ((x-cx)**2 + (y-cy)**2 <=rad**2)]) / len(x[x**2 + y**2 + z**2 <= 1])
        while True:
            #veldisp0.66 = 502+/-83 -> Fig 5 Andreon 2010 [median range: 10-(20)-32; 68% range around central '()' median 13-30, in quadrature 20-12+16] galaxies inside R_200
            #veldisp0.49 = 518+/-99 -> Fig 5 Andreon 2010 [median range: 10-(20)-39; 68% range around central '()' median 13-30, in quadrature 20-12+21] galaxies inside R_200
            #veldisp0.31 = 390+50-60 -> Fig 5 Andreon 2010 [median range: 5-(8)-12; 68% range around central '()' median 6-12, in quadrature 8-4+6] galaxies inside R_200
            #veldisp0.31 = 440+90-80 -> Fig 5 Andreon 2010 [median range: 6-(13)-20; 68% range around central '()' median 8-18, in quadrature 13+/-8.5] galaxies inside R_200
            if zgroup == 0.6588: med = 20; stdinf = 12; stdsup = 16
            if zgroup == 0.4956: med = 20; stdinf = 12; stdsup = 21
            if zgroup == 0.3097: med = 8; stdinf = 4; stdsup = 6
            if zgroup == 0.3097: med = 13; stdinf = 8.5; stdsup = 8.5
            rand = np.random.uniform(0,1,1)[0]
            if rand <= 0.5: x = med - np.abs(np.random.normal(med, stdinf, 1).astype(int)[0] - med) # based on the velocity dispersion - concentration relation above
            else: x = med + np.abs(np.random.normal(med, stdsup, 1).astype(int)[0] - med)
            break
        return x
if mode == "poisson":
    def expected_members():
        while True:
            #x = np.random.poisson(1.0 * all * len(observed120_membersID)/specs, 1)
            x = np.random.poisson(1.0 * all * 11/specs, 1) #PG1115
            #if x >= len(observed120_membersID): break
            if x >= 11: break # PG1115
        return x
if mode == "poisson":
    def expected_members_noprior():
        while True:
            #x = np.random.poisson(1.0 * all * len(observed120_membersID)/specs, 1)
            x = np.random.poisson(1.0 * all * 11/specs, 1) #PG1115
            break
        return x

pdz = np.array([])
for i in range(10000):
    if i % 10 == 0: print i,'/10000'
    #print i
    expected = expected_members()
    if mode == "mcmc": pdz = np.append(pdz,frac * expected)
    if mode == "poisson": pdz = np.append(pdz,expected)
pdznoprior = np.array([])
for i in range(10000):
    if i % 10 == 0: print i,'/10000'
    #print i
    expected = expected_members_noprior()
    if mode == "mcmc": pdznoprior = np.append(pdznoprior,frac * expected)
    if mode == "poisson": pdznoprior = np.append(pdznoprior,expected)

pdzselect = np.array([])
for i in range(samples):
    if mode == "mcmc": missing120_membersID = np.random.choice(a=pool, size=int(frac * expected_members() - len(observed120_membersID)), replace=False)
    if mode == "poisson": missing120_membersID = np.random.choice(a=pool, size=int(expected_members() - len(observed120_membersID)), replace=False)
    missing120_membersra = np.array([])
    missing120_membersdec = np.array([])
    pdzselect = np.append(pdzselect,len(missing120_membersID)+len(observed120_membersID))
    for j in range(len(missing120_membersID)):
        missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==missing120_membersID[j]][0])
        missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==missing120_membersID[j]][0])
    for k in range(len(observed120_membersID)):
        # removing also the known group members
        missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==observed120_membersID[k]][0])
        missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==observed120_membersID[k]][0])
    #if zgroup == 0.6588: np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgrouphandpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
    #if zgroup == 0.4956: np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgroup049handpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
import pylab as plt
plt.clf()
plt.hist(pdz-11,bins=50,normed=True,label='w/ observed number prior')
plt.hist(pdznoprior-11,bins=50,normed=True,label='w/o observed number prior')
#plt.hist(pdzselect,normed=True,alpha = 0.5)
plt.xlabel(r'Expected number of missing members', fontsize=20)
plt.ylabel(r'normalized counts', fontsize=20)
plt.legend(loc="upper left")
plt.show()
print np.percentile(pdz-11,[16,50,84]) # I ran the code several times untill the two distributions match fairly well
