# Run this code in order to sample possible galaxies part of an incomplete spectroscopic group
# no uncertainties assumed on group centroid because they are small enough to ignore
# if I get error ValueError: Cannot take a larger sample than population when 'replace=False' simply run the code again or increase the faintmagspec or photoztolerance

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

selectpdz = False # whether or not to actually select a number of desired samples, not just compute the theoretical distributions
if selectpdz == True: samples = 10 # define any number of desired samples
#mode = "poisson"
mode = "mcmc"
photoztolerance = 1.5 # number of sigmas; doesn't apply to PG1115
#lens = 'PG1115'
lens = 'WFI2033'
#zgroup = 0.6588 # WFI2033lens
zgroup = 0.4956 # WFI2033
if lens == 'PG1115': zgroup = 0.3097
if lens == 'WFI2033': limmag = 22.5
if lens == 'PG1115': limmag = 22.5 # slightly fainter than faintest in Momcheva
if lens == 'WFI2033': faintmagspec = 22.5
if lens == 'PG1115': faintmagspec = 22.5 # PG1115, slightly fainter than faintest in Momcheva
if lens == 'WFI2033': center_lensx = '20:33:42.080'; center_lensy = '-47:23:43.00' # WFI2033
if lens == 'PG1115': center_lensx = '11:18:16.90'; center_lensy = '+07:45:59.00' # PG1115
center_lens = SkyCoord(center_lensx + ' ' + center_lensy, frame='fk5', unit=(u.hourangle, u.deg))
if zgroup == 0.6588:
    center_groupx = '308.43557011'; center_groupy = '-47.37411275'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 26 # converted to arcsec
    virrad = 142 # actually R_200 in arcsec - from Dominique's email on Dec 3 2018
    virrad_err = 29
if zgroup == 0.4956:
    center_groupx = '308.46337200'; center_groupy = '-47.36336725'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 60
    virrad = 164
    virrad_err = 33
if zgroup == 0.3097:
    center_groupx = '169.5681'; center_groupy = '7.7648'
    center_group = SkyCoord(center_groupx + ' ' + center_groupy, frame='fk5', unit=(u.deg, u.deg))
    err_group = 20
    virrad = 180
    virrad_err = 36 # arbitrary 20%, as average between the ones from Dominique

sep_groupx = center_lens.separation(SkyCoord(center_groupx + ' ' + center_lensy, frame='fk5', unit=(u.deg, u.deg))).arcsec
sep_groupy = center_lens.separation(SkyCoord(center_lensx + ' ' + center_groupy, frame='fk5', unit=(u.hourangle, u.deg))).arcsec
# arcsec # using R_200, which is robust against Dominique's and is used in the reference papers
# Each of the 2 groups in WFI2033 have 3 galaxies outside R_200 so I should not count these
if zgroup == 0.6588: observed_members = 16 # inside virial radius; the last column in the table from Dominique
if zgroup == 0.4956: observed_members = 7
if zgroup == 0.3097: observed_members = 13
if lens == 'PG1115':
    file = '/Users/cerusu/Dropbox/Davis_work/code/PG1115/PG1115.cat'
    data = np.loadtxt(file,usecols=[2,3,4,8,11])
    othergals = np.loadtxt('/Users/cerusu/Dropbox/Davis_work/code/PG1115/galotherredshift.cat', usecols=[0])
    groupgals = np.loadtxt('/Users/cerusu/Dropbox/Davis_work/code/PG1115/galgroup.cat')
    ra = 0
    dec = 1
    id = 2
    flux_rad = 3
    r = 4
if lens == 'WFI2033':
    file = "/Users/cerusu/Dropbox/Davis_work/code/WFI2033/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23_withbpzeazylephareclassified_IRACmagslephareclassifiedF160W.cat" # WFI2033
    data = np.loadtxt(file,usecols=[2,3,4,8,28,29,30,40,97])
ra = 0
dec = 1
i = 2
id = 3
z = 4
zinf = 5
zsup = 6
spec = 7
cls = 8

coord = SkyCoord(ra=data[:,ra]*u.degree, dec=data[:,dec]*u.degree, frame='fk5')
sep = coord.separation(center_lens).arcsec
if lens == 'WFI2033': all = len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,cls] >= 0)])
if lens == 'PG1115': all = len(data[:,id][(sep <= 120) & (data[:,r] <= limmag) & (data[:,flux_rad] >= 1.25)])
print "gals: ",all
if lens == 'WFI2033': specs120 = len(data[:,id][(sep <= 120) & (data[:,i] <= limmag) & (data[:,spec] > 0)])
if lens == 'PG1115':
    observed120_membersID = groupgals[:,0]
    specs120 = len(observed120_membersID) + 1 + np.shape(othergals)[0]  # adding the lens to the counts
print "specs inside 120\": ",specs120
if lens == 'WFI2033': observed120_membersID = data[:,id][(data[:,spec] <= zgroup + 0.01) & (data[:,spec] >= zgroup - 0.01) & (sep <= 120) & (data[:,i] <= limmag) & (data[:,cls] >= 0)]
if (lens == 'PG1115') or (zgroup == 0.6588): print 'members inside 120\":',len(observed120_membersID) + 1
else: print 'members inside 120\":',len(observed120_membersID)
if lens == 'WFI2033': pool = data[:,id][(data[:,spec] == -1) & (sep <= 120) & (data[:,cls] >= 0) & (data[:,z] - photoztolerance * (data[:,z] - data[:,zinf]) <= zgroup) & (data[:,z] + photoztolerance * (data[:,zsup] - data[:,z]) >= zgroup) & (data[:,i] <= faintmagspec)]

if lens == 'PG1115':
    pool = data[:,id][(sep <= 120) & (data[:,r] <= limmag) & (data[:,flux_rad] >= 1.25)]
    for i in range(len(othergals)):
        pool = np.delete(pool,np.where(pool == othergals[i]))
    for i in range(np.shape(groupgals)[0]):
        pool = np.delete(pool,np.where(pool == groupgals[:,0][i]))
print 'size of pool: ',len(pool)

#veldisp0.66 = 502+/-83 -> Fig 5 Andreon 2010 [median range: 10-(20)-32; 68% range around central '()' median 13-30, in quadrature 20-12+16] galaxies inside R_200
#veldisp0.49 = 518+/-99 -> Fig 5 Andreon 2010 [median range: 10-(20)-39; 68% range around central '()' median 13-30, in quadrature 20-12+21] galaxies inside R_200
#veldisp0.31 = 390+50-60 -> Fig 5 Andreon 2010 [median range: 5-(8)-12; 68% range around central '()' median 6-12, in quadrature 8-4+6] galaxies inside R_200
#veldisp0.31 = 440+90-80 -> Fig 5 Andreon 2010 [median range: 6-(13)-20; 68% range around central '()' median 8-18, in quadrature 13+/-8.5] galaxies inside R_200
if zgroup == 0.6588: med = 20; stdinf = 12; stdsup = 16
if zgroup == 0.4956: med = 20; stdinf = 12; stdsup = 21
if zgroup == 0.3097: med = 8; stdinf = 4; stdsup = 6 # Wilson
#if zgroup == 0.3097: med = 13; stdinf = 8.5; stdsup = 8.5 # Momcheva

if mode == "mcmc":
    def expected_members():
        while True:
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
            frac = 1.0*len(x[(x**2 + y**2 + z**2 <= 1) & ((x-cx)**2 + (y-cy)**2 <=rad**2)]) / len(x[x**2 + y**2 + z**2 <= 1])
            nr = 0
            while True:
                nr += 1
                rand = np.random.uniform(0,1,1)[0]
                if rand <= 0.5: x = int(round(med - np.abs(np.random.normal(med, stdinf, 1) - med))) # based on the velocity dispersion - concentration relation above
                else: x = int(round(med + np.abs(np.random.normal(med, stdsup, 1) - med)))
                #print x,frac * x, observed_members, len(observed120_membersID)
                if (lens == 'PG1115') or (zgroup == 0.6588):
                    if ((x >= observed_members) and (frac * x >= len(observed120_membersID) + 1)) or (nr == 100): break
                else:
                    if ((x >= observed_members) and (frac * x >= len(observed120_membersID))) or (nr == 100): break
            if (lens == 'PG1115') or (zgroup == 0.6588):
                if ((x >= observed_members) and (frac * x >= len(observed120_membersID) + 1)): break
            else:
                if ((x >= observed_members) and (frac * x >= len(observed120_membersID))): break
        return x,frac
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
        fracnoprior = 1.0*len(x[(x**2 + y**2 + z**2 <= 1) & ((x-cx)**2 + (y-cy)**2 <=rad**2)]) / len(x[x**2 + y**2 + z**2 <= 1])
        while True:
            rand = np.random.uniform(0,1,1)[0]
            if rand <= 0.5: x = int(round(med - np.abs(np.random.normal(med, stdinf, 1) - med))) # based on the velocity dispersion - concentration relation above
            else: x = int(round(med + np.abs(np.random.normal(med, stdsup, 1) - med)))
            break
        return x,fracnoprior
if mode == "poisson":
    def expected_members():
        while True:
            #x = np.random.poisson(1.0 * all * len(observed120_membersID)/specs120, 1)
            if (lens == 'PG1115') or (zgroup == 0.6588):
                x = np.random.poisson(1.0 * all * (len(observed120_membersID) + 1)/specs120, 1)  # adding the lens to the counts
                if x >= len(observed120_membersID)+1: break
            else:
                x = np.random.poisson(1.0 * all * (len(observed120_membersID))/specs120, 1)  # adding the lens to the counts
                if x >= len(observed120_membersID): break
        return x
if mode == "poisson":
    def expected_members_noprior():
        while True:
            #x = np.random.poisson(1.0 * all * len(observed120_membersID)/specs, 1)
            if (lens == 'PG1115') or (zgroup == 0.6588): x = np.random.poisson(1.0 * all * (len(observed120_membersID) + 1)/specs120, 1)
            else: x = np.random.poisson(1.0 * all * len(observed120_membersID)/specs120, 1)
            break
        return x

pdz = np.array([])
theorysamples = 50000
for i in range(theorysamples):
    if i % 10 == 0: print i,'/',theorysamples
    #print i
    expected = expected_members()
    if mode == "mcmc": pdz = np.append(pdz,int(round(expected[1] * expected[0])))
    if mode == "poisson": pdz = np.append(pdz,expected[0])
pdznoprior = np.array([])
for i in range(theorysamples):
    if i % 10 == 0: print i,'/',theorysamples
    #print i
    expected = expected_members_noprior()
    if mode == "mcmc": pdznoprior = np.append(pdznoprior,int(round(expected[1] * expected[0])))
    if mode == "poisson": pdznoprior = np.append(pdznoprior,expected[0])

if selectpdz == True:
    pdzselect = np.array([])
    for i in range(samples):
        expected = expected_members()
        if mode == "mcmc":
            if (lens == 'PG1115') or (zgroup == 0.6588): missing120_membersID = np.random.choice(a=pool, size=int(round(expected[1] * expected[0] - len(observed120_membersID) - 1)), replace=False)
            else: missing120_membersID = np.random.choice(a=pool, size=int(round(expected[1] * expected[0] - len(observed120_membersID))), replace=False)
        if mode == "poisson":
            if (lens == 'PG1115') or (zgroup == 0.6588):  missing120_membersID = np.random.choice(a=pool, size=int(round(expected[0] - len(observed120_membersID) - 1)), replace=False)
            else: missing120_membersID = np.random.choice(a=pool, size=int(round(expected[0] - len(observed120_membersID))), replace=False)
        missing120_membersra = np.array([])
        missing120_membersdec = np.array([])
        #pdzselect = np.append(pdzselect,len(missing120_membersID)+len(observed120_membersID))
        if (lens == 'PG1115') or (zgroup == 0.6588): pdzselect = np.append(pdzselect,len(missing120_membersID)+len(observed120_membersID) + 1)
        else: pdzselect = np.append(pdzselect,len(missing120_membersID)+len(observed120_membersID))
        for j in range(len(missing120_membersID)):
            missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==missing120_membersID[j]][0])
            missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==missing120_membersID[j]][0])
        for k in range(len(observed120_membersID)):
            # removing also the known group members
            missing120_membersra = np.append(missing120_membersra,data[:,ra][data[:,id]==observed120_membersID[k]][0])
            missing120_membersdec = np.append(missing120_membersdec,data[:,dec][data[:,id]==observed120_membersID[k]][0])
        #if zgroup == 0.6588: np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgrouphandpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
        #if zgroup == 0.4956: np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/WFI2033/removelensgroup049handpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')
        #if zgroup == 0.3097: np.savetxt("/Users/cerusu/Dropbox/Davis_work/code/PG1115/removelensgroup"+mode+"handpicked"+str(i)+".cat",np.c_[missing120_membersra,missing120_membersdec],fmt='%.8f %.9f')

import pylab as plt
plt.clf()
if (lens == 'PG1115') or (zgroup == 0.6588):
    plt.hist(pdz - len(observed120_membersID) - 1,bins=50,normed=True,label='w/ observed number prior')
    plt.hist(pdznoprior - len(observed120_membersID) - 1,bins=50,normed=True,label='w/o observed number prior',alpha = 0.5)
    if selectpdz == True: plt.hist(pdzselect - len(observed120_membersID) - 1,bins=50,normed=True,alpha = 0.5)
else:
    plt.hist(pdz - len(observed120_membersID),bins=50,normed=True,label='w/ observed number prior')
    plt.hist(pdznoprior - len(observed120_membersID),bins=50,normed=True,label='w/o observed number prior',alpha = 0.5)
    if selectpdz == True: plt.hist(pdzselect - len(observed120_membersID),bins=50,normed=True,alpha = 0.5)
plt.xlabel(r'Expected number of missing members', fontsize=20)
plt.ylabel(r'normalized counts', fontsize=20)
plt.legend(loc="upper left")
plt.show()
if (lens == 'PG1115') or (zgroup == 0.6588): print np.percentile(pdz - len(observed120_membersID) - 1,[16,50,84]) # I ran the code several times untill the two distributions match fairly well
else: print np.percentile(pdz - len(observed120_membersID),[16,50,84])

# comment out the following lines. In case I want to produce a paper plot, I need to run the following lines in iPython so I can run for mcmc and poisson in turn.
# first run with poisson
#pdz_poisson_lens = pdz
#pdznoprior_poisson_lens = pdznoprior
#len_lens = len(observed120_membersID)
# next run with mcmc
#pdz_mcmc_los = pdz
#pdznoprior_mcmc_los = pdznoprior
#len_los = len(observed120_membersID)
#plt.clf()
#bin = 50
#if (lens == 'PG1115'):
#    plt.hist(pdz_mcmc - len(observed120_membersID) - 1,bins=bin,normed=True,label='Volume-based; w/ observed number prior',color='r')
#    plt.hist(pdznoprior_mcmc - len(observed120_membersID) - 1,bins=bin,normed=True,label='Volume-based; w/o observed number prior',alpha = 0.5,color='r')
#    plt.hist(pdznoprior_poisson - len(observed120_membersID) - 1,bins=bin,normed=True,label='Poisson-based; w/o observed number prior',alpha = 0.5,color='k')
#else:
#    plt.hist(pdz_mcmc - len(observed120_membersID),bins=bin,normed=True,label='Volume-based; w/ observed number prior',color='r')
#    plt.hist(pdznoprior_mcmc - len(observed120_membersID),bins=bin,normed=True,label='Volume-based; w/o observed number prior',alpha = 0.5,color='r')
#    plt.hist(pdznoprior_poisson - len(observed120_membersID),bins=bin,normed=True,label='Poisson-based; w/o observed number prior',alpha = 0.5,color='k')
#plt.xlabel(r'Expected number of missing group members', fontsize=16)
#plt.ylabel(r'normalized counts', fontsize=16)
#plt.legend(loc="upper left")
#plt.title(r'PG1115+080', fontsize=20)
#if (lens == 'PG1115'): print np.percentile(pdz - len(observed120_membersID) - 1,[16,50,84]) # I ran the code several times untill the two distributions match fairly well
#else: print np.percentile(pdz - len(observed120_membersID),[16,50,84])
#plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/PG1115/estimatinggroupmembersPG1115.png', dpi=250, bbox_inches='tight')

#plt.clf()
#bin = 100
#plt.hist(pdz_mcmc_lens - len_lens - 1,bins=bin,normed=True,label='z=0.66 Volume-based; w/ observed number prior',color='k',linestyle='--',histtype='step')
#plt.hist(pdznoprior_mcmc_lens - len_lens - 1,bins=bin,normed=True,label='z=0.66 Volume-based; w/o observed number prior',color='k',linestyle=':',histtype='step')
#plt.hist(pdznoprior_poisson_lens - len_lens - 1,bins=bin,normed=True,label='z=0.66 Poisson-based',color='k',linestyle='-',histtype='step')
#plt.hist(pdz_mcmc_los - len_los,bins=bin,normed=True,label='z=0.49 Volume-based; w/ observed number prior',color='r',linestyle='--',histtype='step')
#plt.hist(pdznoprior_mcmc_los - len_los,bins=bin,normed=True,label='z=0.49 Volume-based; w/o observed number prior',color='r',linestyle=':',histtype='step')
#plt.hist(pdznoprior_poisson_los - len_los,bins=bin,normed=True,label='z=0.49 Poisson-based',color='r',linestyle='-',histtype='step')
#plt.xlabel(r'Expected number of missing group members', fontsize=16)
#plt.ylabel(r'normalized counts', fontsize=16)
#plt.legend(loc="upper left")
#plt.xlim([-20,30])
#plt.savefig('/Users/cerusu/Dropbox/Davis_work/code/WFI2033/estimatingmissinggroupmembersWFI2033.png', dpi=250, bbox_inches='tight')
