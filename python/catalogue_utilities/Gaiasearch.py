# read RA and DEC from a file and output a corresponding DS9 regions file
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

rad = 9.0 # arcsec
#file = "/Volumes/LaCieSubaru/Ciprian_candidates/candidatesJan2019.cat"
file = "/Users/cerusu/OneDrive - Subaru Telescope/Ciprian/candidatesJan2019.cat"
radius = u.Quantity(rad, u.arcsec)
cat = np.loadtxt(file, usecols = [0,1], unpack=True)
#for i in range(1):
for i in range(len(cat[0])):
    coord = SkyCoord(ra=cat[0][i], dec=cat[1][i], unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()
    print cat[0][i], cat[1][i]
    print 'dist           ra             dec            G       astrom_exc_noise parallax/err pmra/err    pmdec/err'
    for j in range(len(r)):
        print r[j]['dist']*3600,r[j]['ra'],r[j]['dec'],r[j]['phot_g_mean_mag'],r[j]['astrometric_excess_noise'],abs(r[j]['parallax_over_error']),np.sqrt((r[j]['pmra']/r[j]['pmra_error'])**2+(r[j]['pmdec']/r[j]['pmdec_error'])**2)
