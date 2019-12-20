# Cone search
import numpy as np
import astropy.units as u

search = "Gaia"
#search = "PS1"
#file = "/Volumes/LaCieSubaru/Ciprian_candidates/candidatesJan2019.cat"
#file = "/Users/cerusu/OneDrive - Subaru Telescope/candidates.cat"
file = "/Users/cerusu/Dropbox/clustermostpromising.txt"
cat = np.loadtxt(file, usecols = [0,1], unpack=True)

if search == "Gaia":
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    #rad = 9.0 # arcsec # used for Gaia
    rad = 20.0 # arcsec
    radius = u.Quantity(rad, u.arcsec)
    #for i in range(1):
    for i in range(len(cat[0])):
        coord = SkyCoord(ra=cat[0][i], dec=cat[1][i], unit=(u.degree, u.degree), frame='icrs')
        j = Gaia.cone_search_async(coord, radius)
        r = j.get_results()
        #print r.columns
        print cat[0][i], cat[1][i]
        print 'dist           ra             dec            G       astrom_exc_noise astrom_exc_noise_sig proper_motion_significance'
        for j in range(len(r)):
            #if j ==1: print cat[0][i], cat[1][i], SkyCoord(cat[0][i], cat[1][i], unit='deg').separation(SkyCoord(r[j]['ra'], r[j]['dec'], unit='deg')).arcsec, r[j]['phot_g_mean_mag'],r[j]['phot_bp_mean_mag'],r[j]['phot_rp_mean_mag']
            print r[j]['dist']*3600,r[j]['ra'],r[j]['dec'],r[j]['phot_g_mean_mag'],r[j]['astrometric_excess_noise'],r[j]['astrometric_excess_noise_sig'],np.sqrt((r[j]['pmra']/r[j]['pmra_error'])**2+(r[j]['pmdec']/r[j]['pmdec_error'])**2)

def panstarrs_query(ra_deg, dec_deg, rad_asec, maxmag=23,
                    maxsources=1): # from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['RAJ2000', 'DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag)},
                    row_limit=maxsources)
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % (rad_asec / 3600.0)),
                               catalog="II/349/ps1")[0]

if search == "PS1":
    rad = 1.0 # arcsec
    from astroquery.vizier import Vizier
    import astropy.coordinates as coord
    for i in range(len(cat[0])):
        try: print(panstarrs_query(cat[0][i], cat[1][i], 2))
        except: pass
