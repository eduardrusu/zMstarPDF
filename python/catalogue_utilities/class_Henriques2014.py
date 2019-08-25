import numpy as np

class Henriques2014():
    """
    Returns the structure of the files of type GGL_los_8_0_0_N_4096_ang_4_Henriques2014_galaxies_on_plane_63_f.images.
    """
    galaxy_struct = np.dtype([
  ('galaxy_id'                 ,'i8'      ), #0_LL $       ,                   , id of galaxy (unique)
  ('halo_id'                   ,'i8'      ), #0_LL $       ,                   , id of (sub)halo, the galaxy belongs to(?)
  ('first_prog_gal'            ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('next_prog_gal'             ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('last_prog_gal'             ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('FOF_central_gal'           ,'i8'      ), #0_LL $       ,                   , id of fof halo the galaxy belong to (i.e. common id for all galaxies in same group or cluster)
  ('file_tree_nr'              ,'i8'      ), #0_LL $       ,                   , id of file containing the merger tree the galaxy belongs to
  ('descendant_gal'            ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('main_leaf_id'              ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('tree_root_id'              ,'i8'      ), #0_LL $       ,                   , id of some other galaxy needed for galaxy merger tree structure
  ('subhalo_id'                ,'i8'      ), #0_LL $       ,                   , id of (sub)halo, the galaxy belongs to(?)
  ('main_subhalo_id'           ,'i8'      ), #0_LL $       ,                   , id of main (sub)halo of fof halo, the galaxy belongs to(?)
  ('peano_key'                 ,'i4'      ), #0L $         ,                   , id of small subcube of simulation cube containing galaxy
  ('redshift'                  ,'f4'      ), #0.0 $        ,                   , redshift of galaxy
  ('type'                      ,'i4'      ), #0L $         ,                   , indicated positional status of galaxy in fof group (0 = central, 1 = satellite with subhalo, 2= satellite without resolved subhalo)
  ('snapshot_number'           ,'i4'      ), #0L $         ,                   , simulation snapshot the galaxy belongs to
  ('group_number'              ,'i4'      ), #0L $         ,                   , yet another id of the fof halo the galaxy belongs to
  ('next_group_number'         ,'i4'      ), #0L $         ,                   , yet another id of the fof halo the galaxy will belong to in the next snapshot
  ('cube_index'                ,'i4'      ), #0L $         ,                   , index of periodic copy of simulation cube the galaxy is located
  ('central_m_vir'             ,'f4'      ), #0.0 $        , [10^10 Msun/h]    , virial mass (as defined by m_crit200) of the FOF group the galaxy resides in.
  ('central_r_vir'             ,'f4'      ), #0.0 $        , [Mpc/h]           , virial radius (as defined by r_crit200) of the FOF group the galaxy resides in
  ('position'                  ,'f4',    3), #fltarr(3) $  , [rad, rad, Mpc/h] , angular position (first two components) and line-of-sight comoving distance (last component) of galaxy
  ('velocity'                  ,'f4',    3), #fltarr(3) $  , [km/s]            , physical peculiar velocity of galaxy (first two components transverse, last component parallel to l.o.s.)
  ('len'                       ,'i4'      ), #0L $         ,                   , number of particle in subhalo associated with galaxy
  ('m_vir'                     ,'f4'      ), #0.0 $        , [10^10 Msun/h]    , virial mass (as defined by m_crit200) of the FOF group this galaxy was in when last it was a type 0 galaxy. I.e. current mass for type 0 galaxies, "infall virial mass" for type 1,2 galaxies.
  ('r_vir'                     ,'f4'      ), #0.0 $        , [Mpc/h]           , comoving virial radius (as defined by r_crit200) of the FOF group this galaxy was in when last it was a type 0 galaxy. I.e. current virial radius for type 0 galaxies, "infall virial radius" for type 1,2 galaxies
  ('v_vir'                     ,'f4'      ), #0.0 $        , [km/s]            , physical virial velocity of the subhalo the galaxy is/was the center of.
  ('v_max'                     ,'f4'      ), #0.0 $        , [km/s]            , physical maximum rotational velocity of the subhalo of which this galaxy is the center, or the last value for satellite galaxies.
  ('gas_spin'                  ,'f4',    3), #fltarr(3) $  , [Mpc/h km/s]      , spin of the cold gas disk of galaxy
  ('stellar_spin'              ,'f4',    3), #fltarr(3) $  , [Mpc/h km/s]      , spin of the stellar disk of galaxy
  ('infall_v_max'              ,'f4'      ), #0.0 $        , [km/s]            , physical maximum rotational velocity of the host halo of this galaxy atinfallSnap.
  ('infall_v_max_peak'         ,'f4'      ), #0.0 $        , [km/s]            , physical maximum past rotational velocity of the host halo of this galaxy.
  ('infall_snap'               ,'f4'      ), #0L $         ,                   , id of snapshot the galaxy lost type = 0 status
  ('infall_hot_gas'            ,'f4'      ), #0.0 $        , [10^10 Msun/h]    , mass in hot gas at the time of infall (same as hotGas for type 0 galaxies).
  ('hot_radius'                ,'f4'      ), #0.0 $        , [Mpc/h]           , radius out to which hot gas extends: rvir for type 0; 0 for type 2; maximum radius out to which hot gas is not stripped for type 1.
  ('ori_merg_time'             ,'f4'      ), #0.0 $        , [yr]              , estimated dyniamical friction time (in years) when the merger clock is set.
  ('merg_time'                 ,'f4'      ), #0.0 $        , [yr]              , estimated remaining merging time (in years). oriMergeTime - time since the merger clock is set.
  ('distance_to_central_gal'   ,'f4',    3), #fltarr(3) $  , [Mpc/h (?)]       , distance between this galaxy and the central galaxy of the fof group
  ('cold_gas'                  ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , Mass in the cold gas disk.
  ('stellar_mass'              ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , total mass in stars in the disk and the bulge together.
  ('bulge_mass'                ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass of stars in the bulge.
  ('disk_mass'                 ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass of stars in the disk.
  ('hot_gas'                   ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in hot gas.
  ('ejected_mass'              ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in the ejected mass component
  ('black_hole_mass'           ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass of the central black hole
  ('ICM'                       ,'f4'      ), #0.0 $        , (?)
  ('metals_cold_gas'           ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in cold gas
  ('metals_bulge_mass'         ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in bulge
  ('metals_disk_mass'          ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in disk
  ('metals_hot_gas'            ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in hot gas
  ('metals_ejected_mass'       ,'f4'      ), #0.0 $        , [10^10 Msun /h]   , mass in metals in the ejected mass component
  ('metals_ICM'                ,'f4'      ), #0.0 $        , (?)
  ('primordial_accretion_rate' ,'f4'      ), #0.0 $        , [Msun/yr]         , Accretion rate of primordial gas.
  ('cooling_rate'              ,'f4'      ), #0.0 $        , [Msun/yr]         , cooling rate
  ('cooling_rate_before_AGN'   ,'f4'      ), #0.0 $        , [Msun/yr]         , cooling rate if there was no AGN feedback.
  ('sfr'                       ,'f4'      ), #0.0 $        , [Msun/yr]         , Star formation rate
  ('sfr_bulge'                 ,'f4'      ), #0.0 $        , [Msun/yr]         , Star formation rate in bulge.
  ('x_ray_lum'                 ,'f4'      ), #0.0 $        , [log10(erg/sec)]  , Log10 of X-Ray luminosity in erg/sec
  ('bulge_size'                ,'f4'      ), #0.0 $        , [Mpc/h]           , Half mass radius of bulge
  ('stellar_disk_radius'       ,'f4'      ), #0.0 $        , [Mpc/h]           , Size of the stellar disk, 3x the scale length.
  ('gas_disk_radius'           ,'f4'      ), #0.0 $        , [Mpc/h]           , Size of the gas disk (?)
  ('cos_inclination'           ,'f4'      ), #0.0 $        , (?)
  ('disrupt_on'                ,'i4'      ), #0L $         ,                   , 0: galaxy merged onto merger center; 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center
  ('merge_on'                  ,'i4'      ), #0L $         ,                   , 0: merger clock not set yet;
  ('cooling_radius'            ,'f4'      ), #0.0 $        , [Mpc/h]           , the radius within which the cooling time scale is shorter than the dynamical timescale
  ('quasar_accretion_rate'     ,'f4'      ), #0.0 $        , [Msun/yr]         , Rate at which cold gas is accreted into the central black hole in the quasar mode.
  ('radio_accretion_rate'      ,'f4'      ), #0.0 $        , [Msun/yr]         , Rate at which hot gas is accreted into the central black hole in the radio mode.
  ('mag'                       ,'f4',   40), #fltarr(xx) $ , [mag]             , observer-frame apparent (AB) magnitude of galaxy (magnification included, dust extinction not included ?)
  ('mag_bulge'                 ,'f4',   40), #fltarr(xx) $ , [mag]             , observer-frame apparent (AB) magnitude of galaxy bulge (magnification included, dust extinction not included ???)
  ('mag_dust'                  ,'f4',   40), #fltarr(xx) $ , [mag]             , observer-frame apparent (AB) magnitude of galaxy (magnification included, dust extinction included ?)
  ('mass_weight_age'           ,'f4'      ), #0.0 $        , [10^9 yr]         , The age of this galaxy, weighted by mass of its components.
  ('rband_weight_age'          ,'f4'      ), #0.0 $        , [10^9 yr]         , The age of this galaxy, weighted by mass of its components.
  ('sfh_ibin'                  ,'i4'      ), #0L $         ,                   , Index of the higest star formation history bin currently in use.
  ('sfh_numbins'               ,'i4'      ), #0L $         ,                   , Number of non-empty star formation history bins.
  ('distortion'                ,'f4',(2,2)), #fltarr(4) $  ,                   , (11, 12, 21, 22)-components of distortion matrix
  ('plane_number'              ,'i4'      )  #0L $         ,                   , index of redshift slice (and lens plane) the galaxy is associated with
    ])

    filter_number_for_c_johnson_U      =  0
    filter_number_for_c_johnson_B      =  1
    filter_number_for_c_johnson_V      =  2
    filter_number_for_c_johnson_rc     =  3
    filter_number_for_c_johnson_ic     =  4
    filter_number_for_vista_johnson_Z  =  5
    filter_number_for_vista_johnson_Y  =  6
    filter_number_for_vista_johnson_J  =  7
    filter_number_for_vista_johnson_H  =  8
    filter_number_for_c_johnson_K      =  9
    filter_number_for_vista_johnson_ks = 10
    filter_number_for_i1_band          = 11
    filter_number_for_i2_band          = 12
    filter_number_for_i3_band          = 13
    filter_number_for_i4_band          = 14
    filter_number_for_u_band_trans     = 15
    filter_number_for_g_band_trans     = 16
    filter_number_for_r_band_trans     = 17
    filter_number_for_i_band_trans     = 18
    filter_number_for_z_band_trans     = 19
    filter_number_for_ACS_WFC_F435W    = 20
    filter_number_for_ACS_WFC_F475W    = 21
    filter_number_for_ACS_WFC_F606W    = 22
    filter_number_for_ACS_WFC_F625W    = 23
    filter_number_for_ACS_WFC_F775W    = 24
    filter_number_for_ACS_WFC_F814W    = 25
    filter_number_for_ACS_WFC_F850_LP  = 26
    filter_number_for_GALEX_FUV        = 27
    filter_number_for_GALEX_NUV        = 28
    filter_number_for_NIC_F110W        = 29
    filter_number_for_NIC_F160W3       = 30
    filter_number_for_VIMOS_U          = 31
    filter_number_for_WFC3_IR_F105W    = 32
    filter_number_for_WFC3_IR_F125W    = 33
    filter_number_for_WFC3_IR_F160W    = 34
    filter_number_for_WFC3_UVIS_F225W  = 35
    filter_number_for_WFC3_UVIS_F275W  = 36
    filter_number_for_WFC3_UVIS_F336W  = 37
    filter_number_for_WFPC2_F300W      = 38
    filter_number_for_WFPC2_F450W      = 39

#    with open("/lfs08/rusucs/0408/GGL_los_8_0_0_N_4096_ang_4_Henriques2014_galaxies_on_plane_1_f.images", mode = 'rb') as file:
#        lower_bound = np.fromfile(file, 'f8', 2)
#        upper_bound = np.fromfile(file, 'f8', 2)
#        plane_angle, = np.fromfile(file, 'f8', 1)
#        redshift, = np.fromfile(file, 'f8', 1)
#        n_galaxies, = np.fromfile(file, 'i8', 1)
#        n_cells = np.fromfile(file, 'i4', 2)
#        galaxy = np.fromfile(file, galaxy_struct, n_galaxies)
#        xi = galaxy['mag'][:,filter_number_for_i_band_trans]
#        x0 = galaxy['redshift']
#        x = np.c_[x0,xi]
