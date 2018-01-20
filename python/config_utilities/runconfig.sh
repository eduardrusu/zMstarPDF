#! /bin/sh
#Builds the template libraries and puts them in $LEPHAREWORK/

LEPHAREDIR=$HOME/lephare_dev
PARA=$HOME/lephare_dev/config/lephare_CFHTLenS.para

LEPHAREWORK=$HOME/lephare_dev/work

export LEPHAREDIR
export LEPHAREWORK

if [ ! -d $LEPHAREWORK  ]; then
    echo "creating $LEPHAREWORK..."
    mkdir $LEPHAREWORK
    mkdir $LEPHAREWORK/filt
    mkdir $LEPHAREWORK/lib_mag
    mkdir $LEPHAREWORK/lib_bin
fi

#Filter transmission curves
$LEPHAREDIR/source/filter -c $PARA -FILTER_LIST cfht/megacam/CFHTLS_u.pb,cfht/megacam/CFHTLS_g.pb,cfht/megacam/CFHTLS_r.pb,cfht/megacam/CFHTLS_i.pb,cfht/megacam/CFHTLS_y.pb,cfht/megacam/CFHTLS_z.pb
#Galaxy templates -------------------------------#
#CWW
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB    LIB_CWW -GAL_SED $LEPHAREDIR/sed/GAL/CWW_KINNEY/CWW_MOD.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_CWW -GAL_LIB_OUT MAG_CWW -MOD_EXTINC 3,10 -EB_V 0.,0.05,0.1,0.15,0.2,0.25,0.3

#NOT FOUND:
#CWW - improved Ilbert et al. (2006) templates + optimized for CFHTLS and interpolated
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB    LIB_CWW_CFHTLS -GAL_SED $LEPHAREDIR/sed/GAL/CFHTLS_WIDE/CE_MOD.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_CWW_CFHTLS -GAL_LIB_OUT MAG_CWW_CFHTLS -MOD_EXTINC 40,70 -EXTINC_LAW SMC_prevot.dat -EB_V 0.000,0.05,0.1,0.15,0.2,0.25

#NOT FOUND:
#CWW - same as obove but slightly more recent
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB    LIB_CWW_CFHTLS_NEW -GAL_SED $LEPHAREDIR/sed/GAL/PHOTO_230506/CE_MOD.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_CWW_CFHTLS_NEW -GAL_LIB_OUT MAG_CWW_CFHTLS_NEW -MOD_EXTINC 38,100 -EXTINC_LAW SMC_prevot.dat -EB_V 0.000,0.05,0.1,0.15,0.2

#CWW - COSMOS
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB    LIB_CWW_COSMOS  -GAL_SED     $LEPHAREDIR/sed/GAL/COSMOS_SED/COSMOS_MOD.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_CWW_COSMOS  -GAL_LIB_OUT MAG_CWW_COSMOS -MOD_EXTINC 13,23,23,31,23,31,23,31 -EXTINC_LAW SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat -EB_V 0.000,0.100,0.200,0.300,0.400,0.500

#Bruzual & Charlot - Atsushi
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB    LIB_BC03 -GAL_SED $LEPHAREDIR/sed/GAL/HYPERZ/HYPERZ_GIS_MOD.list -SEL_AGE $LEPHAREDIR/sed/GAL/BC03_HSC/AGE_GISSEL_HZ.dat
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_BC03 -GAL_LIB_OUT MAG_BC03

#Bruzual & Charlot - phys para
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB  LIB_BC03 -GAL_SED $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_MOD.list -SEL_AGE $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_BC03 -GAL_LIB_OUT MAG_BC03 -EXTINC_LAW calzetti.dat -EB_V 0.,0.05,0.1,0.15,0.2,0.25,0.3 -MOD_EXTINC 0,100

#Bruzual & Charlot - phys para      for me, to be able to run without recompiling
#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB  LIB_BC03 -GAL_SED $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_MOD.list -SEL_AGE $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_BC03 -GAL_LIB_OUT MAG_BC03 -EXTINC_LAW calzetti.dat -EB_V 0.,0.1,0.2,0.3,0.4,0.5 -MOD_EXTINC 0,100

#Bruzual & Charlot - Ilbert et al. 2009
$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB  LIB_BC03_I09 -GAL_SED $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_I09_MOD.list -SEL_AGE $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list
$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_BC03_I09 -GAL_LIB_OUT MAG_BC03_I09 -EXTINC_LAW calzetti.dat -EB_V 0.,0.1,0.2,0.3,0.4,0.5 -MOD_EXTINC 0,20

#$LEPHAREDIR/source/sedtolib -c $PARA -t G -GAL_LIB  LIB_BC03_NO_EXT -GAL_SED $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_MOD.list -SEL_AGE $LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list
#$LEPHAREDIR/source/mag_gal  -c $PARA -t G -GAL_LIB_IN LIB_BC03_NO_EXT -GAL_LIB_OUT MAG_BC03_NO_EXT  -EB_V 0.

#Star templates ---------------------------------#
$LEPHAREDIR/source/sedtolib -c $PARA -t S -STAR_LIB    LIB_STAR -STAR_SED $LEPHAREDIR/sed/STAR/STAR_MOD.list
$LEPHAREDIR/source/mag_star -c $PARA -t Q -STAR_LIB_IN LIB_STAR -STAR_LIB_OUT MAG_STAR

#QSO templates ----------------------------------#
$LEPHAREDIR/source/sedtolib -c $PARA -t Q -QSO_LIB    LIB_QSO -QSO_SED $LEPHAREDIR/sed/QSO/QSO_MOD.list
$LEPHAREDIR/source/mag_gal  -c $PARA -t Q -QSO_LIB_IN LIB_QSO -QSO_LIB_OUT MAG_QSO
