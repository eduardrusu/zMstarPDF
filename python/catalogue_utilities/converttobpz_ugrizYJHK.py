##########################
# The code is used to convert raw matched photometric catalogues (with matched PSFs) into BPZ-expected input files, except that the ID, image/world coordinates, flux radii and flags are still included. The code uses the measured 1-sigma mag limits to treat non-detections (not yet non-exposures). The code applies a MAG_ISO -> MAG_AUTO conversion. The conversion is as follows:
# for detections in i+r: auto_x = iso_x + average(auto_r_noconv - iso_r, auto_i_noconv - iso_i)
# for detections in i: auto_x = iso_x + auto_i_noconv - iso_i
# This expresion is useful to find total mags for photoz (with luminosity priors) and stellar mass computation, but ignores color gradients, as described in Erben et al. 2014
##########################

# Careful, because the columns ID changes depending on the input file; also, change the 1-sigma detection limits as appropriate
# Before running, make sure that the mag_iso_i columns in the convolved frame contains no 99 or -99, because it is used for the conversions; in the case of detections in ir only, also mag_iso_r

import numpy as np 

file = "/Users/eduardrusu/Desktop/WFI2033/WFI2033analysis/rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23.cat"
#file = "/Users/eduardrusu/Desktop/WFI2033/WFI2033analysis/i_detect_i_and_ir_rnoconv_inoconv_ugrizYJHK_detectin_ir_short_potentiallyi23.cat"

# I calculated the 5-sigma limits, now convert to 1-sigma
limit_u = 24.25 + 1.74 # AB; u-band is not updated
limit_g = 23.77 + 1.74 # AB
limit_r = 23.81 + 1.74 # AB
limit_i = 23.18 + 1.74 # AB
limit_z = 22.91 + 1.74 # AB
limit_Y = 21.40 + 1.74 # AB
limit_J = 20.63 + 1.74 # Vega
limit_H = 19.52 + 1.74 # Vega
limit_K = 18.87 + 1.74 # Vega

if "detect_i_and_ir" in file:
    xpix = 0
    ypix = 1
    xwcs = 2
    ywcs = 3
    id = 4 # detections in i
    flx_rad = 8 # measurement in i, detection in i
    flg = 16
    auto_i_noconv_detecti = 12
    auto_i_noconv_detecti_err = 13
    iso_u = 41
    iso_u_err = 42
    iso_g = 45
    iso_g_err = 46
    iso_r = 49
    iso_r_err = 50
    iso_i = 53
    iso_i_err = 54
    iso_z = 57
    iso_z_err = 58
    iso_Y = 61
    iso_Y_err = 62
    iso_J = 65
    iso_J_err = 66
    iso_H = 69
    iso_H_err = 70
    iso_K = 73
    iso_K_err = 74
else:
    xpix = 0
    ypix = 1
    xwcs = 2
    ywcs = 3
    id = 4 # detections in i
    flx_rad = 17 # measurement in i
    flg = 13
    auto_r_noconv = 9
    auto_r_noconv_err = 10
    auto_i_noconv = 18
    auto_i_noconv_err = 19
    iso_u = 24
    iso_u_err = 25
    iso_g = 28
    iso_g_err = 29
    iso_r = 32
    iso_r_err = 33
    iso_i = 36
    iso_i_err = 37
    iso_z = 40
    iso_z_err = 41
    iso_Y = 44
    iso_Y_err = 45
    iso_J = 48
    iso_J_err = 49
    iso_H = 52
    iso_H_err = 53
    iso_K = 56
    iso_K_err = 57

data = np.loadtxt(file,unpack=True)

# m =  99.0, dm = 28.4 -- example of a non-detection: dm = 1-sigma detection limit
data[iso_u_err][data[iso_u] > limit_u] = limit_u
data[iso_u][data[iso_u] > limit_u] = 99.0
data[iso_g_err][data[iso_g] > limit_g] = limit_g
data[iso_g][data[iso_g] > limit_g] = 99.0
data[iso_r_err][data[iso_r] > limit_r] = limit_r
data[iso_r][data[iso_r] > limit_r] = 99.0
data[iso_i_err][data[iso_i] > limit_i] = limit_i
data[iso_i][data[iso_i] > limit_i] = 99.0
data[iso_z_err][data[iso_z] > limit_z] = limit_z
data[iso_z][data[iso_z] > limit_z] = 99.0
data[iso_Y_err][data[iso_Y] > limit_Y] = limit_Y
data[iso_Y][data[iso_Y] > limit_Y] = 99.0
data[iso_J_err][data[iso_J] > limit_J] = limit_J
data[iso_J][data[iso_J] > limit_J] = 99.0
data[iso_H_err][data[iso_H] > limit_H] = limit_H
data[iso_H][data[iso_H] > limit_H] = 99.0
data[iso_K_err][data[iso_K] > limit_K] = limit_K
data[iso_K][data[iso_K] > limit_K] = 99.0

# implement the iso -> auto correction, but only for mags != 99,-99
if "detect_i_and_ir" in file:
    isoauto = data[auto_i_noconv_detecti] - data[iso_i]
    isoauto[np.abs(data[iso_i]) == 99.0] = 0
    data[iso_u][np.abs(data[iso_u]) != 99.0] = data[iso_u][np.abs(data[iso_u]) != 99.0] + isoauto[np.abs(data[iso_u]) != 99.0]
    data[iso_g][np.abs(data[iso_g]) != 99.0] = data[iso_g][np.abs(data[iso_g]) != 99.0] + isoauto[np.abs(data[iso_g]) != 99.0]
    data[iso_r][np.abs(data[iso_r]) != 99.0] = data[iso_r][np.abs(data[iso_r]) != 99.0] + isoauto[np.abs(data[iso_r]) != 99.0]
    data[iso_i][np.abs(data[iso_i]) != 99.0] = data[iso_i][np.abs(data[iso_i]) != 99.0] + isoauto[np.abs(data[iso_i]) != 99.0]
    data[iso_z][np.abs(data[iso_z]) != 99.0] = data[iso_z][np.abs(data[iso_z]) != 99.0] + isoauto[np.abs(data[iso_z]) != 99.0]
    data[iso_Y][np.abs(data[iso_Y]) != 99.0] = data[iso_Y][np.abs(data[iso_Y]) != 99.0] + isoauto[np.abs(data[iso_Y]) != 99.0]
    data[iso_J][np.abs(data[iso_J]) != 99.0] = data[iso_J][np.abs(data[iso_J]) != 99.0] + isoauto[np.abs(data[iso_J]) != 99.0]
    data[iso_H][np.abs(data[iso_H]) != 99.0] = data[iso_H][np.abs(data[iso_H]) != 99.0] + isoauto[np.abs(data[iso_H]) != 99.0]
    data[iso_K][np.abs(data[iso_K]) != 99.0] = data[iso_K][np.abs(data[iso_K]) != 99.0] + isoauto[np.abs(data[iso_K]) != 99.0]
else:
    isoauto = np.mean([data[auto_i_noconv] - data[iso_i],data[auto_r_noconv] - data[iso_r]],axis=0)
    isoauto[np.abs(data[iso_i]) == 99.0] = 0
    isoauto[(np.abs(data[iso_i]) != 99.0) & ((np.abs(data[auto_r_noconv]) == 99.0) | (np.abs(data[iso_r]) == 99.0))]  = data[auto_i_noconv][(np.abs(data[iso_i]) != 99.0) & ((np.abs(data[auto_r_noconv]) == 99.0) | (np.abs(data[iso_r]) == 99.0))] - data[iso_i][(np.abs(data[iso_i]) != 99.0) & ((np.abs(data[auto_r_noconv]) == 99.0) | (np.abs(data[iso_r]) == 99.0))]
    data[iso_u][np.abs(data[iso_u]) != 99.0] = data[iso_u][np.abs(data[iso_u]) != 99.0] + isoauto[np.abs(data[iso_u]) != 99.0]
    data[iso_g][np.abs(data[iso_g]) != 99.0] = data[iso_g][np.abs(data[iso_g]) != 99.0] + isoauto[np.abs(data[iso_g]) != 99.0]
    data[iso_r][np.abs(data[iso_r]) != 99.0] = data[iso_r][np.abs(data[iso_r]) != 99.0] + isoauto[np.abs(data[iso_r]) != 99.0]
    data[iso_i][np.abs(data[iso_i]) != 99.0] = data[iso_i][np.abs(data[iso_i]) != 99.0] + isoauto[np.abs(data[iso_i]) != 99.0]
    data[iso_z][np.abs(data[iso_z]) != 99.0] = data[iso_z][np.abs(data[iso_z]) != 99.0] + isoauto[np.abs(data[iso_z]) != 99.0]
    data[iso_Y][np.abs(data[iso_Y]) != 99.0] = data[iso_Y][np.abs(data[iso_Y]) != 99.0] + isoauto[np.abs(data[iso_Y]) != 99.0]
    data[iso_J][np.abs(data[iso_J]) != 99.0] = data[iso_J][np.abs(data[iso_J]) != 99.0] + isoauto[np.abs(data[iso_J]) != 99.0]
    data[iso_H][np.abs(data[iso_H]) != 99.0] = data[iso_H][np.abs(data[iso_H]) != 99.0] + isoauto[np.abs(data[iso_H]) != 99.0]
    data[iso_K][np.abs(data[iso_K]) != 99.0] = data[iso_K][np.abs(data[iso_K]) != 99.0] + isoauto[np.abs(data[iso_K]) != 99.0]

fileout = file[:-4] + "_forbpz_NEW.cat"

if "detect_i_and_ir" in file:
    str = "X_IMAGE Y_IMAGE X_WORLD Y_WORLD MAG_AUTO_i_noconv_detect_i MAG_AUTO_ERR_i_noconv_detect_i FLUX_RADIUS_i_noconv_detect_i FLAG_i_noconv_detect_i #ID u u_err g g_err r r_err i i_err z z_err Y Y_err J J_err H H_err K K_err"
    dataout = np.c_[data[xpix],data[ypix],data[xwcs],data[ywcs],data[auto_i_noconv_detecti],data[auto_i_noconv_detecti_err],data[flx_rad],data[flg],data[id],data[iso_u],data[iso_u_err],data[iso_g],data[iso_g_err],data[iso_r],data[iso_r_err],data[iso_i],data[iso_i_err],data[iso_z],data[iso_z_err],data[iso_Y],data[iso_Y_err],data[iso_J],data[iso_J_err],data[iso_H],data[iso_H_err],data[iso_K],data[iso_K_err]]
    np.savetxt(fileout,dataout,header=str,fmt='%.4f \t %.4f \t %.9f \t %.9f \t %.2f \t %.2f \t %.3f \t %d \t %d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f')
else:
    str = "X_IMAGE Y_IMAGE X_WORLD Y_WORLD MAG_AUTO_i_noconv MAG_AUTO_ERR_i_noconv FLUX_RADIUS_i_noconv FLAG #ID u u_err g g_err r r_err i i_err z z_err Y Y_err J J_err H H_err K K_err"
    dataout = np.c_[data[xpix],data[ypix],data[xwcs],data[ywcs],data[auto_i_noconv],data[auto_i_noconv_err],data[flx_rad],data[flg],data[id],data[iso_u],data[iso_u_err],data[iso_g],data[iso_g_err],data[iso_r],data[iso_r_err],data[iso_i],data[iso_i_err],data[iso_z],data[iso_z_err],data[iso_Y],data[iso_Y_err],data[iso_J],data[iso_J_err],data[iso_H],data[iso_H_err],data[iso_K],data[iso_K_err]]
    np.savetxt(fileout,dataout,header=str,fmt='%.4f \t %.4f \t %.9f \t %.9f \t %.2f \t %.2f \t %.3f \t %d \t %d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f')
