# creates FOV figure from Birrer et al. 2018

from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
plt.clf()
pix = 1283 # number of pixels in 4 arcmin
scale = 240.0 / pix
maglim1 = 24
maglim2 = 23
z_l = 1.789
x,y,i,photoz,spec,classify = np.loadtxt("catalogues/i24at2sigma_iunconv_igrKconv_detectin_iunconv_corrisoautoerredit_short_withbpzlephareclass.cat",usecols=[0,1,10,30,42,49],unpack=True)
sep = np.sqrt((x - pix/2.0)**2 + (y - pix/2.0)**2) * scale
redshift = np.copy(photoz)
redshift[spec >= 0] = spec[spec >= 0]

spec = spec[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
photoz = photoz[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
x = x[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
y = y[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
classify = classify[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
sep_ = sep[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
i_ = i[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
redshift = redshift[(sep <= 120) & (i <= maglim1) & (redshift < z_l)]
sep = sep_
i = i_

x_spec_bright = x[((classify == 0) | (classify == 1)) & (i <= maglim2)]
y_spec_bright = y[((classify == 0) | (classify == 1)) & (i <= maglim2)]
x_spec_faint = x[((classify == 0) | (classify == 1)) & (i > maglim2)]
y_spec_faint = y[((classify == 0) | (classify == 1)) & (i > maglim2)]
spec_bright = spec[((classify == 0) | (classify == 1)) & (i <= maglim2)]
spec_faint = spec[((classify == 0) | (classify == 1)) & (i > maglim2)]
redshift_spec_bright = redshift[((classify == 0) | (classify == 1)) & (i <= maglim2)]
redshift_spec_faint = redshift[((classify == 0) | (classify == 1)) & (i > maglim2)]

x_group_bright = x_spec_bright[(spec_bright > 0.736) & (spec_bright <= 0.756)]
y_group_bright = y_spec_bright[(spec_bright > 0.736) & (spec_bright <= 0.756)]
x_group_faint = x_spec_faint[(spec_faint > 0.736) & (spec_faint <= 0.756)]
y_group_faint = y_spec_faint[(spec_faint > 0.736) & (spec_faint <= 0.756)]
redshift_group_bright = redshift_spec_bright[(spec_bright > 0.736) & (spec_bright <= 0.756)]
redshift_group_faint = redshift_spec_faint[(spec_faint > 0.736) & (spec_faint <= 0.756)]

x_star_bright = x[(classify < 0) & (i <= maglim2)]
y_star_bright = y[(classify < 0) & (i <= maglim2)]
x_star_faint = x[(classify < 0) & (i > maglim2)]
y_star_faint = y[(classify < 0) & (i > maglim2)]
redshift_star_bright = redshift[(classify < 0) & (i <= maglim2)]
redshift_star_faint = redshift[(classify < 0) & (i > maglim2)]

x_galnospec_bright = x[(classify == 2) & (i <= maglim2)]
y_galnospec_bright = y[(classify == 2) & (i <= maglim2)]
x_galnospec_faint = x[(classify == 2) & (i > maglim2)]
y_galnospec_faint = y[(classify == 2) & (i > maglim2)]
redshift_galnospec_bright = redshift[(classify == 2) & (i <= maglim2)]
redshift_galnospec_faint = redshift[(classify == 2) & (i > maglim2)]

image = fits.getdata("images/J1206_GMOSi_CFHTLSscale_weighted_bkg96_masked.fits")
image[image<0] = 0.0001
zmax = 0.1*round(10*np.max(redshift)) + 0.1

mask = fits.getdata("catalogues/mskforplot.fits")

plt.scatter(x_star_bright,y_star_bright,marker='*',edgecolors='none',s=30,linewidths=1,c=redshift_star_bright,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_star_faint,y_star_faint,marker='*',edgecolors='none',s=10,linewidths=1,c=redshift_star_faint,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_galnospec_bright,y_galnospec_bright,marker='o',edgecolors='none',s=30,linewidths=1,c=redshift_galnospec_bright,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_galnospec_faint,y_galnospec_faint,marker='o',edgecolors='none',s=10,linewidths=1,c=redshift_galnospec_faint,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_spec_bright,y_spec_bright,marker='s',edgecolors='none',linewidths=1,s=30,c=redshift_spec_bright,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_spec_faint,y_spec_faint,marker='s',edgecolors='none',linewidths=1,s=10,c=redshift_spec_faint,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_group_bright,y_group_bright,marker='s',edgecolors='k',linewidths=1,s=30,c=redshift_group_bright,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.scatter(x_group_faint,y_group_faint,marker='s',edgecolors='k',linewidths=1,s=10,c=redshift_group_faint,cmap = plt.cm.get_cmap("CMRmap"),alpha=0.5,vmin=0,vmax=zmax)
plt.colorbar(format='%.1f')#,boundaries=[0,1.8])
plt.imshow(image, cmap='gray_r', norm=LogNorm(), origin='lower', vmin=0.001, vmax=100)
plt.imshow(mask, cmap='Oranges', origin='lower', alpha=0.2)
circle1 = plt.Circle((pix/2.0,pix/2.0),45/120.0*pix/2.0,color='k',fill=False)
circle2 = plt.Circle((pix/2.0,pix/2.0),pix/2.0,color='k',fill=False)
fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)
fig = plt.gca()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.savefig('FOV_J1206.png', dpi=300, bbox_inches='tight')
