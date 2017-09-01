#!/usr/bin/env python
import os
import aplpy
import numpy as np

import matplotlib
matplotlib.use('Agg')



import matplotlib.pyplot as pyplot


from astropy.io.votable import parse_single_table

#datadir="/home/tash/Downloads/"
# Chenoa's computer:
import astropy
from astropy.io import fits
data, header = fits.getdata("Orion_rms.fits", header=True)
header.remove('CRPIX4')
header.remove('CDELT4')
header.remove('CRVAL4')
header.remove('CTYPE4')
fits.writeto("Orion_rms_fix.fits", data, header, clobber=True)
#rm_cube = data[0,:,:,:].reshape((99,2000,2000))

bigfig=pyplot.figure(figsize=(30,30))

#RA=272.0067
#Dec=-29.895
#fov_zoom=5
vmin=0.15
vmax=1.15
fov_mwa=0.5
fov_wise=0.07
s_plot=3500

#NO1posloc=np.array([265.7333,-24.7247])
#NO2posloc=np.array([256.000,-36.8717])
#NO3posloc=np.array([267.6667,-34.1253])
#SHposloc=np.array([264.0417,-25.5642])

gc=aplpy.FITSFigure('Orion_rms_fix.fits', figure=bigfig)
gc.show_colorscale(vmin=vmin, vmax=vmax, stretch='linear', cmap='cubehelix_r')

gc.add_colorbar()
gc.colorbar.set_axis_label_text('RMS (Jy beam$^{-1}$)')
gc.colorbar.set_width(0.8)
gc.colorbar.set_location('right')
gc.colorbar.set_pad(1)
gc.colorbar.set_font(size='40')
gc.colorbar.set_axis_label_font(size='40')
gc.show_contour(levels=(0.3, 0.5), colors='black', linewidths=1, smooth=75, kernel='gauss' )

#gc.show_markers(NO1posloc[0], NO1posloc[1], edgecolor='red', linewidth=2, marker='*', s=s_plot)
#gc.show_markers(NO2posloc[0], NO2posloc[1], edgecolor='red', linewidth=2, marker='*', s=s_plot)
#gc.show_markers(NO3posloc[0], NO3posloc[1], edgecolor='red', linewidth=2, marker='*', s=s_plot)
#gc.show_markers(SHposloc[0], SHposloc[1], edgecolor='red', linewidth=2, marker='d',  s=s_plot)
#gc.set_title('RMS image at 107MHz', fontsize='50')

gc.tick_labels.set_xformat("hh")
gc.tick_labels.set_yformat("dd")
gc.tick_labels.set_font(size='50')

gc.axis_labels.set_xtext('Right Ascension (J2000)')

gc.axis_labels.set_ytext('Declination (J2000)')
gc.axis_labels.set_font(size=50, weight='heavy')
gc.axis_labels.set_xpad(0.5)
#gc.recenter(RA,Dec)
gc.set_tick_color('black')

bigfig.savefig('RMS.png',pad_inches=0.25,bbox_inches='tight')
