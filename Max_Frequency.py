#!/bin/python

#  Max_Frequency.py
#  
#
#  Created by Chenoa Tremblay on 31/08/2017.
#
import os
import argparse
import numpy as np
import astropy
from astropy.io import fits
import matplotlib
from matplotlib import pyplot
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import wcs
import math
import aplpy
matplotlib.use('Agg')
from astropy.io.votable import parse_single_table
import sys
import os.path

mycube=sys.argv[1]
#ra=sys.argv[2]
#dec=sys.argv[3]
#slice=sys.argv[4]
cat_files=sys.argv[2]

file=open(cat_files)
for line in file:
    ra,dec = line.strip().split(",")
    ra=str(ra)
    dec=str(dec)

    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    RA=c.ra.deg
    Dec=c.dec.deg
    Gal=c.galactic

    datacube = fits.open(mycube)
    data = datacube[0].data
    header = datacube[0].header


    rp = datacube[0].header['CRPIX3']
    rf = datacube[0].header['CRVAL3']
    df = datacube[0].header['CDELT3']
    nf = datacube[0].header['NAXIS3']
    frequency=rf + df*(np.arange(nf)-rp)

    #bmaj=header['BMAJ']
    #bmin=header['BMIN']
    #cdelt1=header['CDELT1']
    #cdelt2=header['CDELT2']
    #beam_pix=math.sqrt((bmaj*bmin)/(-cdelt1*cdelt2))
    beam_pix=4.5
    w = wcs.WCS(header, naxis=2)
    
    ypix,xpix=c.to_pixel(w,origin=0,mode='wcs')

    xpix1=xpix-(beam_pix/2)
    xpix2=xpix+(beam_pix/2)
    ypix1=ypix-(beam_pix/2)
    ypix2=ypix+(beam_pix/2)
    
    signal=[]
    for x in range(0, 100):
        value = np.nanmean(data[:,x,xpix1:xpix2,ypix1:ypix2])
        #print rms_number_x
        signal.append(value)

    signal=np.multiply(signal,1.13)
    RMS=np.nanstd(signal)

    max_signal=np.nan_to_num(signal)
    max_value=np.max(max_signal)
    sli=np.argmax(max_signal, axis=0)

    if os.path.exists('freq'):
    	with open('freq','a') as f:
            f.write(str(frequency[sli]/10**6)+"\n")
            f.close()
    print frequency[sli]/10**6
    if os.path.exists('RMS'):
        with open('RMS','a') as p:
            p.write(str(RMS)+"\n")
            p.close()

    print RMS

    if os.path.exists('snr'):
        with open('snr','a') as q:
            q.write(str(max_value/RMS)+"\n")
            q.close()


    print max_value/RMS
