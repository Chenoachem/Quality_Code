# coding: utf-8

# In[1]:

import os
import numpy as np
import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import wcs
import sys
import os.path
import csv
get_ipython().magic(u'matplotlib inline')


# In[4]:

mycube="Cube_0143_Grouped/Orion_I_143_Subtracted.fits"


# In[61]:

datacube = fits.open(mycube)
data = datacube[0].data
header = datacube[0].header


# In[64]:

my_randomsx = random.sample(xrange(500, 1000), 80)
my_randomsy = random.sample(xrange(500, 1000), 80)
print my_randomsx
print my_randomsy


# In[83]:

spectral_rms=[]
for i in range(0,80):
    signal=[]
    for x in range(0, 100):
        value = np.mean(data[:,x,my_randomsy[i]:my_randomsy[i]+1,my_randomsx[i]:my_randomsx[i]+1])
        signal.append(value)
        RMS=np.nanstd(signal)
    spectral_rms.append(RMS)


# In[84]:

print spectral_rms


# In[90]:

local_rms=[]
for slice in range(0,100):
    value = np.nanstd(data[:,slice,300:600,300:600])
    local_rms.append(value)


# In[91]:

print local_rms


# In[92]:

local_rms=np.array(local_rms)
spectral_rms=np.array(spectral_rms)


# In[93]:

ratio=[]
for value in range(0,80):
    rat=(np.divide(local_rms,spectral_rms[value])
    ratio.append(rat)


# In[94]:

print ratio


# In[95]:

print len(ratio)


# In[98]:

print len(np.ravel(ratio))


# In[99]:

ratio_list=np.ravel(ratio)


# In[100]:

passing=np.count_nonzero((0.9 < ratio_list) & (ratio_list < 1.1))
print passing


# In[104]:

print 1205./8000.
