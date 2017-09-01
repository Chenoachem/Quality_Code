import astropy
from astropy.io import fits
from astropy.modeling import models
import numpy as np
import sys
import numpy as np
import matplotlib.pyplot as plt


mycube=sys.argv[1]
slice=sys.argv[2]

datacube = fits.open(mycube)
data = datacube[0].data
header = datacube[0].header
cube = data[:,:,:]

my_slice=int(slice)
counts,bins = np.histogram(cube[:,slice,400:700,400:700],bins=np.arange(-5,5,0.2))
fig2=plt.figure()
ax=fig2.add_subplot(1,2,1)
ax.semilogy(bins[:-1],counts)
#ax.set_xlim([-2,2])
ax.set_title("Pixel Histogram 400:700 - Slice" + slice)
ax2=fig2.add_subplot(1,2,2)
ax2.plot(bins[:-1],counts)
ax2.set_xlim([-4,4])
#ax2.set_title("Pixel Histogram 400:700, 400:700")
plt.tight_layout(pad=0.7)

my_slice=str(my_slice)

fig2.savefig('Histogram_'+my_slice+'.png')
