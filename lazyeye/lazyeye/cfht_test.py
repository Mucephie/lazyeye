import numpy as np
from astropy.table import vstack, join, Table
from astropy.io import fits
import astropy.coordinates as acoords
from astropy.wcs import WCS
import os
from os import path
import pathlib
import time
import datetime

from astropy import units as un
from astroquery.sdss import SDSS
from astroquery.cadc import Cadc

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import cm
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, make_lupton_rgb)

import astroalign as aa

import lazyeye_core as le


START_DATE_TIME = datetime.datetime.now()

print('\n\n\nStarting time: ', START_DATE_TIME)

print('\nFinding data...\n\n')



ra = 248.042385175791
dec = 26.0552014129828

# ra = 178.454228
# dec = 52.326782

# ra = 179.823687
# dec = 52.686

# ra = 179.399838
# dec = 53.374737

# ra = 161.653
# dec = 63.223

# ra = 148.888
# dec = 69.065


# aco = acoords.SkyCoord(sheet_table['ra'], sheet_table['dec'] , unit='deg')
aco = acoords.SkyCoord(ra, dec , unit='deg')



cadc = Cadc()
result = cadc.query_region(aco, collection='CFHT')
# result = cadc.query_region('08h45m07.5s +54d18m00s', collection='CFHT')
resulttab = Table(result.filled(), masked = False)
resulttab.remove_column('position_bounds_samples')
resulttab.remove_column('position_bounds')
resulttab.show_in_browser(jsviewer=True)
 


# images = le.sdss_images(aco[0])
# images = le.sdss_images(aco, band = ['r', 'g', 'u'])
images = le.mega_images(aco)
# images = le.sdss_images(aco, band = ['z', 'i', 'r'])

for image in images:
        print(image)
        image.info()
        plt.imshow(image, origin='lower')
        print('\n')


image =  images[3]
image.info()

plt.imshow(image.data, origin='lower')
hdu = image
wcs = WCS(hdu.header)

num_sigma = 2
data = hdu.data

mean = np.mean(data)
std = np.std(data)
mean_up = mean + num_sigma * std
mean_down = mean - num_sigma * std

plt.subplot(projection=wcs)
# plt.imshow(hdu.data, origin='lower', cmap = 'viridis', vmin = mean_down, vmax = mean_up )

plt.imshow(image, origin='lower')

plt.grid(color='white', ls='solid', linewidth = 0.3)
plt.xlabel('ra')
plt.ylabel('dec')
# cbar = plt.colorbar()
plt.show()

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))