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

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import cm
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize, make_lupton_rgb)

import astroalign as aa

# Home built imports
# import search as ispy
# import SDSS_Qdep as sdssQ
# import photometry as phot
# import YES
# import cosmog as wanda
# import geometry as geom
# import graphing as graf

import lazyeye_core as le


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')

# strdir = os.getcwd()
# print('Data found in ', strdir, '\n')
# lsa_main = Table.read('lsa_150xxxx_catalogue.csv')
# lsaIDs = np.array(Table.read('lsa_cat_mod-2.csv')['Sheet'])

# lsa_non_hc = Table.read('SDSS_DR9_NON_HC.csv')

# u = 13
# # u = 6 -> 1502468
# # u = 13 -> 1502833
# # u = 9 -> 1503132
# sheetu = lsa_main[lsa_main['sheet']==lsaIDs[u]]

# sheet_nonhc = sheetu[sheetu['Hid']>700000]['index']

# for nhc in range(len(sheet_nonhc)):
#         idx = int(sheet_nonhc[nhc] - 1)
#         gal_nhc = lsa_non_hc[lsa_non_hc['ï»¿ID'] == sheetu[idx]['Hid']]
#         sheetu[idx]['ra'] = gal_nhc['ra'] 
#         sheetu[idx]['dec'] = gal_nhc['dec'] 
#         sheetu[idx]['z'] = gal_nhc['z']
#         sheetu[idx]['sID'] = gal_nhc['SpecObjID']
 
# sheetu['z'].name = 'z*'
# sheetu['ra_su'] = sheetu['ra']
# sheetu['dec_su'] = sheetu['dec']
# sheetu.show_in_browser() ##

# lsaID = str(lsaIDs[u])
# print('Sheet: ', lsaID)

# ra = np.array(sheetu['ra_su']).astype(np.float)
# dec = np.array(sheetu['dec_su']).astype(np.float)

# co = acoords.SkyCoord(ra, dec, unit='deg')

# radi = 3
# drn = 14
# sdss_arr = le.sdss_find(co, radi, drn)

# galIDs = np.array(sdss_arr['ObjID']).astype(np.float)
# galSIDs = np.array(sdss_arr['specObjID']).astype(np.float)

# sdss_arr['ra_dr16'] = sdss_arr['ra']
# sdss_arr['dec_dr16'] = sdss_arr['dec']

# sdss_arr['ra'] = np.around(sdss_arr['ra'], decimals=3)
# sdss_arr['dec'] = np.around(sdss_arr['dec'], decimals=3)

# sheetu['ra'] = np.around(sheetu['ra'], decimals=3)
# sheetu['dec'] = np.around(sheetu['dec'], decimals=3)

# sdss_arr['type'].name = 'sdss_type'

# sdss_arr.show_in_browser() ##

# sheet_table = join(sheetu, sdss_arr, keys = ['ra', 'dec'], join_type = 'inner')
# sheet_table.remove_columns(['ra', 'dec', 'specObjID1', 'sID'])
# sheet_table['ra'] = sheet_table['ra_su']
# sheet_table['dec'] = sheet_table['dec_su']
# sheet_table.show_in_browser() ##

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

# images = le.sdss_images(aco[0])
# images = le.sdss_images(aco, band = ['r', 'g', 'u'])
images = le.sdss_images(aco, band = ['i', 'r', 'g'])
# images = le.sdss_images(aco, band = ['z', 'i', 'r'])

image_r =  images[0]
image_g = images[1]
image_b = images[2]



img_g, footprint_g = aa.register(image_g[0].data, image_r[0].data)
img_b, footprint_b = aa.register(image_b[0].data, image_r[0].data)

mean_r = np.mean(image_r[0].data)
mean_g = np.mean(img_g)
mean_b = np.mean(img_b)

img_g = img_g * (mean_r/mean_g)
img_b = img_b * (mean_r/mean_b)  *  0.7

image_01 = images[0]
# print('Image 1 shape: ', image_01)

hdu = image_r[0]
wcs = WCS(hdu.header)



# image = make_lupton_rgb(image_r[0].data, image_g[0].data, image_b[0].data, stretch=0.2)


image = make_lupton_rgb(image_r[0].data, img_g, img_b, stretch=0.4)


num_sigma = 2
data = hdu.data
mean = np.mean(data)
std = np.std(data)
mean_up = mean + num_sigma * std
mean_down = mean - num_sigma * std

plt.subplot(projection=wcs)
# plt.imshow(hdu.data, origin='lower', cmap = 'viridis', vmin = mean_down, vmax = mean_up )
# rgb = make_lupton_rgb(image_r, image_g, image_b, Q=10, stretch=0.5, filename="test.jpeg")
plt.imshow(image, origin='lower')

plt.grid(color='white', ls='solid', linewidth = 0.3)
plt.xlabel('ra')
plt.ylabel('dec')
cbar = plt.colorbar()
plt.show()

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))