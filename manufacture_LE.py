#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from matplotlib.image import imread
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import rotate
from scipy.ndimage import gaussian_filter
import random

#Options
#Decrease brightness
dec_bright = False
factor_dimmer = 10**2

#Rotate image
rotate_image = True

#Move LE around image
move_LE = True

#Load random existing snle image
snle_image_paths = os.listdir('snle_updated_25_1_24/')
snle_image_paths = ['snle_updated_25_1_24/'+x for x in snle_image_paths if x.endswith(".fits")]
snle_filename = snle_image_paths[random.randint(0, len(snle_image_paths)-1)]
snle = fits.open(snle_filename)[0]
snle_data = snle.data

#Load corresponding snle mask
wo_folder = snle_filename.split('/')[1].split('.')[0]+"_"+snle_filename.split('_')[6].split('.')[0]
mask = 'snle_masks_updated_25_1_24/'+wo_folder+'_mask.png'
mask_raw = imread(mask)
mask_image = np.dot(mask_raw[...,:4], [1, 1, 1, 1])
snle_name = mask.split('snle_masks_updated_25_1_24/')[1].split('_mask')[0]

#Load image you want to insert the LE into
crab_image_paths = os.listdir('DF_mod/')
crab_image_paths = ['DF_mod/'+x for x in crab_image_paths if x.endswith(".fits")]
crab_filename = crab_image_paths[random.randint(0, len(crab_image_paths)-1)]
crab = fits.open(crab_filename)[0]
crab_data = crab.data
crab_name = crab_filename.split('DF_mod/')[1].split('_r_')[0]

#Fix errors in header
snle.header['RADESYSa'] = snle.header['RADECSYS']
del snle.header['RADECSYS']
del snle.header['MJD-OBS']
del snle.header['DATE-OBS']

#Removes background from snle image
snle_data[mask_image == 0] = 0
snle_data = snle_data[~np.all(mask_image == 0, axis=1)]
mask = mask_image[~np.all(mask_image == 0, axis=1)]
snle_data = np.delete(snle_data, np.argwhere(np.all(mask == 0, axis=0)), axis=1)

#Pad edges of snle image, so it is the same shape as background image
df_y, df_x = (3190, 4024)
snle_y, snle_x = snle_data.shape
x = np.floor((df_x - snle_x)/2).astype(int)
y = np.floor((df_y - snle_y)/2).astype(int)
if (snle_y + 2*y != df_y) and (snle_x + 2*x != df_x):
    snle_pad = np.pad(snle_data, ((y+1, y), (x+1, x)), 'constant', constant_values=(0, 0))
elif (snle_y + 2*y != df_y):
    snle_pad = np.pad(snle_data, ((y+1, y), (x, x)), 'constant', constant_values=(0, 0))
elif (snle_x + 2*x != df_x):
    snle_pad = np.pad(snle_data, ((y, y), (x+1, x)), 'constant', constant_values=(0, 0))
else:
    snle_pad = np.pad(snle_data, ((y, y), (x, x)), 'constant', constant_values=(0, 0))
    
#Plot Original images
ax1 = plt.subplot(2,2,1, projection=WCS(crab.header))
ax1.imshow(crab_data, cmap='gray', origin='lower', norm='log')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Original Crab Image\n%s' %crab_name)

ax2 = plt.subplot(2,2,2, projection=WCS(snle.header))
ax2.imshow(snle_pad, cmap='gray', origin='lower', norm='log')
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.set_title('LE with Padding\n%s' %snle_name)

#Reprojection
snle.header['CRVAL1'] = crab.header['CRVAL1']
snle.header['CRVAL2'] = crab.header['CRVAL2']
snle.data = snle_pad
array, footprint = reproject_interp(snle, crab.header)
array[np.isnan(array)] = 0
array[(-0.1 < array) & (array < 0.1)] = 0 #Map bad pixels to 0

#Optional augmentations
#Decrease brightness
if dec_bright==True:
    array = array/factor_dimmer

#Rotate LE by random angle
if rotate_image==True:
    angle=random.randint(0, 360)
    array = rotate(array, angle, reshape=False)
    array[(-0.1 < array) & (array < 0.1)] = 0 #Map bad pixels to 0

#Randomly move LE around image
if move_LE==True:
    max_y, max_x = crab_data.shape
    index = np.where(array != 0)
    x_LE_min = np.min(index[1])
    x_LE_max = np.max(index[1])
    y_LE_min = np.min(index[0])
    y_LE_max = np.max(index[0])
    move_dim = [random.randint(-y_LE_min, max_y-y_LE_max), random.randint(-x_LE_min, max_x-x_LE_max)] #y, x
    array = np.roll(array, move_dim, axis=(0, 1))
    
#Smoothing to match resolution
array = gaussian_filter(array, sigma=1)

#Plot reprojected image
ax3 = plt.subplot(2,2,3, projection=WCS(crab.header))
ax3.imshow(array, cmap='gray', origin='lower', norm='log')
ax3.coords['ra'].set_axislabel('Right Ascension')
ax3.coords['dec'].set_axislabel('Declination')
ax3.set_title('Reprojected LE')

#Overlap images
overlap = array + crab_data

#Plot overlapped image
ax4 = plt.subplot(2,2,4, projection=WCS(crab.header))
ax4.imshow(overlap, cmap='gray', origin='lower', norm='log')
ax4.coords['ra'].set_axislabel('Right Ascension')
ax4.coords['dec'].set_axislabel('Declination')
ax4.set_title('Overlapped Images')
plt.subplots_adjust(wspace=0.8, hspace=0.8)
#plot_name = 'Process_%s_%s.png' %(crab_name, snle_name)
#plt.savefig(plot_name)
#plt.close('all')

#Create new mask
new_mask = np.zeros(array.shape)
new_mask[array != 0] = 1
new_mask = np.flip(new_mask, axis=0)

#Plot new mask
m_ax = plt.subplot(1,1,1, projection=WCS(crab.header))
m_ax.imshow(new_mask, cmap='gray', origin='lower')
m_ax.coords['ra'].set_axislabel('Right Ascension')
m_ax.coords['dec'].set_axislabel('Declination')
mask_name = 'mask_%s_%s' %(crab_name, snle_name)
plt.savefig('%s.png' %mask_name)
plt.close('all')

#Save mask
np.save('mask_%s_%s.npy' %(crab_name, snle_name), new_mask)

#Plot final image
ax = plt.subplot(1,1,1, projection=WCS(crab.header))
ax.imshow(overlap, cmap='gray', origin='lower', norm='log')
ax.coords['ra'].set_axislabel('Right Ascension')
ax.coords['dec'].set_axislabel('Declination')
plot_name = 'Final_%s_%s.png' %(crab_name, snle_name)
#plt.savefig(plot_name)
plt.close('all')

#Save new files
#outfile_snle_pad = 'Pad_%s.fits' %(snle_name)
#hdu_snle_pad = fits.PrimaryHDU(snle_pad, header=snle.header)
#hdu_snle_pad.writeto(outfile_snle_pad, overwrite=True)

#outfile_snle = 'Reproject_%s.fits' %(snle_name)
#hdu_snle = fits.PrimaryHDU(array, header=crab.header)
#hdu_snle.writeto(outfile_snle, overwrite=True)

outfile = 'Overlap_%s_%s.fits' %(crab_name, snle_name)
hdu = fits.PrimaryHDU(overlap, header=crab.header)
hdu.writeto(outfile, overwrite=True)