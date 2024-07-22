#!/usr/bin/env python3
import argparse
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.colors as colors

from astropy.io.fits import getdata
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize

from photutils.aperture import CircularAperture

# Read in arguments
filename = 'crab1507_20220225_20180212_r_diff.fits'
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fitsfile', '-f', type=str, help='FITS file name', required=False, default=filename)
args = parser.parse_args()

hdulist = fits.open(args.fitsfile)

# Extract image and header
image = hdulist[0].data
header = hdulist[0].header
(ny, nx) = image.shape

x_cen = nx/2.0
y_cen = ny/2.0

# Instantiate WCS object
w = WCS(header)

# Center pixel position on sky
sky_cen = w.pixel_to_world(x_cen, y_cen)

# Sky position of corner pixel
sky_cor = w.pixel_to_world(0, 0)

# Angular distance between center and corner
ang_sep = sky_cen.separation(sky_cor)

# Create a Vizier object containing columns from UCAC4
# rmag < 20 excludes empty values
v = Vizier(columns=["RAJ2000", "DEJ2000", "rmag"], catalog="I/322A/out", column_filters={'rmag': '<20'})

# Lift the nominal 50-row limit
v.ROW_LIMIT = -1

# Query a region around a central position
tl = v.query_region(sky_cen, radius=1.01*ang_sep)

# data contains (RAJ2000 value in degrees, DEJ2000 value in degrees, rmag from APASS)
# rmag is "nan" when it is undefined/unknown

ra_deg = np.array(tl[0]['RAJ2000'])
dec_deg = np.array(tl[0]['DEJ2000'])
r_mag = np.array(tl[0]['rmag'])

# Transform the RA and Dec values to pixel coordinates in the image and only plot the stars that are within the image
# The image is nx by ny pixels
x_pix, y_pix = w.world_to_pixel_values(ra_deg, dec_deg)

r_mag = r_mag[(x_pix >= 0) & (x_pix < nx) & (y_pix >= 0) & (y_pix < ny)]
ra_deg = ra_deg[(x_pix >= 0) & (x_pix < nx) & (y_pix >= 0) & (y_pix < ny)]
dec_deg = dec_deg[(x_pix >= 0) & (x_pix < nx) & (y_pix >= 0) & (y_pix < ny)]

# Create a list of x_pix and y_pix positions where x_pix is greater than or equal to 0 and less than nx and y_pix is greater than or equal to 0 and less than ny
x_pix_in = x_pix[(x_pix >= 0) & (x_pix < nx) & (y_pix >= 0) & (y_pix < ny)]
y_pix_in = y_pix[(x_pix >= 0) & (x_pix < nx) & (y_pix >= 0) & (y_pix < ny)]

print(np.min(r_mag))
print(np.max(r_mag))

#Radius sizes for the circular apertures
r_xlrg = 18.0
r_lrg = 13.0
r_med = 10.0
r_sml = 7.0
r_xsml = 5.0

#Create circular aperture with radius set as a function of r_mag
# Create a list of positions for the circular apertures for r_mag greater than 17.0 with a radius of 5 pixels
positions_17 = list(zip(x_pix_in[r_mag > 17.0], y_pix_in[r_mag > 17.0]))
apertures_17 = CircularAperture(positions_17, r=r_xsml)
apertures_17_masks = apertures_17.to_mask(method='center')
masks_17 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_17_masks]

# Create a list of positions for the circular apertures for r_mag greater than 16.5 and less than or equal to 17.0 with a radius of 5 pixels
positions_165_17 = list(zip(x_pix_in[(r_mag > 16.5) & (r_mag <= 17.0)], y_pix_in[(r_mag > 16.5) & (r_mag <= 17.0)]))
apertures_165_17 = CircularAperture(positions_165_17, r=r_xsml)
apertures_165_17_masks = apertures_165_17.to_mask(method='center')
masks_165_17 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_165_17_masks]

# Create a list of positions for the circular apertures for r_mag greater than 16.2 and less than or equal to 16.5 with a radius of 5 pixels
positions_162_165 = list(zip(x_pix_in[(r_mag > 16.2) & (r_mag <= 16.5)], y_pix_in[(r_mag > 16.2) & (r_mag <= 16.5)]))
apertures_162_165 = CircularAperture(positions_162_165, r=r_xsml)
apertures_162_165_masks = apertures_162_165.to_mask(method='center')
masks_162_165 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_162_165_masks]

# Create a list of positions for the circular apertures for r_mag greater than 16.0 and less than or equal to 16.2 with a radius of 5 pixels
positions_16_162 = list(zip(x_pix_in[(r_mag > 16.0) & (r_mag <= 16.2)], y_pix_in[(r_mag > 16.0) & (r_mag <= 16.2)]))
apertures_16_162 = CircularAperture(positions_16_162, r=r_xsml)
apertures_16_162_masks = apertures_16_162.to_mask(method='center')
masks_16_162 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_16_162_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.9 and less than or equal to 16.0 with a radius of 5 pixels
positions_159_16 = list(zip(x_pix_in[(r_mag > 15.9) & (r_mag <= 16.0)], y_pix_in[(r_mag > 15.9) & (r_mag <= 16.0)]))
apertures_159_16 = CircularAperture(positions_159_16, r=r_xsml)
apertures_159_16_masks = apertures_159_16.to_mask(method='center')
masks_159_16 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_159_16_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.8 and less than or equal to 15.9 with a radius of 5 pixels
positions_158_159 = list(zip(x_pix_in[(r_mag > 15.8) & (r_mag <= 15.9)], y_pix_in[(r_mag > 15.8) & (r_mag <= 15.9)]))
apertures_158_159 = CircularAperture(positions_158_159, r=r_xsml)
apertures_158_159_masks = apertures_158_159.to_mask(method='center')
masks_158_159 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_158_159_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.7 and less than or equal to 15.8 with a radius of 5 pixels
positions_157_158 = list(zip(x_pix_in[(r_mag > 15.7) & (r_mag <= 15.8)], y_pix_in[(r_mag > 15.7) & (r_mag <= 15.8)]))
apertures_157_158 = CircularAperture(positions_157_158, r=r_xsml)
apertures_157_158_masks = apertures_157_158.to_mask(method='center')
masks_157_158 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_157_158_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.6 and less than or equal to 15.7 with a radius of 5 pixels
positions_156_157 = list(zip(x_pix_in[(r_mag > 15.6) & (r_mag <= 15.7)], y_pix_in[(r_mag > 15.6) & (r_mag <= 15.7)]))
apertures_156_157 = CircularAperture(positions_156_157, r=r_xsml)
apertures_156_157_masks = apertures_156_157.to_mask(method='center')
masks_156_157 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_156_157_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.5 and less than or equal to 15.6 with a radius of 5 pixels
positions_155_156 = list(zip(x_pix_in[(r_mag > 15.5) & (r_mag <= 15.6)], y_pix_in[(r_mag > 15.5) & (r_mag <= 15.6)]))
apertures_155_156 = CircularAperture(positions_155_156, r=r_xsml)
apertures_155_156_masks = apertures_155_156.to_mask(method='center')
masks_155_156 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_155_156_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.4 and less than or equal to 15.5 with a radius of 5 pixels
positions_154_155 = list(zip(x_pix_in[(r_mag > 15.4) & (r_mag <= 15.5)], y_pix_in[(r_mag > 15.4) & (r_mag <= 15.5)]))
apertures_154_155 = CircularAperture(positions_154_155, r=r_xsml)
apertures_154_155_masks = apertures_154_155.to_mask(method='center')
masks_154_155 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_154_155_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.3 and less than or equal to 15.4 with a radius of 5 pixels
positions_153_154 = list(zip(x_pix_in[(r_mag > 15.3) & (r_mag <= 15.4)], y_pix_in[(r_mag > 15.3) & (r_mag <= 15.4)]))
apertures_153_154 = CircularAperture(positions_153_154, r=r_xsml)
apertures_153_154_masks = apertures_153_154.to_mask(method='center')
masks_153_154 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_153_154_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.2 and less than or equal to 15.3 with a radius of 5 pixels
positions_152_153 = list(zip(x_pix_in[(r_mag > 15.2) & (r_mag <= 15.3)], y_pix_in[(r_mag > 15.2) & (r_mag <= 15.3)]))
apertures_152_153 = CircularAperture(positions_152_153, r=r_xsml)
apertures_152_153_masks = apertures_152_153.to_mask(method='center')
masks_152_153 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_152_153_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.1 and less than or equal to 15.2 with a radius of 5 pixels
positions_151_152 = list(zip(x_pix_in[(r_mag > 15.1) & (r_mag <= 15.2)], y_pix_in[(r_mag > 15.1) & (r_mag <= 15.2)]))
apertures_151_152 = CircularAperture(positions_151_152, r=r_xsml)
apertures_151_152_masks = apertures_151_152.to_mask(method='center')
masks_151_152 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_151_152_masks]

# Create a list of positions for the circular apertures for r_mag greater than 15.0 and less than or equal to 15.1 with a radius of 5 pixels
positions_15_151 = list(zip(x_pix_in[(r_mag > 15.0) & (r_mag <= 15.1)], y_pix_in[(r_mag > 15.0) & (r_mag <= 15.1)]))
apertures_15_151 = CircularAperture(positions_15_151, r=r_xsml)
apertures_15_151_masks = apertures_15_151.to_mask(method='center')
masks_15_151 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_15_151_masks]

           
# Create a list of positions for the circular apertures for r_mag greater than 14.9 and less than or equal to 15.0 with a radius of 7 pixels
positions_149_15 = list(zip(x_pix_in[(r_mag > 14.9) & (r_mag <= 15.0)], y_pix_in[(r_mag > 14.9) & (r_mag <= 15.0)]))
apertures_149_15 = CircularAperture(positions_149_15, r=r_sml)
apertures_149_15_masks = apertures_149_15.to_mask(method='center')
masks_149_15 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_149_15_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.8 and less than or equal to 14.9 with a radius of 7 pixels
positions_148_149 = list(zip(x_pix_in[(r_mag > 14.8) & (r_mag <= 14.9)], y_pix_in[(r_mag > 14.8) & (r_mag <= 14.9)]))
apertures_148_149 = CircularAperture(positions_148_149, r=r_sml)
apertures_148_149_masks = apertures_148_149.to_mask(method='center')
masks_148_149 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_148_149_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.7 and less than or equal to 14.8 with a radius of 7 pixels
positions_147_148 = list(zip(x_pix_in[(r_mag > 14.7) & (r_mag <= 14.8)], y_pix_in[(r_mag > 14.7) & (r_mag <= 14.8)]))
apertures_147_148 = CircularAperture(positions_147_148, r=r_sml)
apertures_147_148_masks = apertures_147_148.to_mask(method='center')
masks_147_148 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_147_148_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.6 and less than or equal to 14.7 with a radius of 7 pixels
positions_146_147 = list(zip(x_pix_in[(r_mag > 14.6) & (r_mag <= 14.7)], y_pix_in[(r_mag > 14.6) & (r_mag <= 14.7)]))
apertures_146_147 = CircularAperture(positions_146_147, r=r_sml)
apertures_146_147_masks = apertures_146_147.to_mask(method='center')
masks_146_147 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_146_147_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.4 and less than or equal to 14.6 with a radius of 7 pixels
positions_144_146 = list(zip(x_pix_in[(r_mag > 14.4) & (r_mag <= 14.6)], y_pix_in[(r_mag > 14.4) & (r_mag <= 14.6)]))
apertures_144_146 = CircularAperture(positions_144_146, r=r_sml)
apertures_144_146_masks = apertures_144_146.to_mask(method='center')
masks_144_146 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_144_146_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.2 and less than or equal to 14.4 with a radius of 7 pixels
positions_142_144 = list(zip(x_pix_in[(r_mag > 14.2) & (r_mag <= 14.4)], y_pix_in[(r_mag > 14.2) & (r_mag <= 14.4)]))
apertures_142_144 = CircularAperture(positions_142_144, r=r_sml)
apertures_142_144_masks = apertures_142_144.to_mask(method='center')
masks_142_144 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_142_144_masks]

# Create a list of positions for the circular apertures for r_mag greater than 14.0 and less than or equal to 14.2 with a radius of 7 pixels
positions_14_142 = list(zip(x_pix_in[(r_mag > 14.0) & (r_mag <= 14.2)], y_pix_in[(r_mag > 14.0) & (r_mag <= 14.2)]))
apertures_14_142 = CircularAperture(positions_14_142, r=r_sml)
apertures_14_142_masks = apertures_14_142.to_mask(method='center')
masks_14_142 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_14_142_masks]

# Create a list of positions for the circular apertures for r_mag greater than 13.9 and less than or equal to 14.0 with a radius of 7 pixels
positions_139_14 = list(zip(x_pix_in[(r_mag > 13.9) & (r_mag <= 14.0)], y_pix_in[(r_mag > 13.9) & (r_mag <= 14.0)]))
apertures_139_14 = CircularAperture(positions_139_14, r=r_sml)
apertures_139_14_masks = apertures_139_14.to_mask(method='center')
masks_139_14 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_139_14_masks]

# Create a list of positions for the circular apertures for r_mag greater than 13.7 and less than or equal to 13.9 with a radius of 7 pixels
positions_137_139 = list(zip(x_pix_in[(r_mag > 13.7) & (r_mag <= 13.9)], y_pix_in[(r_mag > 13.7) & (r_mag <= 13.9)]))
apertures_137_139 = CircularAperture(positions_137_139, r=r_sml)
apertures_137_139_masks = apertures_137_139.to_mask(method='center')
masks_137_139 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_137_139_masks]

# Create a list of positions for the circular apertures for r_mag greater than 13.3 and less than or equal to 13.7 with a radius of 7 pixels
positions_133_137 = list(zip(x_pix_in[(r_mag > 13.3) & (r_mag <= 13.7)], y_pix_in[(r_mag > 13.3) & (r_mag <= 13.7)]))
apertures_133_137 = CircularAperture(positions_133_137, r=r_sml)
apertures_133_137_masks = apertures_133_137.to_mask(method='center')
masks_133_137 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_133_137_masks]

# Create a list of positions for the circular apertures for r_mag greater than 13.0 and less than or equal to 13.3 with a radius of 7 pixels
positions_13_133 = list(zip(x_pix_in[(r_mag > 13.0) & (r_mag <= 13.3)], y_pix_in[(r_mag > 13.0) & (r_mag <= 13.3)]))
apertures_13_133 = CircularAperture(positions_13_133, r=r_sml)
apertures_13_133_masks = apertures_13_133.to_mask(method='center')
masks_13_133 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_13_133_masks]

# Create a list of positions for the circular apertures for r_mag greater than 12.5 and less than or equal to 13.0 with a radius of 10 pixels
positions_125_13 = list(zip(x_pix_in[(r_mag > 12.5) & (r_mag <= 13.0)], y_pix_in[(r_mag > 12.5) & (r_mag <= 13.0)]))
apertures_125_13 = CircularAperture(positions_125_13, r=r_med)
apertures_125_13_masks = apertures_125_13.to_mask(method='center')
masks_125_13 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_125_13_masks]

# Create a list of positions for the circular apertures for r_mag greater than 12.0 and less than or equal to 12.5 with a radius of 10 pixels
positions_12_125 = list(zip(x_pix_in[(r_mag > 12.0) & (r_mag <= 12.5)], y_pix_in[(r_mag > 12.0) & (r_mag <= 12.5)]))
apertures_12_125 = CircularAperture(positions_12_125, r=r_med)
apertures_12_125_masks = apertures_12_125.to_mask(method='center')
masks_12_125 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_12_125_masks]

# Create a list of positions for the circular apertures for r_mag greater than 11.0 and less than or equal to 12.0 with a radius of 10 pixels
positions_11_12 = list(zip(x_pix_in[(r_mag > 11.0) & (r_mag <= 12.0)], y_pix_in[(r_mag > 11.0) & (r_mag <= 12.0)]))
apertures_11_12 = CircularAperture(positions_11_12, r=r_med)
apertures_11_12_masks = apertures_11_12.to_mask(method='center')
masks_11_12 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_11_12_masks]

# Create a list of positions for the circular apertures for r_mag greater than 10.0 and less than or equal to 11.0 with a radius of 13 pixels
positions_10_11 = list(zip(x_pix_in[(r_mag > 10.0) & (r_mag <= 11.0)], y_pix_in[(r_mag > 10.0) & (r_mag <= 11.0)]))
apertures_10_11 = CircularAperture(positions_10_11, r=r_lrg)
apertures_10_11_masks = apertures_10_11.to_mask(method='center')
masks_10_11 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_10_11_masks]

# Create a list of positions for the circular apertures for r_mag less than or equal to 10.0 with a radius of 18 pixels
positions_10 = list(zip(x_pix_in[r_mag <= 10.0], y_pix_in[r_mag <= 10.0]))
apertures_10 = CircularAperture(positions_10, r=r_xlrg)
apertures_10_masks = apertures_10.to_mask(method='center')
masks_10 = [aper_ma.to_image((ny, nx)) for aper_ma in apertures_10_masks]

#Mask out bright stars with circular apertures
image = getdata(args.fitsfile)
image[np.logical_or.reduce(masks_17)] = 0.0
print("greater than 17.0")
image[np.logical_or.reduce(masks_165_17)] = 0.0
print("greater than 16.5 and less than or equal to 17.0")
image[np.logical_or.reduce(masks_162_165)] = 0.0
print("greater than 16.2 and less than or equal to 16.5")
image[np.logical_or.reduce(masks_16_162)] = 0.0
print("greater than 16.0 and less than or equal to 16.2")
image[np.logical_or.reduce(masks_159_16)] = 0.0
print("greater than 15.9 and less than or equal to 16.0")
image[np.logical_or.reduce(masks_158_159)] = 0.0
print("greater than 15.8 and less than or equal to 15.9")
image[np.logical_or.reduce(masks_157_158)] = 0.0
print("greater than 15.7 and less than or equal to 15.8")
image[np.logical_or.reduce(masks_156_157)] = 0.0
print("greater than 15.6 and less than or equal to 15.7")
image[np.logical_or.reduce(masks_155_156)] = 0.0
print("greater than 15.5 and less than or equal to 15.6")
image[np.logical_or.reduce(masks_154_155)] = 0.0
print("greater than 15.4 and less than or equal to 15.5")
image[np.logical_or.reduce(masks_153_154)] = 0.0
print("greater than 15.3 and less than or equal to 15.4")
image[np.logical_or.reduce(masks_152_153)] = 0.0
print("greater than 15.2 and less than or equal to 15.3")
image[np.logical_or.reduce(masks_151_152)] = 0.0
print("greater than 15.1 and less than or equal to 15.2")
image[np.logical_or.reduce(masks_15_151)] = 0.0
print("greater than 15.0 and less than or equal to 15.1")
image[np.logical_or.reduce(masks_149_15)] = 0.0
print("greater than 14.9 and less than or equal to 15.0")
image[np.logical_or.reduce(masks_148_149)] = 0.0
print("greater than 14.8 and less than or equal to 14.9")
image[np.logical_or.reduce(masks_147_148)] = 0.0
print("greater than 14.7 and less than or equal to 14.8")
image[np.logical_or.reduce(masks_146_147)] = 0.0
print("greater than 14.6 and less than or equal to 14.7")
image[np.logical_or.reduce(masks_144_146)] = 0.0
print("greater than 14.4 and less than or equal to 14.6")
image[np.logical_or.reduce(masks_142_144)] = 0.0
print("greater than 14.2 and less than or equal to 14.4")
image[np.logical_or.reduce(masks_14_142)] = 0.0
print("greater than 14.0 and less than or equal to 14.2")
image[np.logical_or.reduce(masks_139_14)] = 0.0
print("greater than 13.9 and less than or equal to 14.0")
image[np.logical_or.reduce(masks_137_139)] = 0.0
print("greater than 13.7 and less than or equal to 13.9")
image[np.logical_or.reduce(masks_133_137)] = 0.0
print("greater than 13.3 and less than or equal to 13.7")
image[np.logical_or.reduce(masks_13_133)] = 0.0
print("greater than 13.0 and less than or equal to 13.3")
image[np.logical_or.reduce(masks_125_13)] = 0.0
print("greater than 12.5 and less than or equal to 13.0")
image[np.logical_or.reduce(masks_12_125)] = 0.0
print("greater than 12.0 and less than or equal to 12.5")
image[np.logical_or.reduce(masks_11_12)] = 0.0
print("greater than 11.0 and less than or equal to 12.0")
image[np.logical_or.reduce(masks_10_11)] = 0.0
print("greater than 10.0 and less than or equal to 11.0")
image[np.logical_or.reduce(masks_10)] = 0.0
print("less than or equal to 10.0")

# Plot the image with the stars masked out
plt.figure(figsize=(12,8))
ax = plt.subplot(projection=w)
norm = ImageNormalize(image, interval=ZScaleInterval(),
                      stretch=LinearStretch())
ax.imshow(image, norm=norm, interpolation='none', origin='lower', cmap='gray')
#ax.scatter(ra_deg, dec_deg, transform=ax.get_transform('fk5'), s=10, color='red')
ax.grid()
ax.set_xlabel('RA J2000.0')
ax.set_ylabel('Dec J2000.0')
plt.show()

#Save new image
name = filename.split('.')[0]
original_image = fits.open(filename)[0]
outfile = '%s_m.fits' %(name)
hdu = fits.PrimaryHDU(image, header=original_image.header)
hdu.writeto(outfile, overwrite=True)