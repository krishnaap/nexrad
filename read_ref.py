#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 11:48:32 2023

@author: krishna
"""
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import numpy as np

# Set global font size
mpl.rcParams.update({'font.size': 20})  # Adjust the size as needed

# Path to your NetCDF file
netcdf_file_path = '/media/nexrad/8_13/NEXRAD-selected/KHGX20220813_192627_V06.nc'

# Open the NetCDF file
dataset = nc.Dataset(netcdf_file_path)

# Extract latitude and longitude as scalar values
latitude = dataset.variables['latitude'][0]
longitude = dataset.variables['longitude'][0]

# Extract the reflectivity data for the first elevation level
ref_data = dataset.variables['REF'][0, :, :]  # Adjust indices as needed

# Set up a map projection
projection = ccrs.PlateCarree()

# Create a figure
fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': projection})  # Adjust figure size as needed

# Add features to the map
ax.coastlines()
gl = ax.gridlines(draw_labels=True)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels

# Set the extent of the map to the specified lat-lon bounds
# ax.set_extent([-95.75, -95, 29.5, 30.25], crs=ccrs.PlateCarree())

# Plotting the reflectivity data
reflectivity = ax.imshow(ref_data, origin='lower', extent=[longitude - 1, longitude + 1, latitude - 1, latitude + 1], transform=projection, interpolation='nearest', cmap='jet')

# Add a color bar
cbar = plt.colorbar(reflectivity, ax=ax, orientation='vertical', pad=0.05, aspect=50)
cbar.set_label('Reflectivity (dBZ)')

# Mark a specific point
specific_lat = 29.47190094
specific_lon = -95.07873535 
ax.plot(specific_lon, specific_lat, 'ro', markersize=20, transform=ccrs.Geodetic())

# Show the plot with larger title
plt.title('Radar Reflectivity', fontsize=16)  # Adjust title font size as needed
plt.show()

