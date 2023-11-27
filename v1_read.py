#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 11:48:32 2023

@author: krishna
"""


import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geopy.distance import geodesic
# Path to your NetCDF file
netcdf_file_path = '/media/krishna/Media1/tmp/nexrad/KHGX20220722_232535_V06.nc'
# Open the NetCDF file
dataset = nc.Dataset(netcdf_file_path)

# Extract latitude and longitude as scalar values
latitude = dataset.variables['latitude'][0]
longitude = dataset.variables['longitude'][0]

# Extract the reflectivity data for the first elevation level
ref_data = dataset.variables['REF'][0, :, :]  # Adjust indices as needed

# Set up a map projection
projection = ccrs.AzimuthalEquidistant(central_longitude=longitude, central_latitude=latitude)

# Create a figure
plt.figure(figsize=(10, 8))
ax = plt.axes(projection=projection)

# Add features to the map
ax.coastlines()
ax.gridlines()

# Calculate the extent using geodesic distances
km_in_degrees_lat = geodesic((latitude, longitude), (latitude + 1, longitude)).km
km_in_degrees_lon = geodesic((latitude, longitude), (latitude, longitude + 1)).km

extent_size_km = 250  # 250 km from the radar
delta_lat = extent_size_km / km_in_degrees_lat
delta_lon = extent_size_km / km_in_degrees_lon

ax.set_extent([longitude - delta_lon, longitude + delta_lon, latitude - delta_lat, latitude + delta_lat], crs=ccrs.PlateCarree())

# Plotting the data
ax.imshow(ref_data, origin='lower', extent=[longitude - delta_lon, longitude + delta_lon, latitude - delta_lat, latitude + delta_lat], transform=ccrs.PlateCarree(), interpolation='nearest')

plt.title('Radar Reflectivity')
plt.show()
#%%

import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geopy.distance import geodesic
import numpy as np

# Path to your NetCDF file
netcdf_file_path = '/media/krishna/Media1/tmp/nexrad/KHGX20220722_232535_V06.nc'

# Open the NetCDF file
dataset = nc.Dataset(netcdf_file_path)

# Extract latitude and longitude as scalar values
latitude = dataset.variables['latitude'][0]
longitude = dataset.variables['longitude'][0]

# Extract the reflectivity data for the first elevation level
ref_data = dataset.variables['REF'][0, :, :]  # Adjust indices as needed

# Grid size and resolution (you may need to adjust these)
grid_size_km = 250  # Half the size of the grid in kilometers
grid_resolution_km = 1  # Distance between grid points in kilometers

# Calculate the extent using geodesic distances
km_in_degrees_lat = geodesic((latitude, longitude), (latitude + 1, longitude)).km
km_in_degrees_lon = geodesic((latitude, longitude), (latitude, longitude + 1)).km

delta_lat = grid_size_km / km_in_degrees_lat
delta_lon = grid_size_km / km_in_degrees_lon

extent = [longitude - delta_lon, longitude + delta_lon, latitude - delta_lat, latitude + delta_lat]

# Set up a map projection
projection = ccrs.PlateCarree()

# Create a figure
plt.figure(figsize=(10, 8))
ax = plt.axes(projection=projection)

# Add features to the map
ax.coastlines()
ax.gridlines()

# Set the extent of the plot
ax.set_extent(extent, crs=ccrs.PlateCarree())

# Calculate lat/lon grid
x = np.linspace(extent[0], extent[1], ref_data.shape[1])
y = np.linspace(extent[2], extent[3], ref_data.shape[0])

# Plotting the data
ax.imshow(ref_data, origin='lower', extent=extent, transform=ccrs.PlateCarree(), interpolation='nearest', cmap='viridis')

# Set labels for the axes
ax.set_xticks(x[::50], crs=ccrs.PlateCarree())  # Modify the step as needed
ax.set_yticks(y[::50], crs=ccrs.PlateCarree())  # Modify the step as needed
ax.xaxis.set_major_formatter(ccrs.PlateCarree().xaxis.get_major_formatter())
ax.yaxis.set_major_formatter(ccrs.PlateCarree().yaxis.get_major_formatter())

plt.title('Radar Reflectivity')
plt.show()
#%%
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geopy.distance import geodesic
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Path to your NetCDF file
netcdf_file_path = '/media/krishna/Media1/tmp/nexrad/KHGX20220722_190443_V06.nc'

# Open the NetCDF file
dataset = nc.Dataset(netcdf_file_path)

# Extract latitude and longitude as scalar values
latitude = dataset.variables['latitude'][0]
longitude = dataset.variables['longitude'][0]

# Extract the reflectivity data for the first elevation level
ref_data = dataset.variables['REF'][0, :, :]  # Adjust indices as needed

# Calculate the extent using geodesic distances
km_in_degrees_lat = geodesic((latitude, longitude), (latitude + 1, longitude)).km
km_in_degrees_lon = geodesic((latitude, longitude), (latitude, longitude + 1)).km

delta_lat = 200 / km_in_degrees_lat
delta_lon = 200 / km_in_degrees_lon

extent = [longitude - delta_lon, longitude + delta_lon, latitude - delta_lat, latitude + delta_lat]

# Set up a map projection
projection = ccrs.PlateCarree()

# Create a figure
plt.figure(figsize=(10, 8))
ax = plt.axes(projection=projection)

# Mark a specific point
specific_lat = 29.47190094
specific_lon = -95.07873535 
ax.plot(specific_lon, specific_lat, 'ro', markersize=10, transform=ccrs.Geodetic())

# Mark the specific point with a star
specific_lat = 29.90172
specific_lon = -95.32496
ax.plot(specific_lon, specific_lat, marker='*', color='black', markersize=12, transform=ccrs.Geodetic())


# Add features to the map
ax.coastlines()
ax.gridlines()

# Set the extent of the plot
ax.set_extent(extent, crs=ccrs.PlateCarree())

# Plotting the data
ax.imshow(ref_data, origin='lower', extent=extent, transform=ccrs.PlateCarree(), interpolation='nearest', cmap='viridis')

# Set labels for the axes
xticks = np.linspace(extent[0], extent[1], 6)
yticks = np.linspace(extent[2], extent[3], 6)
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())

plt.title('Radar Reflectivity')
plt.show()
