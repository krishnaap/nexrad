#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 14:52:09 2023

@author: krishna
"""


#%%
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyart
import cartopy.crs as ccrs
# Set global font size
mpl.rcParams.update({'font.size': 20})  # Adjust the size as needed
# Path to your NEXRAD radar file
nexrad_file_path = 'KHGX_SDUS34_N1PHGX_202207221910'
# Read the radar file using Py-ART
radar = pyart.io.read(nexrad_file_path)
# Set up a map projection
projection = ccrs.PlateCarree()
# Create a plot
fig = plt.figure(figsize=(12, 10))  # Adjust figure size as needed
ax = plt.axes(projection=projection)
# Plot the radar data
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi_map('radar_estimated_rain_rate', 0, ax=ax, projection=projection,
                     vmin=0, vmax=2, cmap='pyart_HomeyerRainbow')
# Add coastlines
ax.coastlines('10m')
# Customize gridlines
gl = ax.gridlines(draw_labels=True)
gl.top_labels = False  # Disable top x-axis labels
gl.right_labels = False  # Disable right y-axis labels

# Show the plot with larger title
plt.title('NEXRAD Radar Estimated Rain Rate with Marked Star', fontsize=16)  # Adjust title font size as needed
plt.show()

