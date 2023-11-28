# NEXRAD plotting
Python codes to plot NEXRAD radar reflectivity and precipitation over a map. A brief note about the Next Generation Weather Radar (NEXRAD) system is a network of 160 high-resolution S-band Doppler weather radars jointly operated by the National Weather Service (NWS), the Federal Aviation Administration (FAA), and the U.S. Air Force. The NEXRAD system detects precipitation and wind, and its data can be processed to map precipitation patterns and movement(https://www.ncei.noaa.gov/products/radar/next-generation-weather-radar). There are different sources available to download, and one of the best sources to download data is from the NEXRAD inventory (https://www.ncdc.noaa.gov/nexradinv/). Here I am using two datasets, one for radar reflectivity and the second for precipitation. 

## Requirements
netCDF4
matplotlib
cartopy
geopy
NumPy

Additionally, in the precipitation plotting, I used Pyart (https://github.com/ARM-DOE/pyart)
