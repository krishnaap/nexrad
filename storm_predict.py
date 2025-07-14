#!/usr/bin/env python3
"""Estimate storm movement using ML and GFS wind data."""

from datetime import datetime, timedelta

import cv2
import numpy as np
import pyart
import xarray as xr

# Example sample path for NEXRAD files
DEFAULT_NEXRAD_PATH = 'KHGX_SDUS34_N1PHGX_202207221910'


# ----------------- Radar utilities -----------------

def read_rain_rate(path):
    """Return rain rate array and radar object for a NEXRAD file."""
    radar = pyart.io.read(path)
    field = radar.fields['radar_estimated_rain_rate']['data']
    return field.filled(np.nan), radar


def compute_optical_flow(prev, curr):
    """Compute mean optical flow vector between two 2-D arrays."""
    prev8 = np.nan_to_num(prev).astype(np.float32)
    curr8 = np.nan_to_num(curr).astype(np.float32)
    flow = cv2.calcOpticalFlowFarneback(
        prev8, curr8, None, 0.5, 3, 15, 3, 5, 1.2, 0
    )
    u = np.nanmean(flow[..., 0])
    v = np.nanmean(flow[..., 1])
    return u, v


# ----------------- GFS utilities -----------------


def fetch_gfs_wind(lat, lon, dt):
    """Fetch 850 hPa wind from the GFS open dataset (0.25 deg)."""
    date = dt.strftime('%Y%m%d')
    cycle = f"{dt.hour // 6 * 6:02d}"
    url = (
        f"https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs{date}/" f"gfs_0p25_{cycle}z"
    )
    ds = xr.open_dataset(url)
    ds_sel = ds.sel(time=dt, method='nearest')
    u = ds_sel['ugrdprs'].sel(isobaric1=85000, lat=lat, lon=lon, method='nearest').item()
    v = ds_sel['vgrdprs'].sel(isobaric1=85000, lat=lat, lon=lon, method='nearest').item()
    return float(u), float(v)


# ----------------- Main execution -----------------

if __name__ == '__main__':
    # Example usage with two sequential radar files
    import argparse

    parser = argparse.ArgumentParser(description="Estimate storm motion")
    parser.add_argument('prev_file', nargs='?', default=DEFAULT_NEXRAD_PATH,
                        help='older NEXRAD file')
    parser.add_argument('curr_file', nargs='?', default=DEFAULT_NEXRAD_PATH,
                        help='newer NEXRAD file')
    args = parser.parse_args()

    prev_field, radar_prev = read_rain_rate(args.prev_file)
    curr_field, radar_curr = read_rain_rate(args.curr_file)

    u_pix, v_pix = compute_optical_flow(prev_field, curr_field)

    # Convert pixel displacement to km using radar gate spacing (approx.)
    gate_length = radar_curr.range['data'][1] / 1000.0  # km
    dx_km = u_pix * gate_length
    dy_km = v_pix * gate_length

    time_curr = datetime.strptime(radar_curr.time['units'].split('since')[1].strip(), '%Y-%m-%dT%H:%M:%SZ') + timedelta(seconds=float(radar_curr.time['data'][0]))
    lat = float(radar_curr.latitude['data'][0])
    lon = float(radar_curr.longitude['data'][0])

    print(f"Optical flow storm motion: dx={dx_km:.2f} km, dy={dy_km:.2f} km")

    try:
        ugfs, vgfs = fetch_gfs_wind(lat, lon, time_curr)
        print(f"GFS 850 hPa wind: u={ugfs:.2f} m/s, v={vgfs:.2f} m/s")
    except Exception as exc:
        print(f"Unable to download GFS data: {exc}")

