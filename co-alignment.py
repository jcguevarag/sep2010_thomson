
"""
Co-alignment of SDO/HMI, RHESSI and Ne results data using as inputs:
Ne_v2.npy ----------------------------------------------------------------------------> Densities [cm^2] calculated
(Crops [2260:2410,25:175] in 4096x4096 arrays)
tiempos_v2.npy -----------------------------------------------------------------------> Times from HMI data in datetime format
hsi_image_20170910_161510.fits -------------------------------------------------------> RHESSI image 1
hsi_image_20170910_164110.fits -------------------------------------------------------> RHESSI image 2
hmi.S_90s.20170910_161200_TAI.3I_recentered.medians_4pics.meanvalue.fits -------------> HMI image 1
hmi.S_90s.20170910_164200_TAI.3I_recentered.medians_4pics.meanvalue.fits -------------> HMI image 2

It has been taking into account a roll about 0.21 degrees clockwise about the solar center between RHESSI and HMI according to:
https://arxiv.org/pdf/1511.04161.pdf


Author: Juan Camilo Guevara Gomez
Latest version: 19 December 2018
"""
import matplotlib
import matplotlib.pyplot as plt
import urllib.request
import numpy as np
import math
import peakutils
import math
import time,os,sys
import datetime
import matplotlib.dates as mdates
from os import listdir
from os.path import isfile, join
import numpy.ma as ma
import sunpy
import sunpy.map
from sunpy.timeseries import TimeSeries
from sunpy.time import TimeRange, parse_time
from sunpy.net import hek, Fido, attrs as a
import astropy.units as u
from astropy.coordinates import SkyCoord
import pickle
import warnings
from scipy import stats
import copy
warnings.filterwarnings('ignore')
np.seterr(divide='ignore', invalid='ignore')

#Read densities and time data

Ne = np.load('Ne_v5.npy')
Ne[Ne<0]=np.nan
times = np.load('tiempos_v2.npy')


for idx,item in enumerate(times):
	if item.hour == 16:
		if item.minute == 12:
			Ne_1 = Ne[idx]
			T_1 = datetime.datetime.strftime(item,'%Y-%m-%d %H:%M:%S')
			print('idx_1',idx)
		if item.minute == 42:
			Ne_2 = Ne[idx]
			T_2 = datetime.datetime.strftime(item,'%Y-%m-%d %H:%M:%S')
			print('idx_2',idx)

# # Ne_1 = np.nan_to_num(Ne_1,copy=True)
# Ne_1[np.isnan(Ne_1)]=np.nanmean(Ne_1)-np.nanstd(Ne_1)
# Ne_2[np.isnan(Ne_2)]=np.nanmean(Ne_2)-np.nanstd(Ne_2)
#Read RHESSI data

rh_1 = sunpy.map.Map('hsi_image_20170910_161510.fits')
rh_2 = sunpy.map.Map('hsi_image_20170910_164110.fits')
# rh_1 = rh_1.rotate(+0.21*u.deg)
# rh_2 = rh_2.rotate(+0.21*u.deg)

#Read HMI data

hmi_1_data = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_161200_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
hmi_1_meta = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_161200_TAI.3I_recentered.medians_4pics.meanvalue.fits').meta
hmi_1_data = np.clip(hmi_1_data, 250, 900)
hmi_1 = sunpy.map.Map(hmi_1_data,hmi_1_meta)
hmi_1 = hmi_1.rotate(0.21*u.deg)
hmi_1.plot_settings['cmap'] = plt.get_cmap('rainbow')

hmi_2_data = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_164200_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
hmi_2_meta = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_164200_TAI.3I_recentered.medians_4pics.meanvalue.fits').meta
hmi_2_data = np.clip(hmi_2_data, 250, 900)
hmi_2 = sunpy.map.Map(hmi_2_data,hmi_2_meta)
hmi_2 = hmi_2.rotate(0.21*u.deg)
hmi_2.plot_settings['cmap'] = plt.get_cmap('rainbow')





#Building up headers for densities

Ne_1_hdr = hmi_1.meta.copy()
Ne_2_hdr = hmi_2.meta.copy()

Ne_1_hdr['naxis1'] = Ne_1.shape[0]
Ne_1_hdr['naxis2'] = Ne_1.shape[1]


x_1 = float(hmi_1.meta['crpix1']) - 25  ###numbers taken from cropping pixels 25 and 2260, check thomson_scattering_recursive_v2.py
y_1 = 2260-float(hmi_1.meta['crpix2'])
Ne_1_hdr['crpix1'] = x_1
Ne_1_hdr['crpix2'] = -y_1

x_2 = float(hmi_2.meta['crpix1']) - 25
y_2 = 2260-float(hmi_2.meta['crpix2'])
Ne_2_hdr['crpix1'] = x_2
Ne_2_hdr['crpix2'] = -y_2

Ne_map1 = sunpy.map.Map(Ne_1,Ne_1_hdr)
Ne_map2 = sunpy.map.Map(Ne_2,Ne_2_hdr)

maps_1 = sunpy.map.Map(hmi_1,Ne_map1,rh_1,composite=True)
maps_1.set_levels(2,[50,70,90,95],percent=True)
maps_1.set_levels(1,[10,30,50,70,90,95],percent=True)
maps_1.plot(title='%s'%T_1),plt.show()

maps_2 = sunpy.map.Map(hmi_2,Ne_map2,rh_2,composite=True)
maps_2.set_levels(2,[40,50,70,90,95],percent=True)
maps_2.set_levels(1,[10,40,50,70,90,95],percent=True)
maps_2.plot(title='%s'%T_2),plt.show()
