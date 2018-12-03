
"""
From <Q>, <I>, <U/I>, <Q/I> finds Thomson scattergin using Ne according to figure 5 in:
http://iopscience.iop.org/article/10.1088/2041-8205/786/2/L19/meta#apjl493681f5

Author: Juan Camilo Guevara Gomez
Latest version: 26 October 2018
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
warnings.filterwarnings('ignore')

#Data for basic calculations, extract shape, etc
datatest = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_123600_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
#datatest = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_160000_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
#Sun radius, arcsec/pixel and central positions taken from header
rsun = 953.105652
arpi = 0.5016345536842105
xc = 2048
yc = 2048
inclang = -8.81220944
#Region over the event, position and large in pixels
xpos = 60 #bottom left rectangle
ypos = 2320 # bottom left rectangle
bsz = 40 #height and width rectangle
#rectangle or region widht an inclination angle of 10 to aligned with loop
r = matplotlib.patches.Rectangle((xpos,ypos),bsz,bsz,angle=inclang,linewidth=1,edgecolor='r',facecolor='none')
#Calculation of Ne for each position according to figure 5 from SHP
#Ne = (I_sc/I_o(ref))(1/thomson)
y_int = np.arange(ypos-30, ypos+bsz+30)
x_int = np.arange(xpos-30, xpos+bsz+30)
g = np.meshgrid(x_int, y_int)
coords = list(zip(*(c.flat for c in g))) #virtual grid for getting coordinates of the rectangle r
del y_int,x_int,g
enc_coor = np.vstack([p for p in coords if r.contains_point(p, radius=0)])
distances_rsun = [] ##Distance from region point to limb in solar radii
for ind,item in enumerate(enc_coor):
    x2 = item[0]
    y2 = item[1]
    h =(math.hypot(x2-xc,y2-yc)-rsun)/rsun
    distances_rsun.append(h)
    del x2,y2,h
psh_points = np.genfromtxt('hmi_055_psh.csv', delimiter=',') #Points from figure 5
th_step1 = np.interp(distances_rsun,psh_points[:,0],psh_points[:,1]) #interpolation using psh_points
grid_Th = np.empty(np.shape(datatest)) * np.nan #same size as data but empty, it will Fill with thomson values
for idx,item in enumerate(enc_coor):
    grid_Th[item[1],item[0]] = th_step1[idx]
refw = int(bsz/2)
I_0 = np.nanmean(datatest[yc-refw:yc+refw,xc-refw:xc+refw])
# I_0 = np.nanmean(datatest[600-refw:600+refw,700-refw:700+refw])
#I_0 = 160
grid_Isc = np.empty(np.shape(datatest)) * np.nan #same sizes as data but empty, it will fill with I_sc values
for idx,item in enumerate(enc_coor):
    grid_Isc[item[1],item[0]]=datatest[item[1],item[0]]
grid_Ne = ((grid_Isc/200)/I_0)*(1/grid_Th)
grid_ne = grid_Ne/(0.725*100000000)

# plt.imshow(grid_Ne,origin='lower',cmap='plasma',interpolation='hamming',vmin=1e20,vmax=2.5e20)
# plt.colorbar()
# plt.xlim(57,110)
# plt.ylim(2310,2363)
# plt.title(r'Ne according to Figure 5 [cm$^{-2}$] at 2017-09-10T12:36:56')
# plt.xlabel(r'Position [pixels]')
# plt.ylabel('Position [pixels]')
# plt.show()
#
# plt.imshow(grid_ne,origin='lower',cmap='plasma',interpolation='hamming')
# plt.colorbar()
# plt.xlim(57,110)
# plt.ylim(2310,2363)
# plt.title(r'Electron density map [cm$^{-3}$] at 2017-09-10T12:36:56')
# plt.xlabel(r'Position [pixels]')
# plt.ylabel('Position [pixels]')
# plt.show()
#
# flat = []
# for item in grid_ne:
#     for jitem in item:
#         if ~np.isnan(jitem):
#             flat.append(jitem)
#
# b,bins,_= plt.hist(flat,60,facecolor='green',alpha=0.75)
# plt.title(r'Electron density histogram distribution at 2017-09-10T12:36:56')
# plt.xlabel(r'n$_e$'+r'[cm$^{-3}$]')
# plt.ylabel(r'Counts inside region')
# plt.show()
