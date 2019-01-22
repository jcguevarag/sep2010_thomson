
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
from scipy import stats
warnings.filterwarnings('ignore')

#Data for basic calculations, extract shape, etc
# datatest = sunpy.map.Map('Stokes_IQUV_averages/hmi.S_90s.20170910_123600_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
datatest = sunpy.map.Map('/Users/juancg/Documents/evento_sep10/Stokes_IQUV_averages/hmi.S_90s.20170910_161800_TAI.3I_recentered.medians_4pics.meanvalue.fits').data
#Sun radius, arcsec/pixel and central positions taken from header
rsun = 953.105652
arpi = 0.5016345536842105
xc = 2048
yc = 2048
inclang = -8.81220944
#Region over the event, position and large in pixels
xpos = 60 #bottom left rectangle
ypos = 2310 # bottom left rectangle
bsz = 50 #height and width rectangle


# xpos = 30 #bottom left rectangle
# ypos = 2200 # bottom left rectangle
# bsz = 300 #height and width rectangle

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
psh_points = np.genfromtxt('/Users/juancg/Documents/evento_sep10/sep2010_thomson/hmi_055_psh.csv', delimiter=',') #Points from figure 5 top
psh_points_pol = np.genfromtxt('/Users/juancg/Documents/evento_sep10/sep2010_thomson/polarization_v_055_blob.csv', delimiter=',') #Points from figure 5 bottom
th_step1 = np.interp(distances_rsun,psh_points[:,0],psh_points[:,1]) #interpolation using psh_points
pol_step1 = np.interp(distances_rsun,psh_points_pol[:,0],psh_points_pol[:,1]) #interpolation using fig5_b theoretical pol
grid_Th = np.empty(np.shape(datatest)) * np.nan #same size as data but empty, it will Fill with thomson values
for idx,item in enumerate(enc_coor):
    grid_Th[item[1],item[0]] = th_step1[idx]
grid_Pol = np.empty(np.shape(datatest)) * np.nan #same size as data but empty, it will Fill with theoretical pol values
for idx,item in enumerate(enc_coor):
    grid_Pol[item[1],item[0]] = pol_step1[idx]
grid_Th = grid_Th[2300:2363,57:120]
grid_Pol = grid_Pol[2300:2363,57:120]
refw = int(bsz/2)
# I_0 = np.nanmean(datatest[yc-refw:yc+refw,xc-refw:xc+refw])
# # I_0 = np.nanmean(datatest[600-refw:600+refw,700-refw:700+refw])
# #I_0 = 160
# grid_Isc = np.empty(np.shape(datatest)) * np.nan #same sizes as data but empty, it will fill with I_sc values
# for idx,item in enumerate(enc_coor):
#     grid_Isc[item[1],item[0]]=datatest[item[1],item[0]]
# grid_Ne = ((grid_Isc/200)/I_0)*(1/grid_Th)
# grid_ne = grid_Ne/(0.725*100000000)
del datatest

"""Coupling with previous code"""
"""---------------------------"""
#Datos de GOES
tr = TimeRange(['2017-09-10 12:30', '2017-09-10 21:30'])
results = Fido.search(a.Time(tr), a.Instrument('XRS'))
files = Fido.fetch(results)
goes = TimeSeries(files)
tci = datetime.datetime.strptime('2017-09-10 12:28:41.40','%Y-%m-%d %H:%M:%S.%f')
tcf = datetime.datetime.strptime('2017-09-10 21:22:41.30','%Y-%m-%d %H:%M:%S.%f')
xrsa = goes.data['xrsb']
client = hek.HEKClient()
flares_hek = client.search(hek.attrs.Time(tr.start, tr.end),
                           hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')
xrgoes = np.array([xrsa[i] for i in range(len(xrsa)) if tci<=xrsa.index[i]<=tcf])
xrtiempo = np.array([xrsa.index[i] for i in range(len(xrsa)) if tci<=xrsa.index[i]<=tcf])
### Reading and sorting HMI meanvalue files ###
### Choose between averages or rebining files ###
mypathI = 'Stokes_IQUV_averages/'
# mypathI = 'Stokes_IQUV_rebin8x8/'
filenamesI = [mypathI+f for f in listdir(mypathI) if isfile(join(mypathI, f))]
I_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3I':
		I_files.append(item)
Q_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3Q':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		Q_files.append(item)
U_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3U':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		U_files.append(item)
V_files = []
for item in filenamesI:
	if item.split('_')[5] == 'TAI.3V':
#		if item.split('_')[-1].split('.')[2] == 'rebin8x8':
		V_files.append(item)
I_files.sort()
Q_files.sort()
U_files.sort()
V_files.sort()


###Reading files and computing stuff
rot_pol_ang = np.radians(90. - 8.81220944)

#Reading I files
tiempo = []                     #List of times
I_central_avg = []              #List of average I on solar centre
I_region_scc = []               #List of regions(arrays) per time
for item in I_files:
    print('Making interval %s for I'%item.split('.')[2])
    I_data = sunpy.map.Map(item).data
    I_head = sunpy.map.Map(item).meta
    I_ref = np.array(I_data[yc-refw:yc+refw,xc-refw:xc+refw])
    I_central_avg.append(np.nanmean(I_ref))
    del I_ref
    grid_Isc = np.empty(np.shape(I_data)) * np.nan #same sizes as data but empty, it will fill with I_sc values
    for cc in enc_coor:
        grid_Isc[cc[1],cc[0]]=I_data[cc[1],cc[0]]
    I_region_scc.append(grid_Isc[2300:2363,57:120])
    del grid_Isc
    tiempo.append(I_head['date-obs'])
    del I_data,I_head
I_region_scc_avg = [np.nanmean(item) for item in I_region_scc]           #List of average I in regions
tiempo = [datetime.datetime.strptime(i,'%Y-%m-%dT%H:%M:%S') for i in tiempo]    #times in datetime DateFormatter

#reading Files Q
Q_central_avg = []              #List of average I on solar centre
Q_region_scc = []               #List of regions(arrays) per time
for item in Q_files:
    print('Making interval %s for Q'%item.split('.')[2])
    Q_data = sunpy.map.Map(item).data
    Q_head = sunpy.map.Map(item).meta
    Q_ref = np.array(Q_data[yc-refw:yc+refw,xc-refw:xc+refw])
    Q_central_avg.append(np.nanmean(Q_ref))
    del Q_ref
    grid_Qsc = np.empty(np.shape(Q_data)) * np.nan #same sizes as data but empty, it will fill with I_sc values
    for cc in enc_coor:
        grid_Qsc[cc[1],cc[0]]=Q_data[cc[1],cc[0]]*np.cos(2*rot_pol_ang)         #multiplying by cos(2theta) rotates Q into Q'
    Q_region_scc.append(grid_Qsc[2300:2363,57:120])
    del grid_Qsc
Q_region_scc_avg = [np.nanmean(item) for item in Q_region_scc]                  #List of average Q in region per time
# Q_region_scc_avg = np.array(Q_region_scc_avg)*np.cos(2*rot_pol_ang)             #Q' considering U=0

###Computing I and Q backgrounds (Considering U = 0)

I_bkg = np.sqrt(np.mean(np.square(np.array(I_region_scc_avg[0:30]))))           #RMS first data
Q_bkg = np.sqrt(np.mean(np.square(np.array(Q_region_scc_avg[0:20]))))           #RMS first data

###Computing DI and DQ per time per pixel
# DI = []
# for inter in I_region_scc:
#     grid_di_tmp = np.empty(np.shape(inter)) * np.nan
#     for i in range(len(inter)):
#         for j in range(len(inter[i])):
#             grid_di_tmp[i,j] = math.hypot(0,inter[i,j]-I_bkg)
#     DI.append(grid_di_tmp)
#     del grid_di_tmp

DI = []
for inter in I_region_scc:
    grid_di_tmp = np.empty(np.shape(inter)) * np.nan
    for i in range(len(inter)):
        for j in range(len(inter[i])):
            if inter[i,j] > I_bkg:
                grid_di_tmp[i,j] = inter[i,j]-I_bkg
            elif inter[i,j] < I_bkg:
                grid_di_tmp[i,j] = 0.0
            else:
                grid_di_tmp[i,j] = inter[i,j]
    DI.append(grid_di_tmp)
    del grid_di_tmp

# DQ = []
# for inter in Q_region_scc:
#     grid_dq_tmp = np.empty(np.shape(inter)) * np.nan
#     for i in range(len(inter)):
#         for j in range(len(inter[i])):
#             grid_dq_tmp[i,j] = math.hypot(0,inter[i,j]-Q_bkg)
#     DQ.append(grid_dq_tmp)
#     del grid_dq_tmp


DQ = []
for inter in Q_region_scc:
    grid_dq_tmp = np.empty(np.shape(inter)) * np.nan
    for i in range(len(inter)):
        for j in range(len(inter[i])):
            if inter[i,j] > I_bkg:
                grid_dq_tmp[i,j] = inter[i,j]-Q_bkg
            elif inter[i,j] < Q_bkg:
                grid_dq_tmp[i,j] = 0.0
            else:
                grid_dq_tmp[i,j] = inter[i,j]
    DQ.append(grid_dq_tmp)
    del grid_dq_tmp

DQ_DI = np.array(DQ)/np.array(DI)                                               #Measured DQ_DI

Ne = []
for idx,inter in enumerate(DQ_DI):
    for i in range(len(inter)):
        for j in range(len(inter[i])):
            pol_m = inter[i,j]*100
            pol_t = grid_Pol[i,j]
            if pol_m<pol_t:
                grid_Ne = ((I_region_scc[idx]/I_bkg)/I_central_avg[idx])*(pol_t/pol_m)*(1/grid_Th)
            else:
                grid_Ne = ((I_region_scc[idx]/I_bkg)/I_central_avg[idx])*(1/grid_Th)
    Ne.append(grid_Ne)
Ne = np.array(Ne)
Ne_min = np.nanmin(Ne)
Ne_max = np.nanmax(Ne)


# for idx,item in enumerate(Ne):
#     plt.imshow(item,origin='lower',cmap='rainbow',vmin=Ne_min,vmax=Ne_max,extent=[0,63*arpi,0,63*arpi])
#     plt.title(r'Ne [cm$^{-2}$] at %s'%datetime.datetime.strftime(tiempo[idx],'%H:%M:%S'))
#     plt.colorbar()
#     plt.xlabel(r'arcsec')
#     plt.ylabel(r'arcsec')
#     plt.savefig('Ne_movie/%i.jpg'%idx,dpi=300)
#     plt.close()

# for idx,item in enumerate(DQ_DI):
#     plt.imshow(item,origin='lower',cmap='rainbow',vmin=0,vmax=5,extent=[0,63*arpi,0,63*arpi])
#     plt.title(r'<$\Delta Q/ \Delta I$> at %s'%datetime.datetime.strftime(tiempo[idx],'%H:%M:%S'))
#     plt.colorbar()
#     plt.xlabel(r'arcsec')
#     plt.ylabel(r'arcsec')
#     plt.savefig('DQ_DI_movie/%i.jpg'%idx,dpi=300)
#     plt.close()

# stdi = np.std(np.array(I_region_scc_avg[0:6]))
# bli = np.ones(len(I_r))*rmsi

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
