
"""
Square electronic desnsity, contourf plots
Ne_v6.npy ----------------------------------------------------------------------------> Densities [cm^2] calculated
tiempos_v6.npy -----------------------------------------------------------------------> Times from HMI data in datetime format


Author: Juan Camilo Guevara Gomez
Latest version: 21 January 2018
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
import imageio




import pandas as pd
import seaborn as sns
import matplotlib.animation as animation

warnings.filterwarnings('ignore')
np.seterr(divide='ignore', invalid='ignore')

#Read densities and time data

Ne = np.load('Ne_v6.npy')
Ne[Ne<=0]=np.nan
times = np.load('times_v6.npy')

Ne_log10 = np.array([np.log10(ne) for ne in Ne])

Ne_min = np.nanmin(Ne_log10)
Ne_max = np.nanmax(Ne_log10)
levels = np.linspace(Ne_min,Ne_max,50)
print(Ne_min,Ne_max)

tiempos = np.array([datetime.datetime.strftime(t,'%Y-%m-%d %H:%M:%S') for t in times])

def fig2data ( fig ):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring ( fig.canvas.tostring_rgb(), dtype=np.uint8 )
    buf.shape = ( w, h,3 )
 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    # buf = np.roll ( buf, 3, axis = 2 )
    return buf

# for idx,ne in enumerate(Ne_log10):
# 	CS = plt.contourf(ne, levels, alpha=0.7, cmap=plt.cm.rainbow, origin='lower')
# 	cbar = plt.colorbar(CS)
# 	cbar.ax.set_ylabel(r'Logaritmic square density')
# 	plt.title('%s'%tiempos[idx])
# 	plt.show()

def plot_for_offset(levels,ne,tiempo,i):
    # Data for plotting
    print(i)
    fig, ax = plt.subplots()
    CS = ax.contourf(ne, levels, alpha=0.7, cmap=plt.cm.rainbow, origin='lower')
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel(r'Logaritmic square density')
    ax.set_title('%s'%tiempo)
    ax.grid()
    # # Used to return the plot as an image rray
    # fig.canvas.draw()       # draw the canvas, cache the renderer
    # image = np.frombuffer(fig.canvas.tostring_rgb(),dtype='uint8')
    # # print(np.shape(image),type(image))
    # print(fig.canvas.get_width_height()[::-1]+ (3,))
    # image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    image = fig2data(fig)

    return image

kwargs_write = {'fps':5.0, 'quantizer':'nq'}
imageio.mimsave('./densities.gif', [plot_for_offset(levels, Ne_log10[i],tiempos[i],i) for i in range(len(Ne_log10[:5]))], fps=5)

