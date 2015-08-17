#!/usr/bin/env python

# make GMVs!
# Emily Wolin, USArray data course, Aug 2015
# emilyw@earth.northwestern.edu

import numpy as np
import obspy
from obspy.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib import animation

#
# Set up time window and filter bands
#
client = Client("IRIS")
origintime = UTCDateTime("2011-08-23T17:51:00.000")
starttime = UTCDateTime("2011-08-23T17:57:00.000")
endtime = starttime+300
#starttime = origintime + 200
#endtime = starttime+500
#starttime = UTCDateTime("2013-05-24T05:54:49")
#endtime = UTCDateTime("2013-05-24T06:29:49")
netlist = "TA"
stalist = "*"
loclist = ""
chanlist = "LHZ"
#pre_filt = (0.001, 0.006, 0.2, 0.4)
pre_filt = (0.001, 0.002, 0.2, 0.4)
#fmin = 0.006
fmin = 0.002
fmax = 0.05

# Get waveform data and station inventory
# remove response before writing to make processing easier if re-read later
#print('Requesting waveform data')
#st = client.get_waveforms(netlist, stalist, loclist, chanlist, starttime, endtime, attach_response=True) 
#print('Finished fetching data')
#st.remove_response(output="DISP", pre_filt=pre_filt)
#st.write('virginia.mseed', format='MSEED')
##st.plot()

# Or comment out the above and use the line below to read local data
st = obspy.read('virginia.mseed')
st.filter("bandpass", freqmin=fmin, freqmax=fmax)

#
# Normalize traces and find absolute max value (for plotting)
#
stalist_obtained = []
absmax = -999999
st.normalize()
for tr in st:
    stalist_obtained.append(tr.stats.station)
    tr_abs_max = np.max(np.abs(tr.data))
    if tr_abs_max > absmax:
        absmax = tr_abs_max

# Get station inventory so we have lats/lons
stalist_obtained_string = str.join(",",stalist_obtained)
inventory = client.get_stations(network=netlist, station=stalist_obtained_string, location=loclist, channel=chanlist, starttime=starttime, endtime=endtime)

# Normalize the color map (can adjust to play with saturation)
cmap = cm.bwr
#norm = mpl.colors.Normalize(vmin=-absmax, vmax=absmax)
norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
print(-absmax, absmax)

#
# Get list of station lat/lons
# 
lats = []
lons = []
codes = []
nets = []
for net in inventory:
    for sta in net:
        lats.append(sta.latitude)
        lons.append(sta.longitude)
        codes.append(sta.code)
        nets.append(net.code)
        
# Set up basemap and convert lat/lons to map coords
mapres = 'c'
#mapres = 'l'  # looks nicer but takes longer to draw
map = Basemap(projection = 'merc', llcrnrlat=20, llcrnrlon=-130, urcrnrlat=55, urcrnrlon=-70)
x, y = map(lons, lats)

#
# Set up figure
#
fig = plt.figure()
gs = gridspec.GridSpec(4, 1)
markersize = 15
refstaindex = 223 # seismogram that will be circled in yellow on map and plotted on ax_seis
c = np.zeros(len(st)) # blank array to hold values at each time step
t = range(len(st[refstaindex].data)) # time values for plotting reference seismogram

ax_map = plt.subplot(gs[0:3, :])
ax_seis = plt.subplot(gs[3, :])
map = Basemap(projection='merc', llcrnrlat=20, llcrnrlon=-130, urcrnrlat=55, urcrnrlon=-70, ax=ax_map, resolution=mapres, area_thresh=2000)

#
# Draw frames, one for each sample:
#
print('Drawing {0} frames...'.format(len(st[0].data)))
for i in range(len(st[0].data)):
    map.drawcoastlines()
# Get values at this time step:
    for j in range(len(st)):
        c[j] = st[j].data[i]

    if i%50 == 0:
        print('Frame {0}'.format(i))

# Draw dots on map
    map.scatter(x,y, c=c, cmap=cmap, norm=norm, s=markersize)
    map.scatter(x[refstaindex], y[refstaindex], linewidth=2, edgecolor='yellow', s=markersize)

# Draw seismogram and time marker for reference
    ax_seis.scatter(t, st[refstaindex].data, c=st[refstaindex].data, cmap=cmap, norm = norm, linewidth=0)
    ax_seis.plot((i,i),(-1,1))

# Save with zero-padded name for ease of use w/ffmpeg
    plt.savefig('{0:04d}.png'.format(i))
    ax_map.clear()
    ax_seis.clear()
    c[:] = np.nan


# To Do:
# - Add event location (if in map area) or great-circle path (if teleseism)
# - THREE COMPONENTS!  (w/matplotlib quiver?)
# - make sure seismogram x-axis really is time not just sample number (works for LH* data at 1 sample/sec but won't for anything else...)
# - add scale bar and/or colorbar on y-axis to show displacement
# - oooooh!!!  could add small TauP-style xsection of globe and animate body and surface waves approaching the array 
# - add arrivals w/TauP
# - automatically window for a given event and frequency band?
# - more sophisticated way to make sure we're getting the same time points for data at different sample rates (probably use decimate method?)
# - option to select time step (right now, just plots 1 frame/sample)
# - adjust scaling/clipping of seismograms (currently covers whole data range)
# - add great-circle path and state boundaries!
