#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:40:41 2022

@author: sduran
"""

from obspy.clients.fdsn import Client
client = Client('IRIS')

from obspy import UTCDateTime
from obspy import read

from matplotlib import pyplot as plt
plt.ion()

#USGS Origin:Time of Hunga,Tonga eruption event. 
time = UTCDateTime('2022-01-15T04:14:45.000')
print(time)

#Rounding start time to nearest whole number.
starttime =time + 8*3600
print(starttime)

#Final major eruption detected at 8:31 UTC 
#according to Matoza & Fee paper on Science.com
#Rounding end time to the nearest whole number. 
endtime = starttime + 4*3600
print(endtime)

net = 'AK'
sta = 'HOM'
loc = 'EP'
chan = 'BDF'

#To view data stream
st = client.get_waveforms(net, sta, loc, chan, starttime, endtime, attach_response=True)
st_rem = st.copy()
#print(st_rem)
#print(st)

for tr in st_rem:
    fs_resp = tr.stats.sampling_rate
    pre_filt = [0.0001, 0.0002, fs_resp/2-2, fs_resp/2] #pre-filt for response removal
    tr.remove_response(pre_filt=pre_filt, output='VEL', water_level=None)
#st_rem.remove_response(output = 'VEL')
#st.plot()
#st_rem.plot()

#View Raw, Raw After Time, and Response Removed
#st_rem =st.copy()
#st_rem.remove_response(output = 'VEL', plot=True)

#Data Filtering
st_filt = st_rem.copy()
st_filt.filter('bandpass', freqmin=0.1, freqmax=10)
#st.plot()
#st_filt.plot()

#Response removal; can be VEL, ACC, DISP
#st_filt.remove_response(output = 'VEL')
#st_filt.plot();


#David added code below
tvec = st_filt[0].times(type='matplotlib')
fig = plt.figure()
ax1 = plt.subplot(311)
ax1.plot(tvec, st_rem[0].data, 'k-')
ax1.set_xlim(tvec[0], tvec[-1])
ax1.xaxis_date()
ax1.set_ylabel('Pressure [Pa]')
ax1.set_title(st_rem[0].id)

ax2 = plt.subplot(312)
ax2.plot(tvec, st_filt[0].data, 'k-')
ax2.set_xlim(tvec[0], tvec[-1])
ax2.xaxis_date()

#Spectogram
ax3 = plt.subplot(313)
st_rem[0].spectrogram(log=True)
plt.show()
