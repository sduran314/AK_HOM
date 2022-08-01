#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:40:41 2022

@author: sduran
"""

from obspy.clients.fdsn import Client
client = Client('IRIS')

from obspy import UTCDateTime
from scipy.signal import spectrogram
import numpy as np
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
sta = 'O17K'
loc = 'EP'
chan = 'BDF'

SPEC_WIN = 50 # spectrogram window length

SAVE = False #save figure or not


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

# compute spectrogram using scipy

# get unfiltered trace for spectrogram
tr = st_rem[0]

# compute number of samples per spectrogram window
nper = int(SPEC_WIN*tr.stats.sampling_rate)

# compute spectrogram
f, t, Pspec = spectrogram(tr.data, fs=tr.stats.sampling_rate, window='hann',
                          scaling='density', nperseg=nper, noverlap=nper*.5)

# convert spectrogram to dB, using 20e-6 as a reference pressure
PspecdB = 10 * np.log10( abs(Pspec) / np.power(20e-6, 2))

# define min and max for spectrogram colormap
cmin = np.nanpercentile(PspecdB, 15)
cmax = np.nanpercentile(PspecdB, 99.5)

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
#ax2.set_ylim(-20, 20)
ax2.xaxis_date()

#Spectogram
ax3 = plt.subplot(313)

# plot the spectrogram image
im = ax3.imshow(PspecdB, extent=[tvec[0], tvec[-1], f[0], f[-1]],
                origin='lower', aspect='auto', interpolation=None,
                cmap='magma_r')
im.set_clim(cmin, cmax)
ax3.set_ylabel('Frequency [Hz]')
ax3.set_ylim(0.1, 15)
ax3.set_xlim(tvec[0], tvec[-1])
ax3.xaxis_date()

# make a colorbar for the spectrogram
pos1 = ax3.get_position()
cloc = [pos1.x0+pos1.width+.03, pos1.y0, 0.02, pos1.height]
cbaxes = fig.add_axes(cloc)
hc = plt.colorbar(im, cax=cbaxes)
hc.set_label('PSD [dB re 20\u03bc$Pa^2$/Hz]')

# save the figure or not
if SAVE:
    fig.savefig(tr.id +'_spec.png', bbox_inches='tight', dpi=250)


