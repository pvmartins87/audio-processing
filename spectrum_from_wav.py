# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:34:58 2023

@author: Paulo Martins
Secop Austria GmbH
paulo.martins@secop.com
+43 664 8559 833

PostProcessor for CONTACT MIC AUDIO --------------------------------------

Usage:
    PostProcess audio files into spectra and total values in 1/3 oct. bands

!! IMPORTANT !!

The code will postprocess the values automatically.
Use Spyder IDE to open the results variable "df" and copy it to Excel

"""

print('Importing dependencies, please wait...')
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.io import wavfile
import easygui

print('Defining standard variables and loading file dialog...')
#%% USER INTERFACE - LOAD FILES
lfreq = .1
ufreq = 11220
plot_all_curves = False

if 'paths' not in locals():
    paths = easygui.fileopenbox(msg='Please select the files to calculate octave bands', title='Audio Files', default='*.wav', multiple=True)

if paths is None:
    sys.exit()

print('Processing selected files...')
#%% MAIN PROG

def Aweight(freq):
    Ra=(freq**4*12194**2)/((freq**2+20.6**2)*np.sqrt((freq**2+107.7**2)*(freq**2+737.9**2))*(freq**2+12194**2))
    A = 20*np.log10(Ra)+2
    return A

octs = [100,125,160,200,250,315,400,500,630,800,'1k','1.25k','1.6k','2k','2.5k','3.15k','4k','5k','6.3k','8k','10k']
octl = np.array([89.1,112,141,178,224,282,355,447,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220], dtype=float)
octc = np.array([100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000], dtype=float)

colnames = ['filename']
colnames.extend(octs)
colnames.append('LwA')

if plot_all_curves:
    fg, ax = plt.subplots(figsize=(16,9))

all_spectrograms = []
df = pd.DataFrame(columns=colnames)
for path in paths:
    filename = os.path.basename(path)

    
    sample_rate, samples = wavfile.read(path)
    try:
        samples = samples[:,0]
    except IndexError:
        pass
    
    samples = samples[samples!=0]
    
    freq, time, vals = signal.spectrogram(samples, sample_rate, nperseg=16384) # 65536 16384
    ix_ = (np.array(freq)>=lfreq) & (np.array(freq)<=ufreq)
    freq = freq[ix_]
    vals = vals[ix_,:]
    octaves = np.zeros((len(octc), len(time)))
    
    for f_ in range(len(octl)-1):
        fx_ = (freq>=octl[f_]) & (freq<=octl[f_+1])
        octaves[f_,:] = vals[fx_,:].sum(axis=0)
    
    valsdBA = 20*np.log10(vals) + np.repeat(np.reshape(Aweight(freq),(-1,1)),vals.shape[1],axis=1)
    octsdBA = 20*np.log10(octaves) + np.repeat(np.reshape(Aweight(octc),(-1,1)),vals.shape[1],axis=1)
    
    tots = 20*np.log10((10**(octsdBA/20)).sum(axis=0))
    
    avgOcts = np.ma.masked_invalid(octsdBA).mean(axis=1).data
    avgTots = np.ma.masked_invalid(tots).mean()
    data_to_append = [ filename ]
    data_to_append.extend(avgOcts)
    data_to_append.append(avgTots)
        
    df.loc[len(df)] = data_to_append
    print(f'{filename} - {avgTots:.1f} dB(A)')
    
    all_spectrograms.append({
        'filename':filename,
        'freq':freq,
        'time':time,
        'valsdBA':valsdBA,
        'octsdBA':octsdBA,
        'totals':tots,
        })
    
    if plot_all_curves:
        ax.plot(data_to_append[3:], label=filename)
    
if plot_all_curves:
    ax.set(
           xlabel='freq. [Hz]',
           xticks=range(len(colnames[3:])),
           xticklabels=colnames[3:],
           ylabel='SPL [dB]',
           title='Contact Mic (corrected)')
    ax.grid()
    ax.legend()

print('Writting excel report...')
df.to_excel(os.path.join(os.path.dirname(path), 'Generated_Report.xlsx'))

print('Done! Closing...')
