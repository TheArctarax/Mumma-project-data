from __future__ import division, print_function
import gwmemory
import matplotlib
from gwmemory import utils as utils
matplotlib.use('Agg')
from matplotlib.pyplot import *
import numpy as np
import pandas as pd
import wget
import itertools
from matplotlib.pyplot import rcParams
rcParams["font.family"]="Times New Roman"
rcParams['axes.unicode_minus'] = False
np.seterr(divide='ignore', invalid='ignore')

'''
when finished, this program will calculate and plot asd's for waveforms and their
associated memory. Comparison is then made with noise asd's by superposing plots.
Noise asd's come from noise_curve.py
'''


# GWmemory natively features a frequency domain converter for memory but not for waveforms.
# This method fixes that problem by following the same structure as frequency_domain_memory().

# Here, I make a little modification and defined a new function frequency domain transform. It
# takes in two argument - the timeseries time_domain_waveform and the sampling time series times,
# and return the frequency series with sampling frequencies.

def frequency_domain_transform(time_domain_waveform, times):
    sampling_frequency = 1 / (times[1] - times[0])
    frequencies = None
    frequency_domain_strain = dict()
    for key in time_domain_waveform:
        frequency_domain_strain[key], frequencies =\
            utils.nfft(time_domain_waveform[key], sampling_frequency)

    return frequency_domain_strain, frequencies


# Variable assignments
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0
d=600
M=60
q=1.99


# Sample space definition for the memory's t-axis.
start_time=-0.51
end_time=0.02
times = np.linspace(start_time, end_time, 10001)


# GW waveform with memory definition
# The sub-function waveforms.surrogate.Surrogate generates a surrogate object.
surr = gwmemory.waveforms.surrogate.Surrogate(q=q, name='nrsur7dq2', spin_1=S1, spin_2=S2, total_mass=M, distance=d, times=times)


# GW waveform only definition
# A surrogate object has the following attributes: frequency_domain_memory (a 
# pycbc.frequencyseries object that has both the ['plus'] and ['cross'] sub-arrays), 
# and frequency_domain_memory which also has ['plus'] and ['cross']). Calling
# these attributes returns both the pycbc frequency series and the sampling
# frequencies.
oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)

# Here, I used frequency_domain_transform to transform to the frequency domain:
oscillatory_tilde, times_tilde = frequency_domain_transform(oscillatory, times)
#oscillatory_tilde['cross'], times_tilde = frequency_domain_transform(oscillatory['cross'], times)

# GW memory definition
memory, times = surr.time_domain_memory(inc=inc, phase=pol)

# Same for the memory:
memory_tilde, times_tilde = frequency_domain_transform(memory, times)
#memory_tilde['cross'], times_tilde = frequency_domain_transform(memory['cross'], times)

# Plot of the frequency domain waveform.
fig = figure(figsize=(12, 6))
loglog(times_tilde, abs(oscillatory_tilde['plus']), color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_waveform.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain memory.
fig = figure(figsize=(12, 6))
loglog(times_tilde, abs(memory_tilde['plus'][:] * np.sqrt(times_tilde[:])), color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_mem.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain waveform + memory.
fig = figure(figsize=(12, 6))
loglog(times_tilde, abs((oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(times_tilde[:])), color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_wave_and_mem.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain waveform, memory, and noise curves.
fig = figure(figsize=(12, 6))
loglog(times_tilde, abs((oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(times_tilde[:])), color='r')
loglog(times_tilde, abs(oscillatory_tilde['plus'][:]*np.sqrt(times_tilde[:])), color='g')
loglog(times_tilde, abs(memory_tilde['plus'][:]*np.sqrt(times_tilde[:])), color='b')
# Here, I would like to call noise_curve() for superposition with the memory and waveform curves.
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_wave_and_mem_and_noise.pdf')

tight_layout()
show()
close()

