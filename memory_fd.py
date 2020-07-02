from __future__ import division, print_function
import gwmemory
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
def frequency_domain_oscillatory(
        model=None, q=None, total_mass=None, spin_1=None, spin_2=None,
        distance=None, inc=None, phase=None, **kwargs):
    """
    Calculate the frequency domain oscillatory waveform according to __reference__.
    Parameters
    ----------
    model: str
        Name of the model, this is used to identify waveform approximant,
        e.g., NRSur7dq2, IMRPhenomD, MWM, etc.
    q: float
        Mass ratio of the binary being considered.
    total_mass: float
        Total mass of the binary being considered in solar units.
    spin_1: array
        Dimensionless spin vector of the more massive black hole.
    spin_2: array
        Dimensionless spin vector of the less massive black hole.
    distance: float
        Distance to the binary in MPC.
    inc: float
        Inclination of the binary to the line of sight.
        If not provided, spherical harmonic modes will be returned.
    phase: float
        Binary phase at coalescence.
        If not provided, spherical harmonic modes will be returned.
    kwargs: dict
        Additional model-specific keyword arguments.
    Returns
    -------
    frequency_domain_strain: dict
        Memory frequency series, either in spherical harmonic modes or
        plus/cross polarisations.
    frequencies: array-like
        Frequency series corresponding to the memory waveform.
    """
    time_domain_strain, times = time_domain_oscillatory(
        model=model, q=q, total_mass=total_mass, spin_1=spin_1, spin_2=spin_2,
        distance=distance, inc=inc, phase=phase, **kwargs)
    sampling_frequency = 1 / (times[1] - times[0])

    frequencies = None
    frequency_domain_strain = dict()
    for key in time_domain_strain:
        frequency_domain_strain[key], frequencies =\
            utils.nfft(time_domain_strain[key], sampling_frequency)

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
oscillatory, frequencies = surr.frequency_domain_oscillatory(inc=inc, phase=pol)

# GW memory definition
memory, frequencies = surr.frequency_domain_memory(inc=inc, phase=pol)


# Plot of the frequency domain waveform.
fig = figure(figsize=(12, 6))
plot(times, oscillatory['plus'], color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain memory.
fig = figure(figsize=(12, 6))
plot(times, memory['plus'], color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain waveform + memory.
fig = figure(figsize=(12, 6))
plot(times, oscillatory['plus'][:] + memory['plus'][:], color='r')
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd.pdf')

tight_layout()
show()
close()

# Plot of the frequency domain waveform, memory, and noise curves.
fig = figure(figsize=(12, 6))
plot(times, oscillatory['plus'][:] + memory['plus'][:], color='r')
plot(times, oscillatory['plus'], color='g')
plot(times, memory['plus'], color='b')
# Here, I would like to call noise_curve() for superposition with the memory and waveform curves.
xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd.pdf')

tight_layout()
show()
close()

