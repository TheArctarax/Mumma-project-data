from __future__ import division, print_function
import gwmemory
import matplotlib
from gwmemory import utils as utils
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

# Q1 : Start by looking at the frequency_domain_memory definition you quoted. Here, we want
#      to define a new function called `frequency_domain_transform`, which takes in two
#      arguments - time_domain_waveform (which is basically either the memory, or the oscillatory)
#      and times (the sampling times array), and do a fourier transform to the frequency
#      domain, and finally, returns the frequency_domain_waveform and sampling frequencies
#      to the user (so two products). Modify the code in frequency_domain_memory to attain
#      the above goals.

def frequency_domain_transform(time_domain_waveform, times):
    if :  
      time_domain_strain, times = time_domain_oscillatory(
        model=model, q=q, total_mass=total_mass, spin_1=spin_1, spin_2=spin_2,
        distance=distance, inc=inc, phase=phase, **kwargs)
    elif :
      time_domain_strain, times = time_domain_memory(
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


# GW waveform-only definition
# A surrogate object has the following attributes: frequency_domain_memory (a 
# pycbc.frequencyseries object that has both the ['plus'] and ['cross'] sub-arrays), 
# and frequency_domain_memory which also has ['plus'] and ['cross']). Calling
# these attributes returns both the pycbc frequency series and the sampling
# frequencies.

# Now, basically you can still do the same (using time_domain_oscillatory) here. 
oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)

# Q2 : Now it's time to apply our newly defined function `frequency_domain_transform`. Use
#      this to transform the osciallory time domain waveform to the frequency domain. We
#      use oscillatory_tilde, frequencies to store the products. Complete the following 
#      line of code.
oscillatory_tilde, frequencies = surr.frequency_domain_transform(oscillatory, times)



# GW memory definition
# Same thing happens for the memory
memory, times = surr.time_domain_memory(inc=inc, phase=pol)

# Q3 : Again, apply our new function to transform the memory into the frequency domain.
memory_tilde, frequencies = surr.frequency_domain_transform(memory, times)


# Plot of the frequency domain waveform.

# Q4 : It's usually best to plot these figures in a log-log scale. Figure how to do so
#      and modify the following line. Note that you also have to change the arguments:
#      In particular, you should plot frequencies as the x-axis, and the **magnitude**
#      of the amplitude_tilde times square root of frequencies as the y-axis.
fig = figure(figsize=(6, 6)
loglog(frequencies, oscillatory_tilde['plus'][:]*np.sqrt(frequencies[:]), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_waveform.pdf')

tight_layout()
show()
close()


# Plot of the frequency domain memory.

# Q5 : Do the same for the memory part.
fig = figure(figsize=(6, 6)
loglog(frequencies, memory_tilde['plus'][:]*np.sqrt(frequencies[:]), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_mem.pdf')

tight_layout()
show()
close()


# Plot of the frequency domain waveform + memory.

# Q6 : And of course, again for the waveform + memory.
fig = figure(figsize=(6, 6))
plot(frequencies, (oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(frequencies), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD $[Hz$^{-1/2}$]$')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('waveform_asd_wave_and_mem.pdf')

tight_layout()
show()
close()


# Plot of the frequency domain waveform, memory, and noise curves.

# Q7 : Finally, do the same for the full plot.
fig = figure(figsize=(6, 6))
loglog(frequencies, (oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(frequencies), color='r')
loglog(frequencies, oscillatory_tilde['plus'], color='g')
loglog(frequencies, memory_tilde['plus'], color='b')

# Here, I would like to call noise_curve() for superposition with the memory and waveform curves.
# Q8 : Finally, to plot the noise curve, simply load the txt files and plot them again here.
dfl=pd.read_csv('L1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfh=pd.read_csv('H1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfv=pd.read_csv('V1_O2_Sensitivity_strain_asd.txt',sep="\s+",index_col=False,header=None, names=["frequency","asd"])

loglog(dfl['frequency'], dfl['asd'], color='tab:blue', label='L1')
loglog(dfh['frequency'], dfh['asd'], color='tab:red', label='H1')
loglog(dfv['frequency'], dfv['asd'], color='tab:purple', label='V1')

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
