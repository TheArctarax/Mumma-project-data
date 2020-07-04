from __future__ import division, print_function
import gwpy
import gwmemory
import matplotlib
from gwmemory import utils as utils
from matplotlib.pyplot import *
import numpy as np
import pandas as pd
from matplotlib.pyplot import rcParams
rcParams["font.family"]="Times New Roman"
rcParams['axes.unicode_minus'] = False
np.seterr(divide='ignore', invalid='ignore')


'''
This program calculates and plots asd's for waveforms and their
associated memory. Comparison is then made with noise asd's by superposing plots.
Noise asd's come from noise_curve.py
'''


# GWmemory natively features a frequency domain converter for memory but not for waveforms.
# This method fixes that problem by following the same structure as frequency_domain_memory().
def frequency_domain_transform(time_domain_strain, times):
  sampling_frequency = 1 / (times[1] - times[0])
  frequencies = None
  frequency_domain_strain = dict()

  for key in time_domain_strain:
    frequency_domain_strain[key], frequencies =\
      utils.nfft(time_domain_strain[key], sampling_frequency)

  return frequency_domain_strain, frequencies


# Variable assignments
f_lower=15
f_upper=312
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0
d=420
M=66.2
q=1.16


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


# Now it's time to apply our newly defined function `frequency_domain_transform`. Use
# this to transform the osciallory time domain waveform to the frequency domain. We
# use oscillatory_tilde, frequencies to store the products.
oscillatory_tilde, frequencies = frequency_domain_transform(oscillatory, times)


# GW memory definition
# Same thing happens for the memory
memory, times = surr.time_domain_memory(inc=inc, phase=pol)


# Again, apply our new function to transform the memory into the frequency domain.
memory_tilde, frequencies = frequency_domain_transform(memory, times)


# Plot of the frequency domain waveform.

# It's usually best to plot these figures in a log-log scale. Figure how to do so
# and modify the following line. Note that you also have to change the arguments:
# In particular, you should plot frequencies as the x-axis, and the **magnitude**
# of the amplitude_tilde times square root of frequencies as the y-axis.
fig = figure(figsize=(6, 6))
loglog(frequencies, abs(oscillatory_tilde['plus'][:]*np.sqrt(frequencies[:])), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD [Hz$^{-1/2}$]')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)
grid(False)

savefig('waveform_asd_waveform.pdf')

tight_layout()
#show()
close()


# Plot of the frequency domain memory.

# Q5 : Do the same for the memory part.
fig = figure(figsize=(6, 6))
loglog(frequencies, abs(memory_tilde['plus'][:]*np.sqrt(frequencies[:])), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD [Hz$^{-1/2}$]')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)
grid(False)

savefig('waveform_asd_mem.pdf')

tight_layout()
#show()
close()


# Plot of the frequency domain waveform + memory.

# Q6 : And of course, again for the waveform + memory.
fig = figure(figsize=(6, 6))
loglog(frequencies, abs((oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(frequencies)), color='r')

xlim(0, 5000)
xlabel('Frequency [Hz]')
ylabel(r'ASD [Hz$^{-1/2}$]')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)
grid(False)

savefig('waveform_asd_wave_and_mem.pdf')

tight_layout()
#show()
close()


# Plot of the frequency domain waveform, memory, and noise curves.

# begin by generating truncated frequency array to remove extraneous Fourier modulation outside the region of interest.
f_new = np.zeros(len(frequencies))

for i in range(len(frequencies)):
  if frequencies[i] > f_upper:
    f_new[i]=f_upper
  elif frequencies[i] < f_lower:
    f_new[i]=f_lower
  else:
    f_new[i]=frequencies[i]

# Finally, do the same for the full plot.
fig = figure(figsize=(6, 6))
loglog(f_new, abs((oscillatory_tilde['plus'][:] + memory_tilde['plus'][:])*np.sqrt(frequencies)), color='tab:orange', linewidth=10, zorder=1, label='Waveform + Memory')
loglog(f_new, abs(oscillatory_tilde['plus']*np.sqrt(frequencies)), color='g', linewidth=12, zorder=1, label='Waveform')
loglog(frequencies, abs(memory_tilde['plus']*np.sqrt(frequencies)), color='b', linewidth=12, zorder=1, label='Memory')

# Finally, to plot the noise curve, simply load the txt files and plot them again here.
dfl=pd.read_csv('L1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfh=pd.read_csv('H1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfv=pd.read_csv('V1_O2_Sensitivity_strain_asd.txt',sep="\s+",index_col=False,header=None, names=["frequency","asd"])

loglog(dfl['frequency'], dfl['asd'], color='gwpy:ligo-livingston', zorder=0 , label='L1 noise')
loglog(dfh['frequency'], dfh['asd'], color='gwpy:ligo-hanford', zorder=0 , label='H1 noise')
loglog(dfv['frequency'], dfv['asd'], color='gwpy:virgo', zorder=0, label='V1 noise')

xlim(10, 1000)
ylim(5*(10**(-24)), 10**(-21))
xlabel('Frequency [Hz]')
ylabel(r'ASD [Hz$^{-1/2}$]')
xscale('log')
yscale('log')
legend(loc='upper right', prop={'size':18})
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)
grid(False)

savefig('waveform_asd_wave_and_mem_and_noise.pdf')

tight_layout()
show()
close()

