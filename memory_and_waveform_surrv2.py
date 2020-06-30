from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
from pycbc.waveform import get_td_waveform
import numpy as np
rcParams["font.family"] = "Times New Roman"
np.seterr(divide='ignore', invalid='ignore')


# Variable assignments
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0
d=600
M=60
q=4


# Sample space definition for the memory's t-axis. Purposely set to begin, end, and have the same number of points as the
# original waveform so that superposition of the timeseries is possible.
start_time=-0.08
end_time=0.02
times = np.linspace(start_time, end_time, 10001)


# GW waveform with memory definition
# The sub-function waveforms.surrogate.Surrogate generates a surrogate object.
surr = gwmemory.waveforms.surrogate.Surrogate(q=q, name='nrsur7dq4', spin_1=S1, spin_2=S2, total_mass=M, distance=400, times=times)


# GW waveform only definition
# A surrogate object has the following attributes: time_domain_memory (a 
# pycbc.timeseries object that has both the ['plus'] and ['cross'] sub-arrays), 
# and time_domain_memory which also has ['plus'] and ['cross']). Calling
# these attributes returns both the pycbc timesseries and the sampling time
# arrays (which you store it as times here).
oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)

# GW memory definition
memory, times = surr.time_domain_memory(inc=inc, phase=pol)

# Plot of GW memory
# By using the above attributes, you can easily plot out the memory,
# original waveform and waveform + memory. 
fig = figure(figsize=(12, 6))
plot(times, memory['plus']*(10**22), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.02)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('memory.pdf')

tight_layout()
show()
close()


# Plot of memory-less waveform
fig = figure(figsize=(12, 6))
plot(times, oscillatory['plus']*(10.0**22.0), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.02)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('original_waveform.pdf')

tight_layout()
show()
close()


# plot of oscillatory + memory components
fig = figure(figsize=(12, 6))

fig.add_subplot(2, 1, 2)
plot(times, (oscillatory['plus']+memory['plus'])*(10.0**22.0), color='r', label=r'Waveform $\plus$ Memory')
plot(times, oscillatory['plus']*(10.0**22.0), linestyle='--' , color='tab:purple', label='Original Waveform')
#plot(times, memory['plus']*(10.0**22.0), linestyle='dotted', color='b', label='Memory')
axhline(0, linestyle=':', color='k')
xlim(-0.04, 0.015)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)
#tight_layout()
#show()
#close()


fig.add_subplot(2, 1, 1)
plot(times, (oscillatory['plus'][:] + memory['plus'][:])*(10.0**22.0), color='r', label=r'Waveform $\plus$ Memory')
plot(times, memory['plus']*(10.0**22.0), linestyle='dotted', color='b', label='Memory')
axhline(0, linestyle=':', color='k')
xlim(-0.05, times[-1])
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('combined.pdf')

tight_layout()
show()
close()



