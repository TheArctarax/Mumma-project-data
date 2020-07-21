from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
from pycbc.waveform import get_td_waveform
import numpy as np
rcParams["font.family"] = "Times New Roman"
rcParams['axes.unicode_minus'] = False
np.seterr(divide='ignore', invalid='ignore')


# Variable assignments
sampling_frequency = 4096
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0.0
d=600
M=60
q=1.0


# Sample space definition for the memory's t-axis. Purposely set to begin, end, and have the same number of points as the
# original waveform so that superposition of the timeseries is possible.
start_time=-0.5
end_time=0.02
times = np.linspace(start_time, end_time, sampling_frequency * (end_time-start_time))


# GW waveform with memory definition
# The sub-function waveforms.surrogate.Surrogate generates a surrogate object.
surr = gwmemory.waveforms.surrogate.Surrogate(q=q, name='nrsur7dq2', spin_1=S1, spin_2=S2, total_mass=M, distance=d, times=times)
'''
Got omega_ref = 0.0093 < 0.0161 = omega_0, too small
Got omega_ref = 0.0093 < 0.0161 = omega_0, too small
Got omega_ref = 0.0093 < 0.0161 = omega_0, too small
Got omega_ref = 0.0093 < 0.0161 = omega_0, too small
Got omega_ref = 0.0155 < 0.0161 = omega_0, too small
Got omega_ref = 0.0093 < 0.0185 = omega_0, too small
Got omega_ref = 0.0093 < 0.0167 = omega_0, too small!
m = 104.21 is lim between these errors and this error:
cannot convert float infinity to integer
'''
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
fig = figure(figsize=(6, 6))
plot(times, memory['plus']*(10**22), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.02)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)
rc('axes', labelsize=16)

# savefig('memorysurr.pdf')

tight_layout()
# show()
close()


# Plot of memory-less waveform
fig = figure(figsize=(12, 6))
plot(times, oscillatory['plus']*(10.0**22.0), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.02)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)
rc('axes', labelsize=16)

savefig('original_waveformsurr.pdf')

tight_layout()
show()
close()


# plot of oscillatory + memory components
fig = figure(figsize=(9, 4.5))

fig.add_subplot(2, 1, 2)
plot(times, (oscillatory['plus']+memory['plus'])*(10.0**22.0), color='r', label=r'Waveform $\plus$ Memory')
plot(times, oscillatory['plus']*(10.0**22.0), linestyle='--' , color='tab:purple', label='Original Waveform')
axhline(0, linestyle=':', color='k')
xlim(-0.04, 0.015)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left', prop={'size':13})
rc('xtick', labelsize=11)
rc('ytick', labelsize=11)
rc('axes', labelsize=13)


fig.add_subplot(2, 1, 1)
plot(times, (oscillatory['plus'][:] + memory['plus'][:])*(10.0**22.0), color='r', label=r'Waveform $\plus$ Memory')
plot(times, memory['plus']*(10.0**22.0), linestyle='dotted', color='b', label='Memory')
axhline(0, linestyle=':', color='k')
xlim(-0.5, times[-1])
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left', prop={'size':13})
rc('xtick', labelsize=11)
rc('ytick', labelsize=11)
rc('axes', labelsize=13)

#savefig('combinedsurr.pdf')

tight_layout()
#show()
close()


# GW waveform with memory definition
# The sub-function waveforms.approximant.Approximant generates an Approximant object.
# You can input whatever waveform approximant that is acceptable by PyCBC in the argument
# name. Note that I have made a small change to the waveform.approximant.py script.
start_time=-0.51
end_time=0.02
times=np.linspace(start_time, end_time, 10001)
approx = gwmemory.waveforms.approximant.Approximant(q=q, name="IMRPhenomD", spin_1=S1, spin_2=S2, total_mass=M, distance=d, times=times)

# GW waveform only definition
# A surrogate object has the following attributes: time_domain_memory (a
# pycbc.timeseries object that has both the ['plus'] and ['cross'] sub-arrays),
# and time_domain_memory which also has ['plus'] and ['cross']). Calling
# these attributes returns both the pycbc timesseries and the sampling time
# arrays (which you store it as times here).
oscillatoryphenom, timesphenom = approx.time_domain_oscillatory(inc=inc, phase=pol)

# GW memory definition
memoryphenom, timesphenom = approx.time_domain_memory(inc=inc, phase=pol)

# Plot of Superposed GW memory using IMRPhenomD and NRSur7dq2
fig = figure(figsize=(4.5, 4.5))
plot(times, memory['plus']*(10**22), color='r', label='NRSur7dq2')
plot(timesphenom, memoryphenom['plus']*(10**22), color='b', linestyle='dotted', label='IMRPhenomD')
axhline(0, linestyle=':', color='k')
xlim(-0.08, 0.02)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left', prop={'size': 12})
rc('xtick', labelsize=11)
rc('ytick', labelsize=11)
rc('axes', labelsize=13)

#savefig('imrsurrmemory.pdf')

tight_layout()
#show()
close()

# Plot of superposed IMRPhenomD and NRSurr7dq2 waveforms
fig = figure(figsize=(4.5, 4.5))
plot(times, oscillatory['plus']*(10.0**22.0), color='r', label='NRSur7dq2')
plot(timesphenom, oscillatoryphenom['plus']*(10.0**22.0), linestyle='dotted', color='b', label='IMRPhenomD')
axhline(0, linestyle=':', color='k')
xlim(-.08, 0.02)
xlabel('Time (s)')
legend(loc='lower left', prop={'size': 13})
ylabel(r'$h_\plus$ $[10^{-22}]$')

#savefig('imr_surr_waveforms.pdf')

tight_layout()
#show()
close()

