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
apx = 'IMRPhenomD'
d=600
delta_t=1./4096
f0=15
m1=30
m2=30


# GW memory-less definition
hp, hc = get_td_waveform(approximant=apx, mass1=m1, mass2=m2, spin1x=S1[0], spin2x=S2[0], spin1y=S1[1],
                         spin2y=S2[1], spin1z=S1[2], spin2z=S2[2], inclination=inc, coa_phase=pol, distance=d,
                         delta_t=delta_t, f_lower=f0)


# Sample space definition for the memory's t-axis. Purposely set to begin, end, and have the same number of points as the
# original waveform so that superposition of the timeseries is possible.
start_time=hp.sample_times[0]
end_time=hp.sample_times[-1]
times = np.linspace(start_time, end_time, len(hp.sample_times))


# GW memory definition
memory, times = gwmemory.time_domain_memory(model=apx, q=m2/m1, total_mass=(m1+m2), distance=d, inc=inc, phase=pol, times=times)


# Plot of GW memory
fig = figure(figsize=(12, 6))
plot(times, memory['plus']*(10**22), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.05)
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
plot(hp.sample_times, hp*(10**22), color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.5, 0.05)
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
plot(hp.sample_times, (hp[:] + memory['plus'][:])*(10**22), color='r', label=r'Waveform $\plus$ Memory')
plot(hp.sample_times, hp*(10**22), linestyle='--' , color='tab:purple', label='Original Waveform')
axhline(0, linestyle=':', color='k')
xlim(-0.04, 0.015)
xlabel('Time (s)')
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)


fig.add_subplot(2, 1, 1)
plot(hp.sample_times, (hp[:] + memory['plus'][:])*(10**22), color='r', label=r'Waveform $\plus$ Memory')
plot(times, memory['plus']*(10**22), linestyle='--', color='b', label='Memory')
axhline(0, linestyle=':', color='k')
xlim(-0.3, 0.02)
ylabel(r'$h_\plus$ $[10^{-22}]$')
legend(loc='upper left')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

savefig('combined.pdf')

tight_layout()
show()
close()



