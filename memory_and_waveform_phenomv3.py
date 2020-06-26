from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
from pycbc.waveform import get_td_waveform
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

# Variable assignments
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0
apx = 'IMRPhenomD'
d=400
delta_t=1./4096
f0=20
m1=30
m2=30
N=10001

# GW memory-less definition
hp, hc = get_td_waveform(approximant=apx, mass1=m1, mass2=m2, spin1x=S1[0], spin2x=S2[0], spin1y=S1[1],
                         spin2y=S2[1], spin1z=S1[2], spin2z=S2[2], inclination=inc, coa_phase=pol, distance=d,
                         delta_t=delta_t, f_lower=f0)

start_time=hp.sample_times[0]
end_time=hp.sample_times[-1]

# Check original waveform size
print(hp.sample_times[0])
print(hp.sample_times[-1])
print(len(hp))

# Sample space definition
times = np.linspace(0,2*delta_t, 2)

# GW memory definition
memory, times = gwmemory.time_domain_memory(model=apx, q=m2/m1, total_mass=(m1+m2), distance=d, inc=inc, phase=pol, times=times)

# Check memory array size
print(memory['plus'].size)
print(times.size)

# Plot of GW memory
fig = figure(figsize=(12, 6))
plot(times, memory['plus'], color='b')
axhline(0, linestyle=':', color='k')
xlim(start_time, end_time)
xlabel('Time [s]')
ylabel(r'$h_m$$_e$$_m$$_\plus$')

tight_layout()
show()

savefig('memory.pdf')

# Plot of memory-less waveform
fig = figure(figsize=(12, 6))
plot(times, hp, color='g')
axhline(0, linestyle=':', color='k')
xlim(start_time, end_time)
xlabel('Time [s]')
ylabel(r'$h_\plus$')

tight_layout()
show()

savefig('original_waveform.pdf')


# plot of oscillatory + memory components
fig = figure(figsize=(12, 6))

fig.add_subplot(2, 1, 1)
plot(times, hp, color='r')
axhline(0, linestyle=':', color='k')
xlim(start_time, end_time)
ylabel(r'$h_m$$_e$$_m$ $_\plus$')

fig.add_subplot(2, 1, 2)
plot(times, hp + memory['plus'], color='r')
axhline(0, linestyle=':', color='k')
xlim(start_time, end_time)
ylabel(r'$h_t$$_o$$_t$ $_\plus$')
xlabel('Time [s]')

tight_layout()
show()
close()

savefig('combined.pdf')

