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
d=200
delta_t=1./4096
f0=15
m1=45
m2=45
plot_begin=-1.5
plot_end=0.1

# GW memory-less definition
hp, hc = get_td_waveform(approximant=apx, mass1=m1, mass2=m2, spin1x=S1[0], spin2x=S2[0], spin1y=S1[1],
                         spin2y=S2[1], spin1z=S1[2], spin2z=S2[2], inclination=inc, coa_phase=pol, distance=d,
                         delta_t=delta_t, f_lower=f0)

start_time=hp.sample_times[0]
end_time=hp.sample_times[-1]

# Check original waveform size
# print(hp.sample_times[0])
# print(hp.sample_times[-1])
# print(len(hp))

# Sample space definition
times = np.linspace(start_time, end_time, len(hp.sample_times))

# GW memory definition
memory, times = gwmemory.time_domain_memory(model=apx, q=m2/m1, total_mass=(m1+m2), distance=d, inc=inc, phase=pol, times=times)

# Check memory array size
# print(memory['plus'].size)
# print(times.size)

# Plot of GW memory
fig = figure(figsize=(12, 6))
plot(times, memory['plus'], color='b')
axhline(0, linestyle=':', color='k')
xlim(-3, 0.1)
xlabel('Time [s]')
ylabel(r'$h_m$$_e$$_m$$_\plus$')

savefig('memory.pdf')

tight_layout()
show()
close()

# Plot of memory-less waveform
fig = figure(figsize=(12, 6))
plot(hp.sample_times, hp, color='g')
axhline(0, linestyle=':', color='k')
xlim(plot_begin, plot_end)
xlabel('Time [s]')
ylabel(r'$h_\plus$')

savefig('original_waveform.pdf')

tight_layout()
show()
close()

# plot of oscillatory + memory components
fig = figure(figsize=(12, 6))

fig.add_subplot(2, 1, 1)
plot(hp.sample_times, hp, color='r')
axhline(0, linestyle=':', color='k')
xlim(plot_begin, plot_end)
ylabel(r'$h_m$$_e$$_m$$_\plus$')

fig.add_subplot(2, 1, 2)
plot(hp.sample_times, hp[:] + memory['plus'][:], color='r')
plot(times, memory['plus'], linestyle='--', color='b')
axhline(0, linestyle=':', color='k')
xlim(plot_begin, plot_end)
ylabel(r'$h_t$$_o$$_t$$_\plus$')
xlabel('Time [s]')

savefig('combined.pdf')

tight_layout()
show()
close()



