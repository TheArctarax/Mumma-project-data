from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

q = 2.
S1 = [0., 0.5, 0.2]
S2 = [0., 0., 0.5]

times = np.linspace(-0.51, 0.06, 10001)
surr = gwmemory.waveforms.surrogate.Surrogate(q=q, spin_1=S1, spin_2=S2, total_mass=60, distance=400, times=times)

inc = np.pi / 2
pol = 0

oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)
memory, times = surr.time_domain_memory(inc=inc, phase=pol)

fig = figure(figsize=(12, 6))
fig.add_subplot(2, 1, 1)
plot(times, oscillatory['plus'], linestyle='--', color='b', alpha=0.5)
plot(times, oscillatory['cross'], linestyle='--', color='r', alpha=0.5)
plot(times, memory['plus'], linestyle='-.', color='b', alpha=0.5)
plot(times, memory['cross'], linestyle='-.', color='r', alpha=0.5)
plot(times, oscillatory['plus'] + memory['plus'], color='b')
plot(times, oscillatory['cross'] + memory['cross'], color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.08, 0.0)

fig.add_subplot(2, 1, 2)
plot(times, oscillatory['plus'], linestyle='--', color='b', alpha=0.5)
plot(times, oscillatory['cross'], linestyle='--', color='r', alpha=0.5)
plot(times, memory['plus'], linestyle='-.', color='b', alpha=0.5)
plot(times, memory['cross'], linestyle='-.', color='r', alpha=0.5)
plot(times, oscillatory['plus'] + memory['plus'], color='b')
plot(times, oscillatory['cross'] + memory['cross'], color='r')
axhline(0, linestyle=':', color='k')
xlim(-0.0, 0.02)

tight_layout()
show()
close()

savefig('oscillatory_and_memory_fig.pdf')
