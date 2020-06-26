from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
from pycbc.waveform import get_td_waveform
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

q = 1.
S1 = [0., 0., 0.]
S2 = [0., 0., 0.]

# times = np.linspace(-0.08, 0.02, 10001)
# surr = gwmemory.waveforms.surrogate.Surrogate(q=q, spin_1=S1, spin_2=S2, total_mass=60, distance=400, times=times)

times = np.linspace(-0.08, 0.02, 10001)
hp, hc = get_td_waveform(approximant='IMRPhenomD', mass1=30, mass2=30, spin1x=S1[0], spin2x=S2[0], spin1y=S1[1], spin2y=S2[1], spin1z=S1[2], spin2z=S2[2], distance=400, delta_t=1./4096, f_lower=20)

inc = np.pi / 2
pol = 0

oscillatory, times = hp.time_domain_oscillatory(inc=inc, phase=pol)
memory, times = hp.time_domain_memory(inc=inc, phase=pol)

fig = figure(figsize=(12, 6))
fig.add_subplot(2, 1, 1)
plot(times, oscillatory['plus'] + memory['plus'], color='b')
axhline(0, linestyle=':', color='k')
xlim(-0.08, 0.02)

tight_layout()
show()
close()

savefig('oscillatory_and_memory_fig2.pdf')

