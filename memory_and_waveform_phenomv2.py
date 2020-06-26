from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import *
from pycbc.waveform import get_td_waveform
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

S1 = [0., 0., 0.]
S2 = [0., 0., 0.]
inc = np.pi / 2
pol = 0
apx = 'IMRPhenomD'
delta_t = 1./4096
d=400
f0=20
m1=30
m2=30
q=m2/m1

times = np.linspace(-0.08, 0.02, 10001)

hp, hc = get_td_waveform(approximant=apx, mass1=m1, mass2=m2, spin1x=S1[0], spin2x=S2[0], spin1y=S1[1],
                         spin2y=S2[1], spin1z=S1[2], spin2z=S2[2], inclination=inc, coa_phase=pol, distance=d,
                         delta_t=delta_t, f_lower=f0)

memory, times = gwmemory.time_domain_memory(model=apx, q=q, total_mass=(m1+m2), distance=d, inc=inc, phase=pol)

fig = figure(figsize=(12, 6))
fig.add_subplot(2, 1, 1)
plot(times, hp['plus'] + memory['plus'], color='b')
axhline(0, linestyle=':', color='k')
xlim(-0.08, 0.02)

tight_layout()
show()
close()

savefig('oscillatory_and_memory_fig2.pdf')

