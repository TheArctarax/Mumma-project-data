from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import figure
import numpy as np

fig = figure(figsize=(12, 4))

q = 1
S1 = [0, 0, 0]
S2 = [0, 0, 0]

colours = ['r', 'b', 'g', 'k']

for ii, model in enumerate(['NRSur7dq2', 'IMRPhenomD', 'SEOBNRv4', 'MWM']):
    h_mem, times = gwmemory.time_domain_memory(q=q, S1=S1, S2=S2, MTot=60., distance=400.,
                                             model=model, inc=np.pi/2, pol=0)

    plot(times, h_mem['plus'] - h_mem['plus'][np.argmin(abs(times+0.25))],
         linestyle='-', color=colours[ii], label=model)
    plot(times, h_mem['cross'] - h_mem['cross'][np.argmin(abs(times+0.25))],
         linestyle='--', color=colours[ii])
        
xlabel('$t (s)$')
ylabel('$\delta h$')
legend(loc='upper left', fontsize=20)

xlim(-0.25, 0.02)

tight_layout()
show()
close()
