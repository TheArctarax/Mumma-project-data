from __future__ import division, print_function
import gwmemory
from matplotlib.pyplot import figure
import numpy as np
from matplotlib.pyplot import *
np.seterr(divide='ignore',invalid='ignore')


fig = figure(figsize=(12, 4))

q = 1.
S1 = [0, 0, 0]
S2 = [0, 0, 0]

colours = ['r', 'b', 'g', 'k']

for ii, model in enumerate(['IMRPhenomD']):
	delta_t = 1.0 / 4096
	times = np.linspace(0, delta_t , 2)
	h_mem, times = gwmemory.time_domain_memory(model=model, q=q, spin_1=S1, spin_2=S2, total_mass=60.,distance=400., inc=np.pi/2, phase=0, times=times)
	plot(times, h_mem['plus'] - h_mem['plus'][np.argmin(abs(times+0.25))],linestyle='-', color=colours[ii], label=model)
	plot(times, h_mem['cross'] - h_mem['cross'][np.argmin(abs(times+0.25))],linestyle='--', color=colours[ii])
        
xlabel('$t (s)$')
ylabel('$\delta h$')
legend(loc='upper left', fontsize=20)

xlim(-0.25, 0.02)

savefig("trial.pdf")
