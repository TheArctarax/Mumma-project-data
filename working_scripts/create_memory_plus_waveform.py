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

##########
# Create the memory waveform ##
##########

## A. Deine the sampling space by defining times using np.linspace
times = np.linspace(0,1.0/4096.0 , 2)

h_mem, times = gwmemory.time_domain_memory(model='IMRPhenomD', q=1, total_mass=60,distance=400., inc=np.pi/2, phase=0, times=times)

print(h_mem['plus'].size)
print(times.size)

plot(times, h_mem['plus'])

xlabel('$t (s)$')
ylabel('$\delta h_{+}$')

savefig("memory.pdf")
show()
close()

##########
# Create the memory-less waveform 
##########

hp, hc = get_td_waveform(approximant="IMRPhenomD", mass1=30, mass2 = 30, delta_t = 1.0/4096.0, f_lower= 15.0, inclination=np.pi/2)
print(hp.sample_times[0])
print(hp.sample_times[-1])

print(len(hp))
plot(hp.sample_times, hp, label='Plus Polarization')
xlabel('Time (s)')
legend()
grid()
savefig("original_waveform.pdf")
show()
close()


#########
# Combine the memory with memory-less waveform
#########

strain = hp[:] + h_mem['plus'][:]

plot(hp.sample_times, hp, label='Plus Polarization')
plot(hp.sample_times, strain, label='Plus Polarization (with memory)')
xlabel('Time (s)')
legend()
grid()
savefig("combined.pdf")
show()
close()

