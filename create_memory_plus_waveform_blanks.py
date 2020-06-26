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

#Q# A. Deine the sampling space by defining times using np.linspace

#Q# B. Now generate h_mem and times using gwmemory.time_domain_memory(...). Don't forget to put in the argument times=times


### These are just for checking the sizes of the h_mem and times array
print(h_mem['plus'].size)
print(times.size)

#Q# C. Now write several lines to plot h_mem['plus'] against times. Set the xlabel to be t(s) and ylabel to be "delta h_{+}". Save the figure as "memory.pdf"



##########
# Create the memory-less waveform 
##########

#Q# D. Now we create hp and hc using get_td_waveform(...) from pycbc. Input the correct arguments into this function. Remember to set the delta_t to be the same as that for the memory. 


## Again, these are just for checking.
print(mock_hp.sample_times[0])
print(mock_hp.sample_times[-1])
print(len(mock_hp))


#Q# E. Now, write several lines to plot the 'plus' memory-less waveform. Set the label names correctly, and save the figure as "original_waveform.pdf".



#########
# Combine the memory with memory-less waveform
#########

## Now, I print both the sizes of hp and h_mem['plus'] to show that they have exactly the same size (so they are add-able). 
print(len(hp), len(h_mem['plus']))

#Q# F. Add hp and h_mem together and define this as strain, i.e. strain equals to hp plus h_mem['plus'] (recall how to add up array element by element. A single line will do the job.)


#Q# G. Finally, we plot (1) strain against times (you should be able to get the time array from an attribute of hp) and (2) hp against time. Label the graphs as "Plus polarized wave (with memory)" and "Plus polarized wave". Set the label names correctly. Finally, save the figure as "combined.pdf"







