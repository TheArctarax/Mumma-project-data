from __future__ import division, print_function
from scipy.signal import get_window
import gwmemory
import matplotlib
from gwmemory import utils as utils
from matplotlib.pyplot import *
import numpy as np
import pandas as pd
import wget
import itertools
from matplotlib.pyplot import rcParams

rcParams["font.family"] = "Times New Roman"
rcParams["axes.unicode_minus"] = False
np.seterr(divide="ignore", invalid="ignore")

"""
when finished, this program will calculate and plot asd's for waveforms and their
associated memory. Comparison is then made with noise asd's by superposing plots.
Noise asd's come from GWOSC.
"""


# GWmemory natively features a frequency domain converter for memory but not for waveforms.
# This method fixes that problem by following the same structure as frequency_domain_memory().

# Q1 : Start by looking at the frequency_domain_memory definition you quoted. Here, we want
#      to define a new function called `frequency_domain_transform`, which takes in two
#      arguments - time_domain_waveform (which is basically either the memory, or the oscillatory)
#      and times (the sampling times array), and do a fourier transform to the frequency
#      domain, and finally, returns the frequency_domain_waveform and sampling frequencies
#      to the user (so two products). Modify the code in frequency_domain_memory to attain
#      the above goals.
def frequency_domain_transform(time_domain_strain, times, window_type=None):
    sampling_frequency = 1 / (times[1] - times[0])
    frequencies = None
    frequency_domain_strain = dict()

    for key in time_domain_strain:
        if window_type != None:
            window = get_window(window_type, time_domain_strain[key].size)
            time_domain_strain[key] = time_domain_strain[key] * window
        frequency_domain_strain[key], frequencies = utils.nfft(
            time_domain_strain[key], sampling_frequency
        )

    return frequency_domain_strain, frequencies


# Variable assignments
S1 = [0.0, 0.0, 0.0]
S2 = [0.0, 0.0, 0.0]
inc = np.pi / 2
pol = 0
# d=420
# M=66.2
d = 20
M = 80
q = 1.16
osc_window_type = ("tukey", 0.1)
mem_window_type = ("kaiser", 0.001)

"""
Window Types provided by scipy.signal.windows.get_window
[https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.get_window.html#scipy.signal.windows.get_window]
------------------------------------------------------------------------------

'boxcar': rectangular window = no window (essentially)
'triang': triangular window, nonzero endpoints
'blackman': 3rd order cosine sum, minimizes leakage, almost as good as Kaiser window at doing so
'hamming': single cosine with nonzero endpoints, minimizes first side lobe
'hann': hamming window but with zero endpoints
'bartlett': triangular window but with zero enpoints, used to taper with little fd modulation
'flattop': 5th order cosine sum, used to measure signal amplitude, makes main lobe flat
'parzen': not sure about this one
'bohman': or this one, either
'blackmanharris': generalized hamming = more cosines, hamming but better
'nuttall': similar to blackman-harris
'barthann': combo of bartlett and hann
('kaiser', beta): formed from Bessel functions, beta=0(rect), 5(hamming), 6(hann), 8.6(blackman)
('gaussian', std_dev): use only in special cases
('general_gaussian', power, width): same here
('slepian', width): maximizes power in main lobe
('dpss', norm half-bandwidth): first term is slepian window
('chebwin', attenuation): uses Chebyshev polynomials, kinda complicated
('exponential', decay constant): seems like it will cut power too quickly
('tukey', taper fraction): tf=0(rect), 1(hann)

------------------------------------------------------------------------------
"""

# Sample space definition for the memory's t-axis.
start_time = -0.51
end_time = 0.02
times = np.linspace(start_time, end_time, 10001)


# GW waveform with memory definition
# The sub-function waveforms.surrogate.Surrogate generates a surrogate object.
surr = gwmemory.waveforms.surrogate.Surrogate(
    q=q,
    name="nrsur7dq2",
    spin_1=S1,
    spin_2=S2,
    total_mass=M,
    distance=d,
    times=times,
)


# GW waveform-only definition
# A surrogate object has the following attributes: frequency_domain_memory (a
# pycbc.frequencyseries object that has both the ['plus'] and ['cross'] sub-arrays),
# and frequency_domain_memory which also has ['plus'] and ['cross']). Calling
# these attributes returns both the pycbc frequency series and the sampling
# frequencies.

# Now, basically you can still do the same (using time_domain_oscillatory) here.
oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)

# Q2 : Now it's time to apply our newly defined function `frequency_domain_transform`. Use
#      this to transform the osciallory time domain waveform to the frequency domain. We
#      use oscillatory_tilde, frequencies to store the products. Complete the following
#      line of code.
oscillatory_tilde, frequencies = frequency_domain_transform(
    oscillatory, times, window_type=osc_window_type
)


# GW memory definition
# Same thing happens for the memory
memory, times = surr.time_domain_memory(inc=inc, phase=pol)

# Q3 : Again, apply our new function to transform the memory into the frequency domain.
memory_tilde, frequencies = frequency_domain_transform(
    memory, times, mem_window_type
)


# Plot of the frequency domain waveform.

# Q4 : It's usually best to plot these figures in a log-log scale. Figure how to do so
#      and modify the following line. Note that you also have to change the arguments:
#      In particular, you should plot frequencies as the x-axis, and the **magnitude**
#      of the amplitude_tilde times square root of frequencies as the y-axis.
fig = figure(figsize=(6, 6))
loglog(
    frequencies,
    abs(oscillatory_tilde["plus"][:] * np.sqrt(frequencies[:])),
    color="r",
)

xlim(0, 5000)
xlabel("Frequency [Hz]")
ylabel(r"ASD $[Hz$^{-1/2}$]$")
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)
rc("axes", labelsize=14)

savefig("waveform_asd_waveform.pdf")

tight_layout()
# show()
close()


# Plot of the frequency domain memory.

# Q5 : Do the same for the memory part.
fig = figure(figsize=(6, 6))
loglog(
    frequencies,
    abs(memory_tilde["plus"][:] * np.sqrt(frequencies[:])),
    color="r",
)

xlim(0, 5000)
xlabel("Frequency [Hz]")
ylabel(r"ASD $[Hz$^{-1/2}$]$")
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)
rc("axes", labelsize=14)

savefig("waveform_asd_mem.pdf")

tight_layout()
# show()
close()


# Plot of the frequency domain waveform + memory.

# Q6 : And of course, again for the waveform + memory.
fig = figure(figsize=(6, 6))
loglog(
    frequencies,
    abs(
        (oscillatory_tilde["plus"][:] + memory_tilde["plus"][:])
        * np.sqrt(frequencies)
    ),
    color="r",
)

xlim(0, 5000)
xlabel("Frequency [Hz]")
ylabel(r"ASD $[Hz$^{-1/2}$]$")
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)
rc("axes", labelsize=14)

savefig("waveform_asd_wave_and_mem.pdf")

tight_layout()
# show()
close()


# Plot of the frequency domain waveform, memory, and noise curves.

# Q7 : Finally, do the same for the full plot.
fig = figure(figsize=(6, 6))
loglog(
    frequencies,
    abs(
        (oscillatory_tilde["plus"][:] + memory_tilde["plus"][:])
        * np.sqrt(frequencies)
    ),
    linewidth=8,
    color="tab:orange",
    zorder=1,
)
# loglog(frequencies, abs(oscillatory_tilde['plus']*np.sqrt(frequencies)), color='g', zorder=1)
loglog(
    frequencies,
    abs(memory_tilde["plus"] * np.sqrt(frequencies)),
    color="b",
    linewidth=6,
    zorder=1,
)

# Here, I would like to call noise_curve() for superposition with the memory and waveform curves.
# Q8 : Finally, to plot the noise curve, simply load the txt files and plot them again here.
dfl = pd.read_csv(
    "L1_O2_Sensitivity_strain_asd.txt",
    sep="\t",
    index_col=False,
    header=None,
    names=["frequency", "asd"],
)
dfh = pd.read_csv(
    "H1_O2_Sensitivity_strain_asd.txt",
    sep="\t",
    index_col=False,
    header=None,
    names=["frequency", "asd"],
)
dfv = pd.read_csv(
    "V1_O2_Sensitivity_strain_asd.txt",
    sep="\s+",
    index_col=False,
    header=None,
    names=["frequency", "asd"],
)

loglog(dfl["frequency"], dfl["asd"], color="tab:blue", zorder=0, label="L1")
loglog(dfh["frequency"], dfh["asd"], color="tab:red", zorder=0, label="H1")
loglog(dfv["frequency"], dfv["asd"], color="tab:purple", zorder=0, label="V1")

xlim(10, 1000)
# ylim(10**(-24), 10**(-21))
xlabel("Frequency [Hz]")
ylabel(r"ASD $[Hz$^{-1/2}$]$")
xscale("log")
yscale("log")
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)
rc("axes", labelsize=14)

savefig("waveform_asd_wave_and_mem_and_noise.pdf")

tight_layout()
show()
close()
