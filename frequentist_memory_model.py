from __future__ import division, print_function
from matplotlib.pyplot import *
import numpy as np
import bilby
import gwmemory
from gwmemory import utils as utils
import itertools
from scipy.signal import get_window
from matplotlib.pyplot import rcParams
import pandas as pd

rcParams["font.family"] = "Times New Roman"
rcParams["axes.unicode_minus"] = False
np.seterr(divide="ignore", invalid="ignore")


"""
This code computes the posterior distribution for the memory constant.
Here, we use GW memory to inject a waveform + memory model into noiseless data.
"""


# Set the parameters of the data segment that we're
# going to inject the signal into
duration = 1.0
sampling_frequency = 4096
f_lower = 15.0


# Set up a random seed for result reproducibility. May or may not need this.
np.random.seed(88170235)


def time_domain_window(time_domain_strain, window_type=None):
    if window_type != None:
        window = get_window(window_type, time_domain_strain.size)
        time_domain_strain = time_domain_strain * window

    return time_domain_strain


def memory_time_model(
    times,
    mass_ratio,
    s1x,
    s2x,
    s1y,
    s2y,
    s1z,
    s2z,
    distance,
    total_mass,
    phase,
    inc,
    memory_constant,
):
    # first, we need a linear sample space
    '''
    end_time can only be up to a certain time after merger, which is set in
    geometric units. Conversion from geometric to physical units is given by:
    phys_time = geo_time * (total_mass * m_sun_to_kg)/(c**3/G).
    '''
    GG = 6.674098281543097e-11
    cc = 2.99792458e8
    m_sun_to_kg = 1.98847e+30
    time_lim = 100.00000001920307 # from error message

    start_time = -0.5
    end_time = time_lim * (total_mass * m_sun_to_kg)/(cc**3/GG) - 0.01
    surr_times = np.linspace(
        start_time, end_time, sampling_frequency * (end_time - start_time)
    )

    # Now, to generate an oscillating and secular waveform...
    surr = gwmemory.waveforms.surrogate.Surrogate(
        q=mass_ratio,
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=total_mass,
        distance=distance,
        times=surr_times,
    )
    oscillatory, surr_times = surr.time_domain_oscillatory(inc=inc, phase=phase)
    memory, surr_times = surr.time_domain_memory(inc=inc, phase=phase)

    # ...and add them
    plus_new = oscillatory["plus"] + memory_constant * memory["plus"]
    cross_new = oscillatory["cross"] + memory_constant * memory["cross"]

    # Next, we want to place them in our sample space
    plus = np.zeros(len(times))
    cross = np.zeros(len(times))
    plus[-len(surr_times):] = plus_new
    cross[-len(surr_times):] = cross_new
    
    # Finally, we need to window before applying an fft
    """
    Window types provided by scipy.signal.windows.get_window
    [https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.get_
    window.html#scipy.signal.windows.get_window]
    ------------------------------------------------------------------------------

    'boxcar': rectangular window = no window (essentially)
    'triang': triangular window, nonzero endpoints
    'blackman': 3rd order cosine sum, minimizes leakage, almost as good as Kaiser window at doing so
    'hamming': single cosine with nonzero endpoints, minimizes first side lobe
    'hann': hamming window but with zero endpoints
    'bartlett': triangular window but with zero endpoints, used to taper with little fd modulation
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
    window_type_plus = ("kaiser", 0.1)
    window_type_cross = ("kaiser", 0.1)

    plus = time_domain_window(plus, window_type=window_type_plus)
    cross = time_domain_window(cross, window_type=window_type_cross)

    return {"plus": plus, "cross": cross}


# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(
    total_mass=60.0,
    s1x=0.0,
    s2x=0.0,
    s1y=0.0,
    s2y=0.0,
    s1z=0.0,
    s2z=0.0,
    distance=400,
    mass_ratio=1.0,
    inc=np.pi / 2,
    psi=0.0,
    phase=0.0,
    memory_constant=1.0,
    ra=0.0,
    dec=0.0,
    geocent_time=0.0
)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    time_domain_source_model=memory_time_model,
    start_time=injection_parameters["geocent_time"] - 0.5,
)


# Set up interferometers. In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(["H1", "L1", "V1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 0.5,
)
ifos.inject_signal(
    waveform_generator=waveform, parameters=injection_parameters
)

# Set up a PriorDict, which inherits from dict.
# By default we will sample all terms in the signal models. However, this will
# take a long time for the calculation, so for this example we will set almost
# all of the priors to be equal to their injected values. This implies the
# prior is a delta function at the true, injected value. In reality, the
# sampler implementation is smart enough to not sample any parameter that has
# a delta-function prior.
# The above list does *not* include mass_1, mass_2, theta_jn and luminosity
# distance, which means those are the parameters that will be included in the
# sampler.  If we do nothing, then the default priors get used.
priors = injection_parameters.copy()
priors["memory_constant"] = bilby.core.prior.Uniform(-1, 3, r"$\lambda$")

# Initialise the likelihood by passing in the interferometer data (ifos) and
# the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform
)
likelihood.parameters.update(injection_parameters)

# Plot Frequentist Likelihood
sample_lambda = np.linspace(-1, 3, 11)
likelihood_values = []

for memory_constant in sample_lambda:
    likelihood.parameters['memory_constant'] = memory_constant
    likelihood_values.append(likelihood.log_likelihood())

# Plot Bayesian Likelihood
result_close = bilby.result.read_in_result(filename='/home/darin/bilby_output/test_mid_distance_result.json')
result_mid = bilby.result.read_in_result(filename='/home/darin/bilby_output/test_mid_distance_new_result.json')
result_far = bilby.result.read_in_result(filename='/home/darin/bilby_output/test_far_distance_result.json')

xclose = result_close.posterior['memory_constant']
xmid = result_mid.posterior['memory_constant']
xfar = result_far.posterior['memory_constant']

#yclose = result_close.posterior.log_likelihood.values
yclose=result_close.posterior['log_likelihood']
#ymid = result_mid.posterior.log_likelihood.values
#yfar = result_far.posterior.log_likelihood.values

print(yclose)

plot(xclose, yclose, color='r', label='100 Mpc')
#plot(xclose, yclose, color='b', label='500 Mpc')
#plot(xclose, yclose, color='g', label='1000 Mpc')

xlim(-1, 3)
xlabel(r"memory constant, $\lambda$")
ylabel(r'Likelihood')
legend(loc="upper right", prop={"size": 10})
grid(False)

# plot(sample_lambda, np.exp(likelihood_values))

savefig('/home/darin/bilby_output/three_distances.pdf')
show()
close()

