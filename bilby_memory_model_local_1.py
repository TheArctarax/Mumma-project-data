from __future__ import division, print_function
import matplotlib
import numpy as np
import bilby
import gwmemory
from gwmemory import utils as utils
import itertools
from scipy.signal import get_window

np.seterr(divide="ignore", invalid="ignore")

"""
When completed, this code will compute the posterior distribution for the
memory constant. Here, we inject a waveform + memory model using GWMemory
into LIGO modelled memory.
"""


# Set the parameters of the data segment that we're
# going to inject the signal into
duration = 1.0
sampling_frequency = 4096
f_lower = 15.0

# Specify the output directory and the name of the simulation.
outdir = "/home/darin/bilby_output"
label = "test_fd_windowing"
bilby.core.utils.setup_logger(outdir=outdir, label=label)


# Set up a random seed for result reproducibility. May or may not need this.
np.random.seed(88170235)


def frequency_domain_transform(time_domain_strain, times, window_type=None):
    sampling_frequency = 1 / (times[1] - times[0])
    frequencies = None

    if window_type != None:
        window = get_window(window_type, time_domain_strain.size)
        time_domain_strain = time_domain_strain * window

    frequency_domain_strain, frequencies = utils.nfft(
        time_domain_strain, sampling_frequency
    )
    return frequency_domain_strain, frequencies


def memory_time_model(
    times,
    q,
    s1x,
    s2x,
    s1y,
    s2y,
    s1z,
    s2z,
    d,
    M,
    psi,
    inc,
    geocent_time,
    memory_constant,
):
    start_time = -0.5
    end_time = 0.0

    surr_times = np.linspace(
        start_time, end_time, sampling_frequency * (end_time - start_time)
    )
    surr = gwmemory.waveforms.surrogate.Surrogate(
        q=q,
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=M,
        distance=d,
        times=surr_times,
    )
    oscillatory, surr_times = surr.time_domain_oscillatory(inc=inc, phase=psi)
    memory, surr_times = surr.time_domain_memory(inc=inc, phase=psi)
    plus = np.zeros(len(times))
    cross = np.zeros(len(times))
    tidx = np.where(times == geocent_time)[0][0]
    print(tidx)
    i = 0
    plus_new = oscillatory["plus"] + memory_constant * memory["plus"]
    cross_new = oscillatory["cross"] + memory_constant * memory["cross"]
    plus_new = np.flip(plus_new)
    cross_new = np.flip(cross_new)
    while i < len(plus_new):
        plus[tidx - i] = plus_new[i]
        cross[tidx - i] = cross_new[i]
        i = i + 1
    return {"plus": plus, "cross": cross}


# Model with memory
def memory_frequency_model(
    f, M, q, s1x, s2x, s1y, s2y, s1z, s2z, d, inc, psi, memory_constant
):
    # Sample space definition for the memory's t-axis. Purposely set to begin,
    # end, and have the same number of points as the original waveform so that
    # superposition of the timeseries is possible.
    start_time = -0.5
    end_time = 0.0
    times = np.linspace(
        start_time, end_time, sampling_frequency * (end_time - start_time)
    )

    # GW waveform with memory definition
    # The sub-function waveforms.surrogate.Surrogate generates a surrogate
    # object.
    surr = gwmemory.waveforms.surrogate.Surrogate(
        q=q,
        name="nrsur7dq2",
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=M,
        distance=d,
        times=times,
    )

    # GW waveform only definition
    # A surrogate object has the following attributes: time_domain_memory (a
    # pycbc.timeseries object that has both the ['plus'] and ['cross']
    # sub-arrays), and time_domain_memory which also has ['plus'] and
    # ['cross']). Calling these attributes returns both the pycbc timeseries
    # and the sampling time arrays (which you store it as times here).
    oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=psi)

    # GW memory definition
    memory, times = surr.time_domain_memory(inc=inc, phase=psi)

    """
    Window Types provided by scipy.signal.windows.get_window
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

    # waveform
    plus = oscillatory["plus"][:] + memory_constant * memory["plus"][:]
    cross = oscillatory["cross"][:] + memory_constant * memory["cross"][:]
    plus_tilde_new, frequencies = frequency_domain_transform(
        plus, times, window_type=window_type_plus
    )
    cross_tilde_new, frequencies = frequency_domain_transform(
        cross, times, window_type=window_type_cross
    )
    plus_tilde = np.zeros(len(f))
    cross_tilde = np.zeros(len(f))
    plus_tilde[: len(plus_tilde_new)] = plus_tilde_new[:]
    cross_tilde[: len(cross_tilde_new)] = cross_tilde_new[:]
    print(len(plus_tilde), len(cross_tilde))
    return {"plus": plus_tilde, "cross": cross_tilde}


# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(
    M=60.0,
    s1x=0.0,
    s2x=0.0,
    s1y=0.0,
    s2y=0.0,
    s1z=0.0,
    s2z=0.0,
    d=400,
    q=1.0,
    inc=np.pi / 2,
    psi=0.0,
    memory_constant=1.0,
    ra=0.0,
    dec=0.0,
    geocent_time=0.0,
)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=memory_frequency_model,
    start_time=injection_parameters["geocent_time"] + start_time,
)


# Set up interferometers. In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(["H1", "L1", "V1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] + start_time,
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

# Run sampler.  In this case we're going to use the `dynesty` sampler
result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler="dynesty",
    use_ratio=True,
    plot=True,
    npoints=100,
    sample="rwalk",
    verbose=True,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
)

# Make a corner plot.
result.plot_corner()
