from __future__ import division, print_function
import matplotlib
matplotlib.use("Agg")
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
duration = 4.
sampling_frequency = 4096
f_lower = 15.0

# Specify the output directory and the name of the simulation.
outdir = "/Users/alvinli/bilby_output"
label = "test1"
bilby.core.utils.setup_logger(outdir=outdir, label=label)


# Set up a random seed for result reproducibility. May or may not need this.
np.random.seed(88170235)

def frequency_domain_transform(time_domain_strain, times):
  sampling_frequency = 1 / (times[1] - times[0])
  frequencies = None
  frequency_domain_strain, frequencies = utils.nfft(time_domain_strain, sampling_frequency)
  return frequency_domain_strain, frequencies

def memory_time_model(times, q, s1x, s2x, s1y, s2y, s1z, s2z, d, M, psi, inc, geocent_time, memory_constant):
    surr_times = np.linspace(-0.5, 0.0, 4096*0.5)
    surr = gwmemory.waveforms.surrogate.Surrogate(q=q, spin_1=[s1x, s1y, s1z], spin_2=[s2x, s2y, s2z], total_mass=M, distance=d, times=surr_times)
    oscillatory, surr_times = surr.time_domain_oscillatory(inc=inc, phase=psi)
    memory, surr_times = surr.time_domain_memory(inc=inc, phase=psi)
    plus = np.zeros(len(times))
    cross = np.zeros(len(times))
    tidx = np.where(times==geocent_time)[0][0]
    print(tidx)
    i=0
    plus_new = oscillatory['plus'] + memory['plus']
    cross_new = oscillatory['cross'] + memory['cross']
    plus_new = np.flip(plus_new)
    cross_new = np.flip(cross_new)
    while i<len(plus_new):
       plus[tidx-i] = plus_new[i]
       cross[tidx-i] = cross_new[i]
       i=i+1
    return {'plus': plus, 'cross': cross}	

# Model with memory
def memory_model(
    M, s1x, s2x, s1y, s2y, s1z, s2z, d, inc, psi, memory_constant
):
    # Sample space definition for the memory's t-axis. Purposely set to begin,
    # end, and have the same number of points as the original waveform so that
    # superposition of the timeseries is possible.
    start_time = -0.08
    end_time = 0.02
    times = np.linspace(start_time, end_time, 4097)

    # GW waveform with memory definition
    # The sub-function waveforms.surrogate.Surrogate generates a surrogate
    # object.
    surr = gwmemory.waveforms.surrogate.Surrogate(
        q=1,
        name="nrsur7dq2",
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=M,
        distance=d,
        times=times
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

    # waveform
    plus = oscillatory["plus"][:] + memory_constant * memory["plus"][:]
    cross = oscillatory["cross"][:] + memory_constant * memory["cross"][:]
    plus_tilde, frequencies = frequency_domain_transform(plus, times)
    cross_tilde, frequencies = frequency_domain_transform(cross, times)
    print(len(plus_tilde), len(cross_tilde))
    return {"plus": plus_tilde, "cross": cross_tilde}


# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(
    M=100.0,
    s1x=0.0,
    s2x=0.0,
    s1y=0.0,
    s2y=0.0,
    s1z=0.0,
    s2z=0.0,
    d=50.0,
    q=1.,
    inc=np.pi / 2,
    psi=0.0,
    memory_constant=1.0,
    ra=0.0,
    dec=0.0,
    geocent_time=0.0,
)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    time_domain_source_model=memory_time_model,
    start_time=injection_parameters['geocent_time'] - 0.5)


# Set up interferometers. In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(["H1", "L1"])
ifos.set_strain_data_from_power_spectral_densities(
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
priors["memory_constant"] = bilby.core.prior.Uniform(-1, 1, r"$\lambda$")

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
    npoints=1000,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
)

# Make a corner plot.
result.plot_corner()
