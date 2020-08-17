#!/home/darin.mumma/surf2020_env/bin/python

from __future__ import division, print_function
import matplotlib
import numpy as np
import bilby
import gwmemory
from gwmemory import utils as utils
import itertools
from scipy.signal import get_window
import argparse
from gwpy.timeseries import TimeSeries

np.seterr(divide="ignore", invalid="ignore")

"""
This code computes the posterior distribution for an n-dimensional parameter
space which includes the memory constant. Here, we use GWMemory to inject a
waveform + memory model into REAL data.
"""

# default values are set to GW150914's parameters
def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', help='output directory')
    parser.add_argument('--label', help='file name without extension')
    parser.add_argument('--m',
                        '--total_mass',
                        help='total mass of CBC source',
                        default=70.436052,
    )
    parser.add_argument('--s1x',
                        help='x-projected spin of first binary component',
                        default=0.,
    )
    parser.add_argument('--s2x',
                        help='x-projected spin of second binary component',
                        default=0.,
    )
    parser.add_argument('--s1y',
                        help='y-projected spin of first binary component',
                        default=0.,
    )
    parser.add_argument('--s2y',
                        help='y-projected spin of second binary component',
                        default=0.,
    )
    parser.add_argument('--s1z',
                        help='z-projected spin of first binary component',
                        default=0.10151,
    )
    parser.add_argument('--s2z',
                        help='z-projected spin of second binary component',
                        default=-0.216688,
    )
    parser.add_argument('--d',
                        '--distance',
                        help='luminosity distance to CBC source',
                        default=342.21115,
    )
    parser.add_argument('--q',
                        '--mass_ratio',
                        help='mass ratio (m1/m2) of CBC source',
                        default=0.9003307,
    )
    parser.add_argument('--i',
                        '--inclination',
                        help='inclination of CBC source (0 = face-on, np.pi/2 = edge-on)',
                        default=2.469643,
    )
    parser.add_argument('--psi',
                        '--polarization_angle',
                        help='gravitational wave polarization angle (0 <= psi <= np.pi)',
                        default=0.035035,
    )
    parser.add_argument('--phase',
                        help='gravitational wave phase (0 <= phase <= 2*np.pi)',
                        default=1.972995,
    )
    parser.add_argument('--mc',
                        '--memory_constant',
                        help='scaling factor for the gw memory term',
                        default=1.,
    )
    parser.add_argument('--ra',
                        help='right ascension of CBC source',
                        default=1.157876,
    )
    parser.add_argument('--dec',
                        help='declination of CBC source',
                        default=-1.19108,
    )
    parser.add_argument('--t',
                        '--geocent_time',
                        help='time of merger (max signal amplitude); usually 0',
                        default=1.126259*(10.0**9.0),
    )

    options = parser.parse_args()
    
    if (options.outdir == None):
        raise IOError('You forgot to specify the output directory')

    return options

options = parse_command_line()

# Set the parameters of the data segment that we're
# going to inject the signal into
duration = 1.0
sampling_frequency = 4096
f_lower = 15.0

# Specify the output directory and the name of the simulation.
outdir = options.outdir
label = options.label
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility. May or may not need this.
np.random.seed(88170235)

# Returns a two-dimensional array with lower and upper time bounds as elements. This is done by creates sur object of equal specification to the desired signal and extracting its get_t_lim attribute.
def get_t_0_t_f(
    mass_ratio,
    s1x,
    s1y,
    s1z,
    s2x,
    s2y,
    s2z,
    distance,
    total_mass,
):

    start_time = -0.5  # arbitrarily chosen
    end_time = 0.0  # also arbitrary
    test_surr_times = np.linspace(
        start_time, end_time, sampling_frequency * (end_time - start_time)
    )

    # Now, create a toy model (only total mass matters) from which we can retrieve time bounds
    test_surr = gwmemory.waveforms.surrogate.Surrogate(
        q=mass_ratio,
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=total_mass,
        distance=distance,
        times=test_surr_times,
    )

    new_test_surr_times = test_surr_times / test_surr.t_to_geo
    return test_surr.sur.find_t_0(
        test_surr.q,
        test_surr.S1,
        test_surr.S2,
        MTot=test_surr.MTot,
        distance=test_surr.distance,
        t=new_test_surr_times,
        LMax=test_surr.LMax,
    )

def memory_time_model(
    times,
    distance,
    phase,
    inc,
    memory_constant,
):    
    # First, to generate an oscillating and secular waveform...
    oscillatory = gwmemory.utils.combine_modes(h_lm, inc, phase)
    memory = gwmemory.utils.combine_modes(h_lm_mem, inc, phase)

    # ...and add them (also scale for distance)
    plus_new = (oscillatory["plus"] + memory_constant * memory["plus"]) / distance
    cross_new = (oscillatory["cross"] + memory_constant * memory["cross"]) / distance

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
    # Next, we want to place them in our sample space
    plus = np.zeros(len(times))
    cross = np.zeros(len(times))
    plus[-len(surr_times) :] = plus_new * window
    cross[-len(surr_times) :] = cross_new * window

    return {"plus": plus, "cross": cross}

# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(
    mass_ratio=float(options.q),
    s1x=float(options.s1x),
    s2x=float(options.s2x),
    s1y=float(options.s1y),
    s2y=float(options.s2y),
    s1z=float(options.s1z),
    s2z=float(options.s2z),
    total_mass=float(options.m),
    distance=float(options.d),
    inc=float(options.i),
    psi=float(options.psi),
    phase=float(options.phase),
    memory_constant=float(options.mc),
    ra=float(options.ra),
    dec=float(options.dec),
    geocent_time=float(options.t),
)

extrinsic_injection_parameters = dict(
    distance=float(options.d),
    inc=float(options.i),
    phase=float(options.phase),
    memory_constant=float(options.mc),
)

# retrieves valid template interval
time_lim = get_t_0_t_f(
    mass_ratio=injection_parameters["mass_ratio"],
    s1x=injection_parameters["s1x"],
    s2x=injection_parameters["s2x"],
    s1y=injection_parameters["s1y"],
    s2y=injection_parameters["s2y"],
    s1z=injection_parameters["s1z"],
    s2z=injection_parameters["s2z"],
    distance=injection_parameters["distance"],
    total_mass=injection_parameters["total_mass"],
)
# We need a linear sample space
"""
end_time can only be up to a certain time after merger, which is set in
geometric units. Conversion from geometric to physical units is given by:
phys_time = geo_time * (total_mass * m_sun_to_kg)/(c**3/G).
"""
GG = 6.674098281543097e-11
cc = 2.99792458e8
m_sun_to_kg = 1.98847e30
t_f = time_lim[1] + 0.9  # Surrogate class cuts bound by 1.0s already
start_time = -0.5
end_time = t_f * (injection_parameters['total_mass'] * m_sun_to_kg) / (cc ** 3 / GG)
surr_times = np.linspace(
    start_time, end_time, sampling_frequency * (end_time - start_time)
)
# We want to create a surrogate object
surr = gwmemory.waveforms.surrogate.Surrogate(
    q=injection_parameters["mass_ratio"],
    spin_1=[injection_parameters["s1x"], injection_parameters["s1y"], injection_parameters["s1z"]],
    spin_2=[injection_parameters["s2x"], injection_parameters["s2y"], injection_parameters["s2z"]],
    total_mass=injection_parameters["total_mass"],
    distance=1.0,
    times=surr_times,
    # modes=[(2,2),(2,-2)],
)

h_lm, surr_times = surr.time_domain_oscillatory()
h_lm_mem, surr_times = surr.time_domain_memory()

window_type = ("kaiser", 0.1)
window = get_window(window_type, surr_times.size)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    time_domain_source_model=memory_time_model,
)

times = waveform.time_array

# Set up interferometers. In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
logger = bilby.core.utils.logger
trigger_time = float(options.t)
roll_off = 0.4  # Roll off duration of tukey window in seconds, default is 0.4s
duration = 4  # Analysis segment duration
post_trigger_duration = 2  # Time between trigger time and end of segment
end_time = trigger_time + post_trigger_duration
start_time = end_time - duration

psd_duration = 32 * duration
psd_start_time = start_time - psd_duration
psd_end_time = start_time

# We now use gwpy to obtain analysis and psd data and create the ifo_list
ifo_list = bilby.gw.detector.InterferometerList([])
for det in ["H1", "L1"]:
    logger.info("Downloading analysis data for ifo {}".format(det))
    ifo = bilby.gw.detector.get_empty_interferometer(det)
    data = TimeSeries.fetch_open_data(det, start_time, end_time)
    ifo.strain_data.set_from_gwpy_timeseries(data)

    logger.info("Downloading psd data for ifo {}".format(det))
    psd_data = TimeSeries.fetch_open_data(det, psd_start_time, psd_end_time)
    psd_alpha = 2 * roll_off / duration
    psd = psd_data.psd(
        fftlength=duration,
        overlap=0,
        window=("tukey", psd_alpha),
        method="median"
    )
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=psd.frequencies.value, psd_array=psd.value)
    ifo_list.append(ifo)

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
priors["memory_constant"] = bilby.core.prior.Uniform(-5, 5, r"$\lambda$")
# priors["distance"] = bilby.core.prior.Uniform(80, 120, r"$d_L$")
priors["psi"] = bilby.core.prior.Uniform(0.0, np.pi, r"$\psi$")
priors["phase"] = bilby.core.prior.Uniform(0.0, 2.0 * np.pi, r"$\phi$")
priors["inc"] = bilby.core.prior.Sine(
    name=r"$\iota$",
    latex_label=r"$\iota$",
    unit=None,
    minimum=0.0,
    maximum=np.pi
)

# Initialise the likelihood by passing in the interferometer data (ifos) and
# the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifo_list, waveform_generator=waveform
)

# Run sampler.  In this case we're going to use the `dynesty` sampler
result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler="dynesty",
    use_ratio=True,
    plot=True,
    npoints=1000,
    sample="rwalk",
    walks=20,
    verbose=True,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
)

# Make a corner plot.
result.plot_corner()
