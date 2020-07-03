from __future__ import division, print_function
from pycbc.waveform import get_fd_waveform
import numpy as np
import bilby
import gwmemory
np.seterr(divide='ignore', invalid='ignore')

'''
When completed, this code will compute the posterior distribution for the memory constant.
Here, we inject a waveform + memory model using GWMemory into LIGO modelled memory.
'''


# Set the parameters of the data segment that we're
# going to inject the signal into
duration = 4.
sampling_frequency = 2500.
f_lower=15

# Specify the output directory and the name of the simulation !!!need to change!!!.
outdir = '/home/darin/bilby_output'
label = 'test1'
bilby.core.utils.setup_logger(outdir=outdir, label=label)


# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# Model with memory 
def memory_model(M, S1, S2, d, inc, pol, memory_constant):
  # Sample space definition for the memory's t-axis. Purposely set to begin, end, and have the same number of points as the
  # original waveform so that superposition of the timeseries is possible.
  start_time=-10
  end_time=0.02
  times = np.linspace(start_time, end_time, sampling_frequency*duration)


  # GW waveform with memory definition
  # The sub-function waveforms.surrogate.Surrogate generates a surrogate object.
  surr = gwmemory.waveforms.surrogate.Surrogate(q=1, name='nrsur7dq2', spin_1=S1,
                                                spin_2=S2, total_mass=M,
                                                distance=d, times=times)


  # GW waveform only definition
  # A surrogate object has the following attributes: time_domain_memory (a 
  # pycbc.timeseries object that has both the ['plus'] and ['cross'] sub-arrays), 
  # and time_domain_memory which also has ['plus'] and ['cross']). Calling
  # these attributes returns both the pycbc timesseries and the sampling time
  # arrays (which you store it as times here).
  oscillatory, times = surr.time_domain_oscillatory(inc=inc, phase=pol)

  # GW memory definition
  memory, times = surr.time_domain_memory(inc=inc, phase=pol)

  # waveform
  plus = oscillitory['plus'][:] + memory_constant*memory['plus'][:]
  cross = oscillitory['cross'][:] + memory_constant*memory['cross'][:]

  return {plus, cross}

# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(M=60., S1=[0., 0., 0.], S2=[0., 0., 0.],
                            d=600., inc=np.pi/2, pol=0., geocent_time=0.,
                            memory_constant=1.)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    time_domain_source_model=memory_model)

# Set up interferometers.  In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(['H1', 'L1', 'V1'])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)

# Set up a PriorDict, which inherits from dict.
# By default we will sample all terms in the signal models.  However, this will
# take a long time for the calculation, so for this example we will set almost
# all of the priors to be equal to their injected values.  This implies the
# prior is a delta function at the true, injected value.  In reality, the
# sampler implementation is smart enough to not sample any parameter that has
# a delta-function prior.
# The above list does *not* include mass_1, mass_2, theta_jn and luminosity
# distance, which means those are the parameters that will be included in the
# sampler.  If we do nothing, then the default priors get used.
priors = injection_parameters.copy()
priors['memory_constant'] = bilby.core.prior.Uniform(0, 1, r'$\lambda$')

# Initialise the likelihood by passing in the interferometer data (ifos) and
# the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator)

# Run sampler.  In this case we're going to use the `dynesty` sampler
result = bilby.run_sampler(
    likelihood=likelihood, priors=priors, sampler='dynesty', npoints=500,
    injection_parameters=injection_parameters, outdir=outdir, label=label)

# Make a corner plot.
result.plot_corner()
