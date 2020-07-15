from __future__ import division, print_function
import matplotlib
matplotlib.use("agg")
import numpy as np
import bilby
import matplotlib.pyplot as plt
import gwmemory
from gwmemory import utils as utils
import itertools
import scipy.optimize as opt
from scipy.signal import get_window

#np.seterr(divide="ignore", invalid="ignore")


'''
When completed, this code will handle 1- and 2-D GW parameter estimation by
plotting the likelihood and identifying each parameters' value at the peak.

Outline
------------------------------------------------------------------------------
i. define linspace for memory constant
ii. define model (pull from bilby code)
iii. define data (also pull from bilby code)
iv. compute likelihood (need likelihood as a function of memory constant)
v. plot likelihood (easy)
vi. find values at peak (also easy)
------------------------------------------------------------------------------
'''

# i.
sampling_frequency = 4096
a = -1 # left endpoint
b = 2 # right endpoint
duration = 1.0
sample_lambda = np.linspace(a, b, sampling_frequency * (b-a))

# ii.
def time_domain_window(time_domain_strain, window_type=None):
    if window_type != None:
        window = get_window(window_type, time_domain_strain.size)
        time_domain_strain = time_domain_strain * window

    return time_domain_strain

def get_t_0_t_f(
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

    GG = 6.674098281543097e-11
    cc = 2.99792458e8
    m_sun_to_kg = 1.98847e+30
    start_time = -0.5
    end_time = 0.0
    test_surr_times = np.linspace(
        start_time, end_time, sampling_frequency * (end_time - start_time)
    )
    
    test_surr = gwmemory.waveforms.surrogate.Surrogate(
        q=mass_ratio,
        spin_1=[s1x, s1y, s1z],
        spin_2=[s2x, s2y, s2z],
        total_mass=total_mass,
        distance=distance,
        times=test_surr_times,
    )

    new_test_surr_times = test_surr_times / test_surr.t_to_geo
    return test_surr.sur.find_t_0(test_surr.q,test_surr.S1,test_surr.S2,MTot=test_surr.MTot,distance=test_surr.distance,t=new_test_surr_times,LMax=test_surr.LMax)

#time_lim = get_t_0(injection_parameters.items())

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

    GG = 6.674098281543097e-11
    cc = 2.99792458e8
    m_sun_to_kg = 1.98847e+30
    start_time = -0.5
    end_time = 0.0
    #time_lim # from error message
    t_f = time_lim[1]
    start_time = -0.5
    end_time = (t_f+0.9) * (total_mass * m_sun_to_kg)/(cc**3/GG)
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

    plus_new = oscillatory["plus"] + memory_constant * memory["plus"]
    cross_new = oscillatory["cross"] + memory_constant * memory["cross"]

    # Next, we want to place them in our sample space
    plus = np.zeros(len(times))
    cross = np.zeros(len(times))
    plus[-len(surr_times):] = plus_new
    cross[-len(surr_times):] = cross_new

    window_type_plus = ("kaiser", 0.1)
    window_type_cross = ("kaiser", 0.1)

    plus = time_domain_window(plus, window_type=window_type_plus)
    cross = time_domain_window(cross, window_type=window_type_cross)

    return {"plus": plus, "cross": cross}

# iii.
injection_parameters = dict(
    total_mass=60.0,
    s1x=0.0,
    s2x=0.0,
    s1y=0.0,
    s2y=0.0,
    s1z=0.0,
    s2z=0.0,
    distance=200,
    mass_ratio=1.0,
    inc=np.pi / 2,
    psi=0.0,
    phase=0.0,
    memory_constant=1.0,
    ra=0.0,
    dec=0.0,
    geocent_time=0.0
)

time_lim = get_t_0_t_f(mass_ratio=injection_parameters['mass_ratio'],
    s1x=injection_parameters['s1x'],
    s2x=injection_parameters['s2x'],
    s1y=injection_parameters['s1y'],
    s2y=injection_parameters['s2y'],
    s1z=injection_parameters['s1z'],
    s2z=injection_parameters['s2z'],
    distance=injection_parameters['distance'],
    total_mass=injection_parameters['total_mass'],
    phase=injection_parameters['phase'],
    inc=injection_parameters['inc'],
    memory_constant=injection_parameters['memory_constant'])

waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    time_domain_source_model=memory_time_model,
    start_time=injection_parameters["geocent_time"] - 0.5,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters
)

ifos = bilby.gw.detector.InterferometerList(["H1", "L1", "V1"])
ifos.set_strain_data_from_zero_noise(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 0.5,
)
ifos.inject_signal(
    waveform_generator=waveform, parameters=injection_parameters
)

priors = injection_parameters.copy()
#priors["memory_constant"] = bilby.core.prior.Uniform(-1, 3, r"$\lambda$")

likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform
)

likelihood.parameters.update(injection_parameters)


# iv.
memory_constant = np.arange(-1., 3., 0.5)
logL = []
for Lambda in memory_constant:
   likelihood.parameters["memory_constant"] = Lambda
   logL.append(likelihood.log_likelihood())	
   print("Now working on lambda = %s"%Lambda)

def logL_function(Lambda):
   likelihood.parameters["memory_constant"] = Lambda
   return likelihood.log_likelihood()

# v.
max_lambda = opt.fminbound(lambda Lambda: -logL_function(Lambda), -1, 3)
print(max_lambda, logL_function(max_lambda))

#fig = figure(figsize=(6, 6))
plt.plot(memory_constant, logL, color="r")
plt.plot(max_lambda, logL_function(max_lambda), 'bo', label="Max log L")
plt.xlim(-1, 3)
plt.xlabel(r'$lambda$')
plt.ylabel(r'log Likelihood')
plt.legend()
#rc("xtick", labelsize=12)
#rc("ytick", labelsize=12)
#rc("axes", labelsize=14)
plt.grid(True)
plt.savefig('likelihood_v_lambda.pdf')

#tight_layout()
#show()
plt.close()

# vi.
max_lambda = opt.fminbound(lambda Lambda: -logL_function(Lambda), -1, 3)


peak = likelihood.np().argmax()
memory_constant = likelihood.sample_lambda[peak]
two_sigma = np.std(likelihood) * 2

print(r'$\lambda$ = {} $\pm$ {}'.format(memory_constant, two_sigma))

