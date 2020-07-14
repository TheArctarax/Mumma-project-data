from __future__ import division, print_function
import matplotlib
import numpy as np
import gwmemory
from gwmemory import utils as utils
import itertools
from scipy.signal import get_window

np.seterr(divide="ignore", invalid="ignore")


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
sample_lambda = np.linspace(a, b, sampling_frequency * (b-a))

# ii.

# iii.

# iv.

# v.
fig = figure(figsize=(6, 6))
plot(sample_lambda, likelihood, color="r")

xlim(-1, 2)
xlabel(r'$\lambda$')
ylabel(r'Likelihood($\lambda$)')
rc("xtick", labelsize=12)
rc("ytick", labelsize=12)
rc("axes", labelsize=14)
grid(False)

savefig('likelihood_v_lambda')

tight_layout()
show()
close()

# vi.
peak = likelihood.np().argmax()
memory_constant = likelihood.sample_lambda[peak]
two_sigma = np.std(likelihood) * 2

print(r'$\lambda$ = {} $\pm$ {}'.format(memory_constant, two_sigma))

