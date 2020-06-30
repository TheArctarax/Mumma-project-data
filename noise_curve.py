#from __future__ import division, print_function
from matplotlib.pyplot import *
import requests
rcParams["font.family"] = "Times New Roman"


# download L1 strain data from O2 off LIGO DCC
url = 'https://dcc.ligo.org/public/0156/G1801952/001/2017-08-06_DCH_C02_L1_O2_Sensitivity_strain_asd.txt'
r = requests.get(url, allow_redirects=True)
open('l1o2strain.txt', 'wb').write(r.content)

# read the resulting .txt file
f = open('l1o2strain.txt','r')
ldata = f.read()
f.close()

# plot the data
fig = figure(figsize=(6, 6))
plot(*ldata, color='b')
xlim(20, 1400)
ylim(5e-24, 1e-20)
xlabel('frequency (Hz)')
ylabel(r'noise strain [$Hz^{-1/2}$]')
xscale('log')
yscale('log')
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=14)

