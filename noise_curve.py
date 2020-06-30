# Here I am importing pandas, an often-used data handling python
# package.
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import requests
import pandas as pd
import wget
#rcParams["font.family"] = "Times New Roman"


# download L1 strain data from O2 off LIGO DCC
url = 'https://dcc.ligo.org/public/0156/G1801952/001/2017-08-06_DCH_C02_L1_O2_Sensitivity_strain_asd.txt'
wget.download(url, './L1_O2_Sensitivity_strain_asd.txt')

#r = requests.get(url, allow_redirects=True)
#open('l1o2strain.txt', 'wb').write(r.content)

# read the resulting .txt file
df = pd.read_csv('L1_O2_Sensitivity_strain_asd.txt', sep=" ", header=None, names=["frequency", "asd"])

#f = open('l1o2strain.txt','r')
#ldata = f.read()
#f.close()

# plot the data
#fig = figure(figsize=(6, 6))
plt.plot(df['frequency'], df['asd'], color='b')
plt.xlim(20, 1400)
plt.ylim(5e-24, 1e-20)
plt.xlabel('frequency (Hz)')
plt.ylabel(r'noise strain [$Hz^{-1/2}$]')
plt.xscale('log')
plt.yscale('log')
#plt.rc('xtick', labelsize=12)
#plt.rc('ytick', labelsize=12)
#plt.rc('axes', labelsize=14)
plt.savefig("noise_curve.pdf")
