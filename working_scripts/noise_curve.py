#HereIamimportingpandas,anoften-useddatahandlingpython
#package.
#importssl
#context=ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
import gwpy
import matplotlib.pyplot as plt
import pandas as pd
import wget
import itertools
from matplotlib.pyplot import rcParams
rcParams["font.family"]="Times New Roman"
rcParams['axes.unicode_minus'] = False


# download L1 strain data from O2 off LIGO DCC.This is only required to be run once. Otherwise, duplicate files will be downloaded.
#
#urls=('https://dcc.ligo.org/public/0156/G1801952/001/2017-08-06_DCH_C02_L1_O2_Sensitivity_strain_asd.txt',
#     'https://dcc.ligo.org/public/0156/G1801950/001/2017-06-10_DCH_C02_H1_O2_Sensitivity_strain_asd.txt',
#     'https://dcc.ligo.org/public/0157/P1800374/001/Hrec_hoft_V1O2Repro2A_16384Hz.txt')
#names=('./L1_O2_Sensitivity_strain_asd.txt', './H1_O2_Sensitivity_strain_asd.txt', './V1_O2_Sensitivity_strain_asd.txt')
#
#for files, fnames in itertools.izip(urls, names):
#   wget.download(files, fnames)


# read the resulting .txt file using pandas. I renamed the columns as frequencies and asd.
dfl=pd.read_csv('L1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfh=pd.read_csv('H1_O2_Sensitivity_strain_asd.txt',sep="\t",index_col=False,header=None,names=["frequency","asd"])
dfv=pd.read_csv('V1_O2_Sensitivity_strain_asd.txt',sep="\s+",index_col=False,header=None, names=["frequency","asd"])


# plot the data
fig = plt.figure(figsize=(4.5, 4.5))
plt.loglog(dfl['frequency'], dfl['asd'], color='gwpy:ligo-livingston', label='Livingston')
plt.loglog(dfh['frequency'], dfh['asd'], color='gwpy:ligo-hanford', label='Hanford')
plt.loglog(dfv['frequency'], dfv['asd'], color='gwpy:virgo', label='Virgo')
plt.xlim(20, 5000)
#plt.ylim(5e-24,1e-20)
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'Noise Strain [Hz$^{-1/2}$]')
plt.grid(False)
plt.legend(loc='upper right', prop={'size':12})
plt.rc('xtick',labelsize=10)
plt.rc('ytick',labelsize=10)
plt.rc('axes',labelsize=12)

plt.savefig("noise_curve.pdf")

plt.tight_layout()
plt.show()
plt.close()

