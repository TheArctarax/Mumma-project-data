import pandas as pd
import matplotlib.pyplot as plt

dfv=pd.read_csv('V1_O2_Sensitivity_strain_asd.txt',sep=" ",index_col=False,header=None, names=["frequency","asd"])

plt.plot(dfv['frequency'], dfv['asd'])
plt.xscale('log')
plt.yscale('log')

plt.savefig('v1_test.pdf')

plt.show()
plt.close()
