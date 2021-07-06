import pandas as pd
import matplotlib.pyplot as plt

# rate_1 = pd.read_csv('./x50_0.01.txt',sep='\t')
rate_2 = pd.read_csv('./x50_1Hz.txt',sep='\t')
rate_3 = pd.read_csv('./x50_4Hz.txt',sep='\t')
rate_4 = pd.read_csv('./x50_5Hz.txt',sep='\t')
#rate_5 = pd.read_csv('./x50_10Hz_filtered.txt',sep='\t')

# rate_1['hours'] = rate_1['t']/3600
rate_2['hours'] = rate_2['t']/3600
rate_3['hours'] = rate_3['t']/3600
rate_4['hours'] = rate_4['t']/3600
#rate_5['hours'] = rate_5['t']/3600

# plt.plot(rate_1['hours'],rate_1['protein'],label='0.01 burst/s')
plt.plot(rate_2['hours'],rate_2['protein'],label='1 burst/s')
plt.plot(rate_3['hours'],rate_3['protein'],label='4 bursts/s')
plt.plot(rate_4['hours'],rate_4['protein'],label='5 bursts/s')
#plt.plot(rate_5['hours'],rate_5['protein'],label='10 bursts/s')

plt.legend()
plt.show()