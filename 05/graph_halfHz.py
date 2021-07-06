import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def obtain_result(location, timestamp):
    file_name = "./protein_profile_0.5Hz/halfHz_"+str(timestamp)+"hrs_"+str(location)+".txt"
    full_data = pd.read_csv(file_name,sep="\t")
    protein_level = full_data.iloc[-1]['protein']
    return protein_level

x = [0,20,40,60,80]
t = [48, 144, 240, 336, 432]
[X,T] = np.meshgrid(x,t)

obtainresult=np.vectorize(obtain_result)
Z = obtainresult(X,T)

plt.figure(figsize=(9.5, 6))
cp = plt.contourf(X, T, Z, 10, cmap = plt.cm.RdBu_r)
#the 9000 is to normalize the function. The 55 ensures the time axis is accurate
cbar = plt.colorbar(cp)
#clabel(cp, inline=True, fontsize=8)
plt.xlabel('distance from soma (μm)')
plt.ylabel('time [hours]')
plt.title('0.5 burst/s protein profile')
plt.show()