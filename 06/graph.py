import core_functions as cf
import matplotlib.pyplot as plt
import numpy as np

x = [0,20,40,60,80,100]
t = [50,100,150,250,500,1000,1250,1500]
[X,T] = np.meshgrid(x,t)

def get_proteinconcentration(x,t):
    return cf.obtain_result(f"./protein_profile/{t}hrs_{x}.txt")


obtainresult=np.vectorize(get_proteinconcentration)
Z = obtainresult(X,T)


plt.figure(figsize=(9.5, 6))
cp = plt.contourf(X, T, Z, cmap = plt.cm.RdBu_r)
cbar = plt.colorbar(cp)

plt.xlabel('distance from soma (μm)')
plt.ylabel('time [hours]')
plt.title('0.05 burst/s protein profile')
plt.show()
