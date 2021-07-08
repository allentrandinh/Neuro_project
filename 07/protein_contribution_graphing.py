import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('./data/protein_contribution.txt',sep="\t")

axis = data['t'][:40*60]
x_0 = data['x=0'][:40*60]
x_50 = data['x=50'][:40*60]
x_100 = data['x=100'][:40*60]

plt.plot(axis/60,x_0,label="x @ 0")
plt.plot(axis/60,x_50,label="x @ 50")
plt.plot(axis/60,x_100,label="x @ 100")
plt.ylabel("Protein contribution")
plt.xlabel("Time (min)")
plt.legend()
plt.show()

'''Chosen threshold: 0.00010185865710213 = 1*10^-4 (t=787, x=100)'''

