import core_functions as cf
import math
import numpy as np
import matplotlib.pyplot as mp

#default setting as appear in the html translational file
dP = 0.243798
vP = 0
protein_halflife_days = 6.6
LifTP = protein_halflife_days*24*60*60
kP=math.log(2,math.e)/(LifTP)
burst_duration = 10
sigma=0.5
t0 = 0
a1 = 1
a2=1/(burst_duration/2*60)
b1=0.1
b2=1/(burst_duration/2*60)

'''
vtf=np.vectorize(cf.tf)
timearray= np.double(np.linspace(0, max(1/a2,1/a1)*50, 100))
Tmaxi=max(vtf(timearray,t0,a1,a2,b1,b2))
'''
Tmaxi = 99

'''
Settings
'''
ave_burst_per_sec = 0.5 #bursts per second
ave_burst_per_min = ave_burst_per_sec*60
length = 100 #micrometer
max_time_display = 10*24 #hours = 10  days
tmax = int(max_time_display*60*60) #in second
#time_gap = 60 #how spaced out are bursts
number_data_point = 100

burst_history = cf.poisson_burst_history(ave_burst_per_sec,length,tmax)

'''plotting distribution frequency'''
'''
x_axis = []
y_axis = []

distribution = {}
for location in range(0,length+1):
    distribution[location] = 0

for burst_time in burst_history.keys():
    for burst_location in burst_history[burst_time]:
        distribution[burst_location] += 1

x = distribution.keys()
y = distribution.values()

frequency = mp.figure(num=1)
mp.bar(x,y)'''
'''end of plotting frequency'''

'''plotting protein profile'''
#generate all points and time of interest
x = 50
time_steady=[]
highest=0
with open('x_50.txt',"w") as f:
    f.write("Time_stamp\t Current_superposition\t Recent_increment\n")
    for t in range(0,tmax+1,1800):
        print(f"Current call with {t}")
        temp = cf.superposition(burst_history,t,x)
        f.write(f"{t}\t{temp}\t{temp-highest}\n")
        if temp>highest:
            highest=temp
        else:
            if len(time_steady)<=5:
                time_steady.append(t)
with open('x_50_timesteady.txt','w') as f:
    f.write(f"With x=50, 0.5 burst/s the time needs to reach steady state is {time_steady}")


'''
[X,T]=np.meshgrid(x,t)

#vectorize superposition function
superpos = np.vectorize(cf.superposition)

Z=(2.1/60)/99*superpos(burst_history,T,X)

#rendering image
plot_2= mp.figure(figsize=(9.5, 6),num=2)
cp = mp.contourf(X, T/60, Z, 10, cmap = 'RdBu_r')
#the 9000 is to normalize the function. The 55 ensures the time axis is accurate
cbar = mp.colorbar(cp)
#clabel(cp, inline=True, fontsize=8)
mp.xlabel('distance from soma (μm)')
mp.ylabel('time [minutes]')
#end of plotting protein profile

mp.show()
'''