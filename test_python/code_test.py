import core_functions as cf
import math
import random
import numpy as np
import matplotlib.pyplot as mp
import pandas as pd


def superposition(burst_history,t_of_interest,x_of_interest,sigma=0.50,a1=1,a2=1/300,b1=0.1,b2=1/300,dP=0.243798,kP=math.log(2,math.e)/(6.6*24*60*60),vP=0):
    #find the total protein concentration at location x_of_interest, time t_of_interest by summing up all component bursts
    print(f"Function being called with {t_of_interest} and {x_of_interest}")
    total = 0
    # sum up all burst in a given time
    for burst_time in burst_history.keys():
        #only sum bursts events that happen strictly before time of interest
        if burst_time < t_of_interest:
            #sum up bursts in all positions
            for burst_position in burst_history[burst_time]:
                total += cf.translational_burst(x_of_interest,burst_position,t_of_interest,burst_time,sigma,a1,a2,b1,b2,dP,kP,vP)
    return total

'''
#default setting as appear in the html translational file
dP = 0.243798
vP = 0
protein_halflife_days = 6.6
LifTP = protein_halflife_days*24*60*60
kP=math.log(2,math.e)/(LifTP)
burst_duration = 10

#max_time_display = 100
#xrange = 50
#x0 = 300

sigma = 0.50

t0 = 0
a1 = 1
a2=1/(burst_duration/2*60)
b1=0.1
b2=1/(burst_duration/2*60)

vtf=np.vectorize(cf.tf)
timearray= np.double(np.linspace(0, max(1/a2,1/a1)*50, 100))
Tmaxi=max(vtf(timearray,t0,a1,a2,b1,b2))

#burst_history = {0:[300,200]} #test_!
#burst_history = {0:[300,200],60:[300,200]} #test 2
burst_history = {0:[330,300,200,250]} #test_3
x = [50,100,200,250,300,350,400,500]
t = [time for time in range(0,40*60,60)]
[X,T]=np.meshgrid(x,t)

superpos = np.vectorize(superposition)

Z=(2.1/60)/Tmaxi*superpos(burst_history,T,X)

#rendering image
plot_2= mp.figure(figsize=(9.5, 6),num=2)
cp = mp.contourf(X, T/60, Z, 10, cmap = 'RdBu_r')
#the 9000 is to normalize the function. The 55 ensures the time axis is accurate
cbar = mp.colorbar(cp)
#clabel(cp, inline=True, fontsize=8)
mp.xlabel('distance from soma (μm)')
mp.ylabel('time [minutes]')
mp.show()
'''

df = pd.read_csv("./cleaned_x50_50Hz.txt",sep="\t")

df['hours']=df['t']/3600
mp.scatter(df['hours'],df['protein'])
mp.xlabel('time [hours]')
mp.ylabel('protein accumulation')
mp.show() 

