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

#max_time_display = 100
#xrange = 50
#x0 = 300

sigma = 0.50

t0 = 0
a1 = 1
a2=1/(burst_duration/2*60)
b1=0.1
b2=1/(burst_duration/2*60)


def random_pos(burst_per_sec, total_length):
    list = []
    while len(list) < burst_per_sec:
        temp = np.random.randint(0, total_length)
        if temp not in list:
            list.append(temp)
    return list

def burst_history_generator(burst_per_sec,total_length,time_length,time_gap):

    #generate a history of burst over time: structure of dictionary
    #key = time point; values = corresponding positions of burst
    history = {}
    for t in range(0,time_length,time_gap):
        history[t] = random_pos(burst_per_sec,total_length)
    return history

def poisson_burst_history(ave_burst_per_sec,total_length,time_length):
    '''
    :param ave_burst_per_sec:
    :param total_length:
    :param time_length: in seconds
    :return:
    '''
    burst_count = np.random.poisson(ave_burst_per_sec,time_length)
    history = {}
    for t in range(0, time_length):
        history[t] = random_pos(burst_count[t],total_length)
    return(history)

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

vtf=np.vectorize(cf.tf)
timearray= np.double(np.linspace(0, max(1/a2,1/a1)*50, 100))
Tmaxi=max(vtf(timearray,t0,a1,a2,b1,b2))

'''
Settings
'''
ave_burst_per_sec = 5 #bursts per second
length = 100 #micrometer
max_time_display = 15 #mins
tmax = int(max_time_display*60) #in second
#time_gap = 60 #how spaced out are bursts
number_data_point = 10

'''
End of settings
'''
#structure of burst_history dictionary
# time_of_burst: list(location_of_burst)
#burst_history = burst_history_generator(burst_per_sec,length,tmax,time_gap)
burst_history = poisson_burst_history(ave_burst_per_sec,length,tmax)

print(burst_history)

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
mp.bar(x,y)

#generate all points and time of interest
x = np.double(np.linspace(0,length,number_data_point))
t = np.double(np.linspace(0, tmax, number_data_point))
[X,T]=np.meshgrid(x,t)

#vectorize superposition function
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
