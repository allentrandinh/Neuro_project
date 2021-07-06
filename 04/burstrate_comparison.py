import concurrent.futures
import core_functions as cf
import math

'''Variables to change'''
max_time = 600 #hours
rate_list = [0.01,0.1,1.0,10]
burst_duration = 10
'''End of variables to change'''

dP = 0.243798
vP = 0
protein_halflife_days = 6.6
LifTP = protein_halflife_days*24*60*60
kP=math.log(2,math.e)/(LifTP)
sigma = 0.50
a1 = 1
a2=1/(burst_duration/2*60)
b1=0.1
b2=1/(burst_duration/2*60)

tmax = int(max_time*3600) #in seconds

def x50_temporal_profile(burst_rate):
    time_to_100_bursts = int(100/burst_rate)
    filename = "x50_"+str(burst_rate)+".txt"
    current_protein=0
    with open(filename, 'w') as f:
        f.write(f"t\tprotein\n")
        for t in range(0, tmax, time_to_100_bursts):
            print(f"Functioned being called with {t + time_to_100_bursts}")
            for position in range(0, 100):
                # calculation a step = gap ahead
                current_protein += cf.translational_burst(50, position, t + time_to_100_bursts, 0, sigma, a1, a2, b1, b2, dP, kP, vP)
            f.write(f"{t + time_to_100_bursts}\t{current_protein}\n")

if __name__=='__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for rate in rate_list:
            executor.submit(x50_temporal_profile,rate)


