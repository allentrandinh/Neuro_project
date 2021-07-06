import numpy as np
import core_functions as cf
import math
import concurrent.futures

burst_timestamp_5 = np.load('./data/timestamp_5Hz.npy')
burst_location_5 = np.load('./data/location_5Hz.npy')
# burst_timestamp_half = np.load('./data/timestamp_halfHz.npy')
# burst_location_half = np.load('./data/location_halfHz.npy')

#interested time point: 48 hrs, 144, 240, 336, 432
#interested location 0,20,40,60,80,100

def accumulative_protein(time_arr, location_arr,interested_t, interested_x,curr_index, curr_protein,filename):
    burst_duration = 10
    dP = 0.243798
    vP = 0
    protein_halflife_days = 6.6
    LifTP = protein_halflife_days * 24 * 60 * 60
    kP = math.log(2, math.e) / (LifTP)
    sigma = 0.50
    a1 = 1
    a2 = 1 / (burst_duration / 2 * 60)
    b1 = 0.1
    b2 = 1 / (burst_duration / 2 * 60)
    protein = curr_protein
    index = curr_index
    max_t = interested_t*3600 #inseconds
    with open(filename, 'a') as f:
        f.write(f"t\tprotein\n")
        for i in range(index, len(time_arr)):
            print(f"Called with index {i}")
            if time_arr[i] > max_t:
                break
            else:
                protein += cf.translational_burst(interested_x, location_arr[i], max_t, time_arr[i], sigma, a1, a2, b1, b2, dP, kP, vP)
                f.write(f"{i}\t{protein}\n")

def helper(arg):
    accumulative_protein(arg[0],arg[1],arg[2],arg[3],arg[4],arg[5],arg[6])

#interested time point: 48 hrs, 144, 240, 336, 432

if __name__=='__main__':
    arg_names = [
        (burst_timestamp_5, burst_location_5, 240, 0, 1485091, 43253044.4190018,"5Hz_240hrs_0.txt"),
        (burst_timestamp_5, burst_location_5, 240, 20, 1491711, 43610181.9305894, "5Hz_240hrs_20.txt"),
        (burst_timestamp_5, burst_location_5, 240, 40, 1538257, 45482492.4019447, "5Hz_240hrs_40.txt"),
        (burst_timestamp_5, burst_location_5, 240, 60, 1513896, 44523072.5330487, "5Hz_240hrs_60.txt"),
        (burst_timestamp_5, burst_location_5, 240, 80, 1508368, 44252928.5756047, "5Hz_240hrs_80.txt"),
        (burst_timestamp_5, burst_location_5, 240, 100, 1505951, 44052414.8354982, "5Hz_240hrs_100.txt")
    ]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0,len(arg_names)):
            executor.submit(helper,arg_names[i])

# accumulative_protein(burst_timestamp_half,burst_location_half,48,0,0,0,"halfHz_48hrs_0.txt")

# print(burst_timestamp_half[0:50])
# print(burst_location_half[0:50])
