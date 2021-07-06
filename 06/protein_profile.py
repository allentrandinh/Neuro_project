import numpy as np
import core_functions as cf
import concurrent.futures

burst_timestamp = np.load('./data/timestamp_0.05Hz.npy')
burst_location = np.load('./data/location_0.05Hz.npy')

'''
need to manually enter in intersted time and location
total_hours of history is 4000 hrs
50, 100, 150,// 250,500,750,1000,1250,1500,
'''

def helper(arg):
    cf.accumulative_protein(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6],protein_halflife_days=3.3)

if __name__ == '__main__':
    arg_names = [
        (burst_timestamp, burst_location, 1000, 0, 0, 0, "./protein_profile/1000hrs_0.txt"),
        (burst_timestamp, burst_location, 1000, 20, 0, 0, "./protein_profile/1000hrs_20.txt"),
        (burst_timestamp, burst_location, 1000, 40, 0, 0, "./protein_profile/1000hrs_40.txt"),
        (burst_timestamp, burst_location, 1000, 60, 0, 0, "./protein_profile/1000hrs_60.txt"),
        (burst_timestamp, burst_location, 1000, 80, 0, 0, "./protein_profile/1000hrs_80.txt"),
        (burst_timestamp, burst_location, 1000, 100, 0, 0, "./protein_profile/1000hrs_100.txt")
    ]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0, len(arg_names)):
            executor.submit(helper, arg_names[i])