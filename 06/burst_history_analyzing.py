import numpy as np
import collections
import matplotlib.pyplot as plt

burst_timestamp = np.load('./data/timestamp_0.05Hz.npy')
burst_location = np.load('./data/location_0.05Hz.npy')

tmax = 1500*3600 #480hrs*3600s

def count_burst_everymin(burst_timestamp_array,tmax):
    current_index = 0
    count_each_sec =[]
    for each_sec in range(0,tmax,60):
        burst_num = 0
        for index in range(current_index,len(burst_timestamp_array)):
            if burst_timestamp_array[index] > each_sec+60:
                break
            elif each_sec <= burst_timestamp_array[index] <= each_sec+60:
                current_index += 1
                burst_num += 1
        count_each_sec.append(burst_num)
    # keys are array elements and values are their corresponding frequencies
    return(collections.Counter(count_each_sec))


'''
#plotting of burst frequency
# keys are array elements and values are their corresponding frequencies
burst_freq = count_burst_everymin(burst_timestamp,tmax)

plt.bar(burst_freq.keys(),burst_freq.values())
plt.xlabel('number of bursts per min')
plt.ylabel('total occurrence')
plt.title('0.05Hz bursts frequencies')
plt.show()


#plotting of location frequencies

location_freq = collections.Counter(burst_location)
plt.bar(location_freq.keys(),location_freq.values())
plt.show()'''
