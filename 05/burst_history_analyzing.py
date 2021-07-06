import numpy as np
import collections
import matplotlib.pyplot as plt

burst_timestamp_5 = np.load('./data/timestamp_5Hz.npy')
burst_location_5 = np.load('./data/location_5Hz.npy')
burst_timestamp_half = np.load('./data/timestamp_halfHz.npy')
burst_location_half = np.load('./data/location_halfHz.npy')

tmax = 480*3600 #480hrs*3600s

def count_burst_everysecond(burst_timestamp_array,tmax):
    current_index = 0
    count_each_sec =[]
    for each_sec in range(0,tmax):
        burst_num = 0
        for index in range(current_index,len(burst_timestamp_array)):
            if burst_timestamp_array[index] > each_sec+1:
                break
            elif each_sec <= burst_timestamp_array[index] <= each_sec+1:
                current_index += 1
                burst_num += 1
        count_each_sec.append(burst_num)
    # keys are array elements and values are their corresponding frequencies
    return(collections.Counter(count_each_sec))

'''
#testing
print(burst_timestamp_half[0:20])
print(count_burst_everysecond(burst_timestamp_half,13))'''

'''
#plotting of burst frequency
# keys are array elements and values are their corresponding frequencies
halfHz_freq = count_burst_everysecond(burst_timestamp_half,tmax)
fiveHz_freq = count_burst_everysecond(burst_timestamp_5,tmax)
plt.bar(fiveHz_freq.keys(),fiveHz_freq.values())
plt.xlabel('number of bursts per second')
plt.ylabel('total occurrence')
plt.title('5Hz bursts frequencies')
plt.show()'''

'''
#plotting of location frequencies

location_freq_half = collections.Counter(burst_location_half)
plt.bar(location_freq_half.keys(),location_freq_half.values())
plt.show()'''

location_freq_five = collections.Counter(burst_location_5)
plt.bar(location_freq_five.keys(),location_freq_five.values())
plt.xlabel('distance from soma')
plt.ylabel('number of bursts')
plt.title('5Hz burst location frequency')
plt.show()