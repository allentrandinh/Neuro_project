import numpy as np

#burst rate = 5 part

burst_rate = 5
total_length = 100
hours_total = 480 #=20 days
tmax = hours_total*3600 #1hrs = 3600s

#average interval between two bursts
ave_interval = 1/burst_rate
#number of interval within the interested duration
interval_num = int(tmax/ave_interval)

#add 0 to get the loop started
burst_timestamp = [0]
burst_location = []


#divide into chunk of 10,000 data points, times 10 then divide by 10 to avoid integer
chunk_num = int(interval_num/10000)
print(chunk_num)
for j in range(0,chunk_num):
    print(f"Call with {j}")
    interval = np.around(np.random.poisson(ave_interval*10,10000)/10,decimals=1)
    position = np.random.randint(0,total_length,10000)
    burst_location.extend(position)
    for i in range(0,10000):
        burst_timestamp.append(np.around(burst_timestamp[-1]+interval[i],decimals=1))

burst_timestamp.remove(0)

np.save('./data/timestamp_5Hz',burst_timestamp)
np.save('./data/location_5Hz',burst_location)
# np.random.randint(0,total_length,burst_per_sec)
# np.random.poisson(ave_burst_per_sec,time_length)''''''

'''
burst_rate = 0.5
total_length = 100
hours_total = 480 #=20 days
tmax = hours_total*3600 #1hrs = 3600s

#average interval between two bursts
ave_interval = 1/burst_rate
#number of interval within the interested duration
interval_num = int(tmax/ave_interval)

#add 0 to get the loop started
burst_timestamp = [0]
burst_location = []


#divide into chunk of 10,000 data points, times 10 then divide by 10 to avoid integer
chunk_num = int(interval_num/10000)
print(chunk_num)
for j in range(0,chunk_num):
    print(f"Call with {j}")
    interval = np.around(np.random.poisson(ave_interval*10,10000)/10,decimals=1)
    print(interval)
    position = np.random.randint(0,total_length,10000)
    burst_location.extend(position)
    for i in range(0,10000):
        burst_timestamp.append(np.around(burst_timestamp[-1]+interval[i],decimals=1))

burst_timestamp.remove(0)

np.save('./data/timestamp_halfHz',burst_timestamp)
np.save('./data/location_halfHz',burst_location)'''
# np.random.randint(0,total_length,burst_per_sec)