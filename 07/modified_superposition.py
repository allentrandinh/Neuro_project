import core_functions as cf
import numpy as np
import json
import math

'''Setting'''
burst_frequency = 5 # number of burst per second
window_time = 200*3600 #in seconds
dendritic_length = 100
time_gap = 600
half_life = 6.6 #days

'''Code'''
burst_history = json.load(open('../burst_history/200_2Hz_history'))

with open('./data/x50_proteinprofile.txt','a') as f:
    f.write('t\tprotein\n')
    for i in range(0,window_time,time_gap):
        print(f"Currently working on time {i}")
        f.write(f"{i}\t{cf.modified_superposition(burst_history,i,50,kP=math.log(2,math.e)/(half_life*24*60*60))}\n")
