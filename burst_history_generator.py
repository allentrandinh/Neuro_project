import core_functions as cf
import json

'''Setting'''
burst_frequency = 2 # number of burst per second
window_time = 200 #in hours
dendritic_length = 100

'''Code'''
burst_history = cf.poisson_burst_history_nonrepeat(burst_frequency,dendritic_length,window_time*3600)
json.dump(burst_history, open(f"./burst_history/{window_time}_{burst_frequency}Hz_history",'w'))