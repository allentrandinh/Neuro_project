import core_functions as cf
import os

os.mkdir('./data')
cf.poisson_interval_history_generator(burst_rate=0.05,total_length=100,hours_total=4000,path='./data')
