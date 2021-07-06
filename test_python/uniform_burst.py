import core_functions as cf
import math
import concurrent.futures
'''
0.5 burst/sec -> 30 bursts/minute -> round to 100 bursts/3mins
5 bursts/sec -> 100 bursts/20s
'''
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

def burst_continue(file_name,gap,start_time, start_protein,tmax):
    file = file_name+'.txt'
    with open(file, 'w') as f:
        f.write(f"t\tprotein\n")
        for t in range(start_time, tmax, gap):
            print(f"Functioned being called with {t + gap}")
            for position in range(0, 100):
                # calculation a step = gap ahead
                start_protein += cf.translational_burst(50, position, t + gap, 0, sigma, a1, a2, b1, b2, dP, kP, vP)
            f.write(f"{t + gap}\t{start_protein}\n")

def helper(args):
    burst_continue(args[0], args[1], args[2], args[3], args[4])

if __name__=='__main__':
    arg_names = [('x50_5hz_2',20,2001660,428352594.683651,2160000),
                 ('x50_10hz_2',10,1010950,773168846.585799,2160000),
                 ('x50_4hz_2',25,1257600,322912844.297582,2160000)]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0,len(arg_names)):
            executor.submit(helper,arg_names[i])
