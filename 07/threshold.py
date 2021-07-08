import math
import core_functions as cf

'''Setting'''
sigma=0.50
a1=1
a2=1/300
b1=0.1
b2=1/300
dP=0.243798
kP=math.log(2,math.e)/(6.6*24*60*60)
vP=0

'''
Generating a series of calculation on protein contribution over time for one burst
Max_distance from burst site = 100 um
initial attempt at x=0,50,and 100


hours = 10

with open('./data/protein_contribution.txt',"w") as f:
    f.write("t\tx=0\tx=50\tx=100\n")
    for i in range(0,hours*3600):
        f.write(f"{i}\t{cf.translational_burst(0,0,i,0,sigma,a1,a2,b1,b2,dP,kP,vP)}\t"
                f"{cf.translational_burst(50,0,i,0,sigma,a1,a2,b1,b2,dP,kP,vP)}\t"
                f"{cf.translational_burst(100,0,i,0,sigma,a1,a2,b1,b2,dP,kP,vP)}\n")

'''

hours = 100

with open('./data/full_proteincontribution.txt','w') as f:

    #adding title
    title = "t"
    for i in range(0,101):
        title += f"\tx={i}"
    title += "\n"
    f.write(title)

    #adding protein concentration along time
    for i in range(0,hours*3600):
        line = f"{i}"
        print(f"Working on t={i}")
        #adding protein concentration at each position
        for k in range(0,101):
            line += f"\t{cf.translational_burst(k,0,i,0,sigma,a1,a2,b1,b2,dP,kP,vP)}"
        line += "\n"
        f.write(line)

'''taking too slow'''