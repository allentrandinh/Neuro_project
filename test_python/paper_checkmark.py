import core_functions as cf
import math

'''
Data for CaMKIIalpha
#v_R : 10**(-3) -> 1.3
#D_P: 4.5 -> 10**(-2)
'''
beta_R = 0.001
beta_P = 0.021
k_R = 1.2*(10**(-5)) #half life 16h
k_P = 1.21*(10**(-6)) #half life 6.64h

#upper_limit_rna
upper_D_R = 3.8*(10**(-3))
upper_v_R = 1.3

#lower limit rna
lower_D_R = 3.8*(10**(-3))
lower_v_R = 5.8*(10**(-3))

'''
#benchmark1: data in graph 1.C
upper0 = cf.steady_mRNA(0,upper_D_R,upper_v_R,k_R,beta_R)
upper50 = cf.steady_mRNA(50,upper_D_R,upper_v_R,k_R,beta_R)
upper100 = cf.steady_mRNA(100,upper_D_R,upper_v_R,k_R,beta_R)
print(upper50/upper0*100,upper100/upper0*100) #99.95 and 99.9

lower0 = cf.steady_mRNA(0,lower_D_R,lower_v_R,k_R,beta_R)
lower50 = cf.steady_mRNA(50,lower_D_R,lower_v_R,k_R,beta_R)
lower100 = cf.steady_mRNA(100,lower_D_R,lower_v_R,k_R,beta_R)
print(lower50/lower0*100, lower100/lower0*100) #90.18 and 81.3
'''

average_D_R = 3.4*(10**(-3))
average_v_R = 6.7*(10**(-3))

S_mRNA = 0.557
z = 0.22
'''
#benchmark2: data in graph 1D, off by 2%
D_P = 4.5
high_D_P_0 = cf.steady_protein(S_mRNA,z,0,D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
high_D_P_50 = cf.steady_protein(S_mRNA,z,50,D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
high_D_P_100 = cf.steady_protein(S_mRNA,z,100,D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
print(high_D_P_50/high_D_P_0*100,high_D_P_100/high_D_P_0*100) #98.33 instead of 98.5

low_D_P = 0.023
low_D_P_0 = cf.steady_protein(S_mRNA,z,0,low_D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
low_D_P_50 = cf.steady_protein(S_mRNA,z,50,low_D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
low_D_P_100 = cf.steady_protein(S_mRNA,z,100,low_D_P,0,k_P,beta_P,lower_D_R,lower_v_R,k_R,beta_R)
print(low_D_P_50/low_D_P_0*100,low_D_P_100/low_D_P_0*100) #0.68.9 instead of 0.71
'''
average_D_P = 0.24
v_P = 0

'''
#benchmark 3: 3A
baseline = cf.steady_mRNA(100,average_D_R,average_v_R,k_R,beta_R)
half_vR = cf.steady_mRNA(100,average_D_R,average_v_R/2,k_R,beta_R)
baseline400 = cf.steady_mRNA(400,average_D_R,average_v_R,k_R,beta_R)
half_vR_400 = cf.steady_mRNA(400,average_D_R,average_v_R/2,k_R,beta_R)
print(half_vR/baseline,half_vR_400/baseline400)
'''

'''
#benchmark 4: 3B
baseline = cf.steady_protein(S_mRNA,z,0,average_D_P,0,k_P,beta_P,average_D_R,average_v_R,k_R,beta_R)
half_vR = cf.steady_protein(S_mRNA,z,0,average_D_P,0,k_P,beta_P,average_D_R,average_v_R/2,k_R,beta_R)
baseline400 = cf.steady_protein(S_mRNA,z,400,average_D_P,0,k_P,beta_P,average_D_R,average_v_R,k_R,beta_R)
half_vR_400 = cf.steady_protein(S_mRNA,z,400,average_D_P,0,k_P,beta_P,average_D_R,average_v_R/2,k_R,beta_R)
print(half_vR/baseline,half_vR_400/baseline)
'''

burst_last = 15 #minutes
peak_translation = 0.035 #protein/s

sigma = 0.5


#print(cf.translational_burst(60,50,0.3,0,40,2,1,1/300,0.1,1/300,0.035,average_D_P,k_P,0))
print(cf.translational_burst(60,50,0.1,0,0.5,1,1/300,0.1,1/300,average_D_P,k_P,0))