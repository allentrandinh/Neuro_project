import core_functions as cf
import math
import numpy as np

'''Finding interval whose protein contribution greater than threshold of all position'''

sigma=0.50
a1=1
a2=1/300
b1=0.1
b2=1/300
dP=0.243798
kP=math.log(2,math.e)/(6.6*24*60*60)
vP=0


'''
check_middle checks if we have any transition in the middle, not currently being used
recursive_index_finder finds the index around the transition with specified setting, works if we are sure they are in range t_low to t_high


'''

def check_middle(interested_position,t_low,t_high, threshold,setting,number_of_points=20):
    '''
    :param interested_position: x position along the dendrite being considered
    :param t_high: t upper bound
    :param t_low: t lower bound
    :param threshold: threshold
    :param setting: "low" or "high". For example if protein contribution at both t_low and t_high are smaller than threshold, setting will be low
    :return: whether or not we have transition along time. True if have, False otherwise
    '''
    check_indices = np.linspace(t_low,t_high,number_of_points)

    if setting == "low":
        for index in check_indices:
            if cf.translational_burst(interested_position,0,index,0,sigma,a1,a2,b1,b2,dP,kP,vP) > threshold:
                return (True,index)

    if setting == "high":
        for index in check_indices:
            if cf.translational_burst(interested_position,0,index,0,sigma,a1,a2,b1,b2,dP,kP,vP) < threshold:
                return (True,index)

    return (False,None)

def recursive_index_finder(interested_position,t_low,t_high, threshold,setting):

    if setting=="low2high":
        #if protein at mid > threshold -> transition in first half range and vice versa
        if t_high == t_low+1:
            if t_high >= threshold:
                return (t_low,t_high)
            else:
                return (None,None)

        mid_t = int(math.ceil((t_low + t_high) / 2))
        protein_at_mid = cf.translational_burst(interested_position, 0, mid_t, 0, sigma, a1, a2, b1, b2, dP, kP, vP)

        if protein_at_mid > threshold:
            return(recursive_index_finder(interested_position,t_low,mid_t,threshold,"low2high"))
        elif protein_at_mid == threshold:
            return(mid_t,mid_t+1)
        else:
            return(recursive_index_finder(interested_position,mid_t,t_high,threshold,"low2high"))

    if setting == "high2low":
        # if protein at mid < threshold -> transition in first half range and vice versa
        if t_high == t_low+1:
            if t_high <= threshold:
                return (t_low,t_high)
            else:
                return (None,None)

        mid_t = int(math.ceil((t_low + t_high) / 2))
        protein_at_mid = cf.translational_burst(interested_position, 0, mid_t, 0, sigma, a1, a2, b1, b2, dP, kP, vP)

        if protein_at_mid < threshold:
            return (recursive_index_finder(interested_position, t_low, mid_t, threshold, "low2high"))
        elif protein_at_mid == threshold:
            return(mid_t,mid_t-1)
        else:
            return (recursive_index_finder(interested_position, mid_t, t_high, threshold, "low2high"))

def first_higher(interested_position,t_low,t_high,threshold,gap):
    for i in range(t_low,t_high,gap):
        if cf.translational_burst(interested_position,0,i,0,sigma, a1, a2, b1, b2, dP, kP, vP) >= threshold:
            return i
    return None

def first_lower(interested_position,t_low,t_high,threshold,gap):
    for i in range(t_low,t_high,gap):
        if cf.translational_burst(interested_position, 0, i, 0, sigma, a1, a2, b1, b2, dP, kP, vP) <= threshold:
            return i
    return None

def get_result(index1,index2,var,threshold, interested_position,sigma, a1, a2, b1, b2, dP, kP, vP):
    if cf.translational_burst(interested_position, 0, index1, 0, sigma, a1, a2, b1, b2, dP, kP, vP) > threshold:
        var = index1
    elif cf.translational_burst(interested_position, 0, index2, 0, sigma, a1, a2, b1, b2, dP, kP, vP) > threshold:
        var = index2

def comprehensive_finder(interested_position,t_low,t_high,threshold,gap):
    low2high = 0
    high2low = 0
    pro_low = cf.translational_burst(interested_position,0,t_low,0,sigma, a1, a2, b1, b2, dP, kP, vP)

    #starts out higher, only go down
    if pro_low > threshold:
        temp1, temp2 = recursive_index_finder(interested_position,t_low,t_high,threshold,setting="high2low")
        get_result(temp1,temp2,high2low,threshold,interested_position,sigma, a1, a2, b1, b2, dP, kP, vP)

    else:
        first_high = first_higher(interested_position,t_low,t_high,threshold,gap)
        if first_high:
            temp1, temp2 = recursive_index_finder(interested_position,t_low,first_high,threshold,setting="low2high")
            get_result(temp1,temp2,low2high,threshold,interested_position,sigma, a1, a2, b1, b2, dP, kP, vP)
            first_low = first_lower(interested_position,first_high,t_high,threshold,gap)
            if first_low:
                temp3,temp4 = recursive_index_finder(interested_position,first_high,first_low,threshold,setting="high2low")
                get_result(temp3,temp4,high2low,threshold,interested_position,sigma, a1, a2, b1, b2, dP, kP, vP)

    return(low2high,high2low)

