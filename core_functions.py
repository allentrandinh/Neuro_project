from mpmath import *
import numpy as np
import math
import pandas as pd
'''
Contains core functions in mmc2 paper
'''

def steady_mRNA(x,D_R, v_R, k_R, beta_R):
    '''
    :param x: dendritic distance
    :param D_R: diffusion coefficient for mRNA
    :param v_R: transport velocity of mRNA
    :param k_R: degradation rate of mRNA
    :param beta_R: transcription rate of mRNA
    :return: steady mRNA level
    '''
    lamda_R = (sqrt((v_R**2)+4*D_R*k_R)-v_R)/(2*D_R)
    return (beta_R*lamda_R/k_R)*exp(-(lamda_R*x))

def steady_dendritic_protein(x,D_P, v_P, k_P, beta_P,D_R, v_R, k_R, beta_R):
    '''
    :param x: dendritic distance
    :param D_P: diffusion coefficient for protein
    :param v_P: transport velocity for protein
    :param k_P: degradation rate for protein
    :param beta_P: translation rate
    :param D_R: diffusion coefficient for mRNA
    :param v_R: transport velocity for mRNA
    :param k_R: degradation rate for mRNA
    :param beta_R: transcription rate of mRNA
    :return: steady dendritic protein level
    '''
    lamda_R = (sqrt((v_R ** 2) + 4 * D_R * k_R) - v_R) / (2 * D_R)
    lamda_P = (sqrt((v_P**2)+4*D_P*k_P)-v_P)/(2*D_P)
    first_term = (beta_P*beta_R*lamda_R)/(k_R*(D_P*(lamda_R**2)+(v_P*lamda_R)-k_P))
    second_term = (-(exp(-(lamda_R*x))))+((D_P*lamda_P*lamda_R+v_P*lamda_P)/k_P)*exp(-(lamda_P*x))
    return first_term*second_term

def steady_protein(S_mRNA, z, x,D_P, v_P, k_P, beta_P,D_R, v_R, k_R, beta_R):
    '''
    :param S_mRNA: fraction of mRNA found in soma
    :param z: fraction of somatically synthesized protein transported to the dendrites
    :param x: dendritic distance
    :param D_P: diffusion coefficient for protein
    :param v_P: transport velocity for protein
    :param k_P: degradation rate for protein
    :param beta_P: translation rate
    :param D_R: diffusion coefficient for mRNA
    :param v_R: transport velocity for mRNA
    :param k_R: degradation rate for mRNA
    :param beta_R: transcription rate of mRNA
    :return: steady protein level
    '''
    dendritic_protein = steady_dendritic_protein(x,D_P, v_P, k_P, beta_P,D_R, v_R, k_R, beta_R)
    lamda_P = (sqrt((v_P**2)+4*D_P*k_P)-v_P)/(2*D_P)
    somatic_protein = (exp(-(lamda_P*x)))*(beta_P*beta_R*lamda_P)/(k_R*k_P)
    return ((S_mRNA*z*somatic_protein) + (1-S_mRNA)*dendritic_protein)

def Eerf(alpha, beta, gamma):
    return exp(2 * alpha * beta) * erf(alpha * gamma + beta / gamma) + exp(-2 * alpha * beta) * erf(
        alpha * gamma - beta / gamma)


def Iminus12(a, b, A, B):
    return 1 / 2 * sqrt(pi) / sqrt(A) * (Eerf(sqrt(A), sqrt(B), sqrt(b)) - Eerf(sqrt(A), sqrt(B), sqrt(a)))


def Iplus12(a, b, A, B):
    return -1. / A * (sqrt(b) * exp(-A * b - B / b) - sqrt(a) * exp(-A * a - B / a)) + 1. / (4 * A) * sqrt(pi) / sqrt(
        A) * (Eerf(sqrt(A), sqrt(B), sqrt(b)) - Eerf(sqrt(A), sqrt(B), sqrt(a))) + 1. / (2 * A) * sqrt(pi) * sqrt(B) * (
                       Eerf(sqrt(B), sqrt(A), sqrt(a ** (-1))) - Eerf(sqrt(B), sqrt(A), sqrt(b ** (-1))))


def translational_burst(x,x0,t,t0,sigma,a1,a2,b1,b2,DP,kP,v_P):
    '''
    :param x: point of interest
    :param x_0: position of translational burst
    :param t: time of interest
    :param t_0: time of translational burst (start?)
    :param T_max:
    :param delta: width of translational burst?
    :param a1: ?
    :param a2: ?
    :param b1: ?
    :param b2: ?
    :param beta_b: amplitude of translational burst/translational rate during the burst
    :param D_P: protein diffusion coefficient
    :param k_P: degradation rate for protein
    :param v_P: transport velocity for protein
    :return:
    '''
    etaP = kP + (v_P**2)/(4*DP)
    muP = v_P/(2*DP)
    A1 = (etaP-a2)/(4*DP)
    A2 = (etaP-b2)/(4*DP)
    P = 0
    if t>t0:
        P = sqrt(sigma) * exp(sigma * muP ** 2 / 2 + muP * (x - x0) + sigma * kP / (4 * DP)) * 1. / 16 * DP ** (-2) * (
                a1 * exp(-a2 * (sigma / (4 * DP) + t - t0)) * (
                    Iminus12(sigma, sigma + 4 * DP * (t - t0), A1, (sigma * muP + 2 * (x - x0)) ** 2 / 4) * (
                        4 * DP * (t - t0) + sigma) - Iplus12(sigma, sigma + 4 * DP * (t - t0), A1,
                                                             (sigma * muP + 2 * (x - x0)) ** 2 / 4)) + b1 * exp(
            -b2 * sigma / (4 * DP) - b2 * (t - t0)) * (
                            -Iminus12(sigma, sigma + 4 * DP * (t - t0), A2, (sigma * muP + 2 * (x - x0)) ** 2 / 4) * (
                                4 * DP * (t - t0) + sigma) + Iplus12(sigma, sigma + 4 * DP * (t - t0), A2,
                                                                     (sigma * muP + 2 * (x - x0)) ** 2 / 4)));
    return (P.real)

def tf(t,t0,a1,a2,b1,b2):
   if t>t0:
       return (t-t0)*(a1*exp(-a2*(t-t0))-b1*exp(-b2*(t-t0)))
   else:
       return 0

def nonrepeat_random_pos(burst_per_sec, total_length):
    list = []
    while len(list) < burst_per_sec:
        temp = np.random.randint(0, total_length)
        if temp not in list:
            list.append(temp)
    return list

def random_pos(burst_per_sec, total_length):
    return np.random.randint(0,total_length,burst_per_sec)

def burst_history_generator_nonrepeat(burst_per_sec,total_length,time_length,time_gap):

    #generate a history of burst over time: structure of dictionary
    #key = time point; values = corresponding positions of burst
    history = {}
    for t in range(0,time_length,time_gap):
        history[t] = nonrepeat_random_pos(burst_per_sec,total_length)
    return history

def burst_history_generator(burst_per_sec,total_length,time_length,time_gap):

    #generate a history of burst over time: structure of dictionary
    #key = time point; values = corresponding positions of burst
    history = {}
    for t in range(0,time_length,time_gap):
        history[t] = random_pos(burst_per_sec,total_length)
    return history

def poisson_burst_history_nonrepeat(ave_burst_per_sec,total_length,time_length):
    '''
    :param ave_burst_per_sec:
    :param total_length:
    :param time_length: in seconds
    :return: burst history where average number of bursts per sec is a poisson distribution + non repeat location
    '''
    burst_count = np.random.poisson(ave_burst_per_sec,time_length)
    history = {}
    for t in range(0, time_length):
        history[t] = nonrepeat_random_pos(burst_count[t],total_length)
    return(history)

def poisson_burst_history(ave_burst_per_sec,total_length,time_length):
    '''
    :param ave_burst_per_sec:
    :param total_length:
    :param time_length: in seconds
    :return: burst history where average number of bursts per sec is a poisson distribution
    '''
    burst_count = np.random.poisson(ave_burst_per_sec,time_length)
    history = {}
    for t in range(0, time_length):
        history[t] = random_pos(burst_count[t],total_length)
    return(history)

def superposition(burst_history,t_of_interest,x_of_interest,sigma=0.50,a1=1,a2=1/300,b1=0.1,b2=1/300,dP=0.243798,kP=math.log(2,math.e)/(6.6*24*60*60),vP=0):
    #find the total protein concentration at location x_of_interest, time t_of_interest by summing up all component bursts
    total = 0
    # sum up all burst in a given time
    for burst_time in burst_history.keys():
        #only sum bursts events that happen strictly before time of interest
        if burst_time < t_of_interest:
            #sum up bursts in all positions
            for burst_position in burst_history[burst_time]:
                total += translational_burst(x_of_interest,burst_position,t_of_interest,burst_time,sigma,a1,a2,b1,b2,dP,kP,vP)
    return total

def modified_superposition(burst_history,t_of_interest,x_of_interest,threshold = 10**(-4),sigma=0.50,a1=1,a2=1/300,b1=0.1,b2=1/300,dP=0.243798,kP=math.log(2,math.e)/(6.6*24*60*60),vP=0):
    #find the total protein concentration at location x_of_interest, time t_of_interest by summing up all component bursts
    total = 0
    # sum up all burst in a given time
    for burst_time in burst_history.keys():
        #only sum bursts events that happen strictly before time of interest
        if int(burst_time) < t_of_interest:
            #sum up bursts in all positions
            for burst_position in burst_history[burst_time]:
                addition = translational_burst(x_of_interest,int(burst_position),t_of_interest,int(burst_time),sigma,a1,a2,b1,b2,dP,kP,vP)
                if addition > threshold:
                    total += addition
        else:
            break
    return total

def poisson_interval_history_generator(burst_rate,total_length,hours_total,path):
    '''
    :param burst_rate:
    :param total_length:
    :param hours_total:
    :param path: path to contain the history, not including final /
    :return:
    '''
    tmax = hours_total * 3600  # 1hrs = 3600s

    # average interval between two bursts
    ave_interval = 1 / burst_rate
    # number of interval within the interested duration
    interval_num = int(tmax / ave_interval)

    # add 0 to get the loop started
    burst_timestamp = [0]
    burst_location = []

    # divide into chunk of 10,000 data points, times 10 then divide by 10 to avoid integer
    chunk_num = int(interval_num / 10000)
    print("Start generating history...")
    for j in range(0, chunk_num):
        interval = np.around(np.random.poisson(ave_interval * 10, 10000) / 10, decimals=1)
        position = np.random.randint(0, total_length, 10000)
        burst_location.extend(position)
        for i in range(0, 10000):
            burst_timestamp.append(np.around(burst_timestamp[-1] + interval[i], decimals=1))

    burst_timestamp.remove(0)

    np.save(f"{path}/timestamp_{burst_rate}Hz", burst_timestamp)
    np.save(f"{path}/location_{burst_rate}Hz", burst_location)
    print(f"Burst timestamp and location successfully saved into {path}")

def accumulative_protein(time_arr, location_arr, interested_t, interested_x, curr_index, curr_protein, filename,protein_halflife_days=6.6):
    burst_duration = 10
    dP = 0.243798
    vP = 0
    LifTP = protein_halflife_days * 24 * 60 * 60 #days -> seconds
    kP = math.log(2, math.e) / (LifTP)
    sigma = 0.50
    a1 = 1
    a2 = 1 / (burst_duration / 2 * 60)
    b1 = 0.1
    b2 = 1 / (burst_duration / 2 * 60)
    protein = curr_protein
    index = curr_index
    max_t = interested_t * 3600  # inseconds
    with open(filename, 'a') as f:
        f.write(f"t\tprotein\n")
        for i in range(index, len(time_arr)):
            print(f"Called with index {i}")
            if time_arr[i] > max_t:
                break
            else:
                protein += translational_burst(interested_x, location_arr[i], max_t, time_arr[i], sigma, a1, a2,
                                                  b1, b2, dP, kP, vP)
                f.write(f"{i}\t{protein}\n")

def obtain_result(result_file_path):
    full_data = pd.read_csv(result_file_path,sep="\t")
    protein_level = full_data.iloc[-1]['protein']
    return protein_level