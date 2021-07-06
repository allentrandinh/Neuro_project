from mpmath import *
import numpy as np
import math
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
    :return:
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
    :return:
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