#COMMON FUNCTIONS USED IN CORRELATION COMPUTING SCRIPTS


import os
# import sys
import time
import math
import random
import numpy as np
from copy import copy
from array import array
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec




# Computes the minimum of a a pair of arrays [ar0, ar1] and outpus the index of the array with the lowest min
# the value of the min and the pair of mins of each array
def spike_check(pcc, N):
    # mins = [min(p[N//2 - 10 : N//2 + 10]) for p in pcc]
    mins = [min(p) for p in pcc]

    spike = False

    m = min(mins)
    i = mins.index(m)

    return i, m, mins





# OBSOLETE - used for computing the min and plotting.
def key_bit(pcc, step, axs, N):

    if(axs[1]):
        axs[0][step].grid(True, linewidth=0.15, color='gray', linestyle='-')
        axs[0][step].plot(pcc[0], color="red", label="bit 0", linewidth=0.5)
        axs[0][step].plot(pcc[1], color="blue", label="bit 1", linewidth=0.5)
        axs[0][step].legend(['Step ' + str(step)], loc="upper right")

    return spike_check(pcc, N)





# Used to updated the values of the spikes in a depth search once a bit guess has been made
def update_spikes(spikes, bit, newrange):
    out = [[0]*2**(c+1) for c in range(newrange)]
    for r in range(newrange - 1):
        for x in range(2**(r+1)):
            out[r][x] = spikes[r+1][x + bit * (2**(r+1))]

    return out


# Used to update plots in a depth search once a bit guess has been made
def update_plots(axs, fig, gs, N, bit, newrange):
    axsnew = [[0 for a in range(2**(newrange))] for b in range(newrange)]

    for r in range(newrange - 1):
        for x in range(2**(r+1)):
            start = 2**(newrange-1-r) - 1
            shift = 2*start

            axsnew[r][x] = fig.add_subplot(gs[r, start + (shift+2)*x : start + (shift+2)*x + 2])
            axsnew[r][x].set_title("bit guess " + bin(x)[2:].zfill(r+1))
            axsnew[r][x].axis(ymin=-0.8,ymax=0.8)
            axsnew[r][x].set_xticks(np.arange(0, N, N//20))
            axsnew[r][x].set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
            axsnew[r][x].plot(axs[r+1][x + bit * (2**(r+1))].get_lines()[0].get_ydata())
            axsnew[r][x].grid(True, linewidth=0.15, color='gray', linestyle='-')
            plt.setp(axsnew[r][x].get_xticklabels(), visible=False)

    for x in range(2**(newrange)):
        axsnew[newrange-1][x] = fig.add_subplot(gs[newrange-1, 2*x : 2*x + 2])
        axsnew[newrange-1][x].set_title("bit guess " + bin(x)[2:].zfill(newrange))
        axsnew[newrange-1][x].axis(ymin=-0.8,ymax=0.8)
        axsnew[newrange-1][x].set_xticks(np.arange(0, N, N//20))
        axsnew[newrange-1][x].set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
        plt.setp(axsnew[newrange-1][x].get_xticklabels(), visible=False)


    return axsnew




# Used to generate keys
# If number of keys is one, then the key indexed "one_key" from input is provided
# If number of keys is 100, 250, or 460 (which is maximum for Tests 0, 1 and 2) then all keys are given
# Otherwise nk random keys are provided a the range which is determined based on nk.

def genkeys(nk, one_key):
    if(nk == 1):
        return [one_key]

    if(nk == 100):
        return range(100)

    if(nk == 250):
        return range(250)
    
    if nk == 460:
        return range(460)

    MAX = 460 if nk > 250 else 250 if nk > 100 else 100

    L = 0
    out = []
    while (L < nk):
        k = random.randint(0, MAX - 1)
        if k not in out:
            out.append(k)
            L += 1
    return out






# Transposes an array of arrays
def transpose(M):
    return list(map(list, zip(*M)))








# Used to read Hamming weights/distances file from disk
# These files are written by the "hamming_weight_computation" c script
def read_weights(i, KEYS, WORDS, INDICES):
    """
    @input i the indice of the ith steps
    """
    nk = len(KEYS)
    nw = len(WORDS)
    ni = len(INDICES)

    input_file = open("weights/bit_" + str(i).zfill(3), 'rb')
    hw = array('i')
    hw.frombytes(input_file.read())
    hw = [list([ hw[k*2*10*3*14 + b*10*3*14 + mp_addfast_index*3*14 + t*14 + w] for mp_addfast_index in INDICES for k in range(nk)]) for b in range(2) for t in range(3) for w in WORDS]
    #hammwt[bit*3*nw + t*nw + w][mp_addfast_index * nk + key]


    # OUTDATED

    # original_stdout = sys.stdout
    # with open('filenameall' +str(i)+ '.txt', 'w') as f:
    #     sys.stdout = f # Change the standard output to the file we created.
    #     for k in range(nk):
    #         for w in WORDS:
    #             print(hw[0*3*nw + 2*nw + w][k], "\t", hw[1*3*nw + 2*nw + w][k])
    #     sys.stdout = original_stdout # Reset the standard output to its original value


    return hw





# Used to read the power traces from disk
# Only power traces of keys associated to KEYS and mp_addfasts associated to NDICES are read
def read_traces(step, KEYS, INDICES, tot_keys, tracelength, test_directory, N):
    nk = len(KEYS)
    ni = len(INDICES)

    input_file = open(test_directory + "Traces/mp_" + str(step), 'rb')
    out = array('d')
    out.frombytes(input_file.read())
    out = [list(out[(mp_addfast_index*tot_keys + k)*tracelength : (mp_addfast_index*tot_keys + k + 1) * tracelength][:N]) for mp_addfast_index in INDICES for k in KEYS]
  
    return transpose(out)



#Calls the c script to compute hamming weights and to update points in the montgomery triple
def new_points_and_hamming_weights(start_step, end_step, traces_directory, KEYS, BITS):
    command = "./hamming_weight_computation " + str(len(KEYS)) + " " + " ".join([str(k) for k in KEYS]) + " " + str(start_step) + " " + str(end_step) + " " + traces_directory + " " + BITS
    os.system(command)

    return





#Calls functions to read the traces, compute hamming weights, updates the points of the montgomery triple and
# then reads the associated Hamming weights
def read_data(start_step, end_step, KEYS, BITS, INDICES, WORDS, traces_directory, N):
    if (traces_directory == "Test_2/") or (traces_directory == "Test_new/"):
        tot_keys = 250
        tracelength = 629
    elif (traces_directory == "Test_3/"):
        tot_keys = 460
        tracelength = 5000
    elif (traces_directory == "Test_0/"):
        tot_keys = 460
        tracelength = 629
    else:
        tot_keys = 100
        tracelength = 1000

    traces = read_traces(end_step, KEYS, INDICES, tot_keys, tracelength, traces_directory, N)
    new_points_and_hamming_weights(start_step, end_step, traces_directory, KEYS, BITS)
    weights = read_weights(end_step, KEYS, WORDS, INDICES)

    return traces, weights





#Used to compute the Pearson correlation coefficient
#Spearman correlation is also possible if corrtype is set to "spearman"
#Inputs are the starting Montgomery triple (S,T,U) which is at the start_step of the three point ladder
#This triple is updated up to end_step by using the BITS provided in input
#The power traces at end_step are read.
#For both cases of bit value of the "end_step" bit, a new triple of points is computed, together with the
# associated hamming weights.
#Then the correlation between traces and hamming weights is computed for both cases (bit guess = 0,1)
#Pair of correlation is returned
#Traces are also returned (only used for plotting purposes)
def correlate(start_step, end_step, KEYS, BITS, positioning_data, corrtype, test_directory, N):

    INDICES, WORDS, SPIKES, RANGE_NT = positioning_data

    nw = len(WORDS)
    nk = len(KEYS)
    ni = len(INDICES)
    nt = len(RANGE_NT)

    traces, hamming_weights = read_data(start_step, end_step, KEYS, BITS, INDICES, WORDS, test_directory, N)
    #traces[trace]               [mp_addfast_index * nk + key]
    #hammwt[bit*nt*nw + t*nw + w][mp_addfast_index * nk + key]

    if (corrtype == "spearman"):
        traces = [stats.rankdata(a) for a in traces]
        hamming_weights = [stats.rankdata(a) for a in hamming_weights]

    elif(corrtype != "pearson"):
        print("Correlation type unknown!")
        exit()

    PCC = np.corrcoef(hamming_weights, traces)
    pcc = [0,0]

    for b in range(2):
        pcc[b] = [ sum([PCC[b*3*nw + t*nw + w][2*3*nw + ((tr + SPIKES[t][WORDS[w]] - N//2)%N)] for t in RANGE_NT for w in range(nw)])/(nt*nw) for tr in range(N)]

    return pcc, transpose(traces)








#Used to create the starting montgomery triple Q,Q-P,P
#This is done for each alice public key in KEYS
def start_points(KEYS, keys_directory):
    start_points = "./make_starting_points " + " ".join([str(k) for k in KEYS]) + " " + keys_directory
    os.system(start_points)

    return






#Once a guess for a bit has been made, the montgomery point triple (S,T,U) is updated
#The points for bit guess that we made (bit_val) are renamed and saved
#The other points ( !bit_val ) are removed
def update_points(i, bit_val):
    correct_point_path = "points/pt_" + str(i).zfill(3) + "_" + str(bit_val)
    wrong_point_path = "points/pt_" + str(i).zfill(3) + "_" + str(1-bit_val)
    new_path = "points/pt_" + str(i).zfill(3)
    cmd1 = "mv " + correct_point_path + " " + new_path
    cmd2 = "rm " + wrong_point_path
    os.system(cmd1)
    os.system(cmd2)

    return



