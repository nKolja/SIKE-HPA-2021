
#Guessing the key using the traces of the file capture.py and the C function key_generator

import os
import time
import math
import random
import numpy as np
from common import *
from array import array
from scipy import stats

#plotting dependencies
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

#progress bar dependencies
from tqdm import tqdm








#Used for setting starting parameters
#Now is obsolete since we use all the parameters that we have, but at the beginning of our experiments we were not able to use all of them

#INDICES are the indices of the mp_addfast which are used
# There are 10 of them out of 22.
# These are indices at positions 3, 7, 14, 15, 16, 17, 18, 19, 20, 21
# However they are now saved at indices 0, 1, 2, ..., 9
# by changing INDICES we choose which ones to use in correlation computation
# INDICES provide "verticality" to our attack

#RANGE_NT is used for selecting the types of arrays to be correlated with power traces
# there are three types
# 0 - input array 1 of mp_addfast
# 1 - input array 2 of mp_addfast
# 2 - output array of mp_addfast
# each one of these consists of 14 32-bit words.
# Set RANGE_NT to [0], [1], [2], [0,1], [0,2], [1,2] or [0,1,2] depending on the types of arrayes to be used

#WORDS is used to select which words will be used in correlations
# There are 14 words for each input or output array
# They are indexed from 0 to 13

#i is the current bit which is being guessed

#N is the length of the trace

#SPIKES is the array of locations of the spikes within the trace, for each of the 14 words
# and for each of the three input/output arrays.
# The spikes depend on the traces, and on their size
# They have been precomputed and divided into three sets depending on the sampling rate and the length of the trace.

def read_positioning_data(i, N):

    INDICES = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    if (i >= 4): #This is due the the first points containing too many zeros in the second coordinate [xQ:1], [xQP:1], [xP:1]. This leads to correlations in mp_addfast_9 being uncomputable. This is why it is skipped until the first bit swap. This does not influence the attack in a meaningful way.
        INDICES.append(9)

    RANGE_NT = [0, 1, 2]
    # RANGE_NT = [0, 1]
    # RANGE_NT = [2]

    WORDS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    # WORDS = [0, 3, 4, 7, 8, 11, 12, 13]
    # WORDS = [0, 3, 4]
    # WORDS = [0]



    # N = 1000, Test_all_mp/
    if (N==1000):
        SPIKES = [[70, 74, 78, 82, 146, 150, 154, 158, 222, 226, 230, 234, 298, 302],\
        [90, 94, 98, 102, 166, 170, 174, 178, 242, 246, 250, 254, 310, 314],\
        [125, 133, 137, 145, 202, 210, 214, 222, 278, 286, 290, 298, 330, 334]]

    # N = 629, Test_0/
    elif(N == 629):
        SPIKES = [[78, 84, 88, 92, 160, 162, 166, 170, 238, 243, 249, 252, 320, 324],\
        [100, 104, 108, 112, 181, 185, 189, 193, 260, 264, 268, 272, 337, 337],\
        [139, 147, 150, 159, 219, 227, 230, 238, 299, 307, 310, 319, 354, 358]]

    else:
        SPIKES = [[N//2]*14]*3
        # SPIKES = [[912, 946, 981, 777, 1556, 1599, 1619, 1420, 2200, 2242, 2284, 2063, 2707, 2775],\
        # [912, 946, 981, 777, 1556, 1599, 1619, 1420, 2200, 2242, 2284, 2063, 2707, 2775],\
        # [1091, 1145, 1198, 1251, 1725, 1781, 1838, 1894, 2368, 305, 312, 2537, 2808, 2910]]


    return [INDICES, WORDS, SPIKES, RANGE_NT]














def main(n, num_keys_used, corrtype, one_key):
    t0 = time.time()

    #Secret key which we're guessing
    sk = bin(int("5b17afd2e60369c5f7bcfafae8ade3e00ca644eaaccfecbc8a1841",16))[2:].zfill(217)
    #Guessed values saved in BITS
    BITS = ""
    #Alice public keys which are used
    KEYS = genkeys(num_keys_used, one_key)


    # Set dir_num = "0" for the test with filtered traces. Tot 460 alice pk available, each trace is of length 629, all traces for all 217 bits of each of the 460 measurements are available.
    # Set dir_num = "1" for test with Chipwhisperer measurements. Tot 100 alice pk available, each trace is of length 1000, only the traces for first 6 bits of the 100 measurements are available.
    dir_num = "0"

    #Depth of the depth search
    #Complexity grows exponentially
    #When it is equal to 0 then only two correlations are computed at each step
    #Reccomended to set depth_search_step to 1
    depth_search_steps = 1

    test_directory = "Test_" + dir_num + "/"
    N = 629 if (dir_num == "0") else 5000 if (dir_num == "3") else 1000

    #Set second value to 1 if you wish to plot -- plotting = [0, 1]
    #Otherwise set plotting = [0, 0]
    plotting = [0, 0]


    #Generate the starting points Q, Q-P, P of the montgomery ladder and write them in /points
    start_points(KEYS, test_directory)





    # Plotting preprocessing
    if plotting[1]:
        fig = plt.figure(constrained_layout=True, figsize=(18,9))
        gs = fig.add_gridspec(depth_search_steps+1, 2**(depth_search_steps+1+1))
        axs = [[0 for a in range(2**(depth_search_steps+1))] for b in range(depth_search_steps+1)]
        for j in range(depth_search_steps+1):
            for k in range(2**(j+1)):
                start = 2**(depth_search_steps+1-1-j) - 1
                shift = 2*start
                axs[j][k] = fig.add_subplot(gs[j, start + (shift+2)*k : start+(shift+2)*k + 2])
                axs[j][k].set_title("bit guess " + bin(k)[2:].zfill(j+1))
                axs[j][k].axis(ymin=-0.8,ymax=0.8)
                axs[j][k].set_xticks(np.arange(0, N, N//20))
                axs[j][k].set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
                plt.setp(axs[j][k].get_xticklabels(), visible=False)




    # Depth search preprocessing
    # Precompute first depth_search_steps-1 steps
    max_step = min(depth_search_steps, n - 1)
    spikes = [[0]*(2**(c+1)) for c in range(max_step + 1)]

    for end_step in range(0, max_step):
        positioning_data = read_positioning_data(end_step, N)

        for newbits in [bin(c)[2:].zfill(end_step)[:end_step] for c in range(2**end_step)]:
            
            pcc, _ = correlate(0, end_step, KEYS, BITS + newbits, positioning_data, corrtype, test_directory, N)                
            _, _, mins = spike_check(pcc, N)

            spikes[end_step][int(newbits + "0", 2)] = mins[0]
            spikes[end_step][int(newbits + "1", 2)] = mins[1]

            if(plotting[1]):
                axs[end_step][int(newbits + "0", 2)].grid(True, linewidth=0.15, color='gray', linestyle='-')
                axs[end_step][int(newbits + "1", 2)].grid(True, linewidth=0.15, color='gray', linestyle='-')
                axs[end_step][int(newbits + "0", 2)].plot(pcc[0])
                axs[end_step][int(newbits + "1", 2)].plot(pcc[1])



    # Start of the actual algorithm
    # For each (of the 2^max_bits) length "max_bits" binary combination 
    # we compute the associated elliptic curve points, their Hamming weights and correlations
    # and add them all up. Select the one with the strongest correlation.

    #EXAMPLE: Say we use depth search of level 2, so we try all the 3-bit combinations (000,001,010,...)
    # Suppose the guess is 010
    # Then we compute elliptic curve points in case the first bit of the sk is 0, and the associated Hamming weights
    # and correlations.
    # Then we do the same if the second bit of the sk is 1
    # Then we do the same if the third bit of the sk is 0
    # We take all three the correlations that we computed and we add them all up.
    # We obtain a single correlation value for the bit combination "010". 
    # If this is the strongest correlation, we set our first guess of the secret key bit to be 0, and move the 
    # window one step forward, now guessing the next three bits, and re-using the computed correlations for the first two.

    start_step = 0
    while start_step < n:
        # Bestspike is the value of the strongest correlation (minimum of all corellation mins)
        # Bestbits is the associated bit(s) guess
        bestspike = depth_search_steps + 1
        bestbits = ""

        # max_step is the depth of the depth search. This is the step at which we compute the correlations
        # It is controlled never to be greater than bit size of Bobs private key
        # Positioning data is used to generate parameters for correlation computation
        max_step = min(depth_search_steps + 1, n - start_step) -  1
        positioning_data = read_positioning_data(start_step + max_step, N)

        for newbits in [bin(c)[2:].zfill(max_step)[:max_step] for c in range(2**max_step)]:

            #For each "max_step" length binary combination, compute the 
            # elliptic curve points, the associated hamming weights and
            # the correlations. Ouptut PCC - the correlation coefficient
            pcc, _ = correlate(start_step, start_step + max_step, KEYS, BITS + newbits, positioning_data, corrtype, test_directory, N)                
            
            #Used to find the strongest correlation (it is negative, thus we search the min)
            #bit = the bit guess associated to the stronger correlation \in {0,1}
            #m = the value of the strongest correlation
            #mins = the mins for both bit guesses
            bit, m, mins = spike_check(pcc, N)

            #Update this binary combination with the associated correlation mins
            # spikes holds the value of the strongest correlation for each guess
            spikes[max_step][int(newbits + "0", 2)] = mins[0]
            spikes[max_step][int(newbits + "1", 2)] = mins[1]
            newmin = m

            #update the correlation for this "max_step" bit combination
            # with the correlation values from previously computed correlations, i.e. for previous bits
            for step in range(max_step):
                newmin += spikes[step][int(newbits[:step+1], 2)]

            #Set the best guess to be the one with the strongest (min) correlation
            if (newmin < bestspike):
                bestspike = newmin
                bestbits = newbits + str(bit)

            if(plotting[1]):
                axs[max_step][int(newbits + "0", 2)].grid(True, linewidth=0.15, color='gray', linestyle='-')
                axs[max_step][int(newbits + "1", 2)].grid(True, linewidth=0.15, color='gray', linestyle='-')
                axs[max_step][int(newbits + "0", 2)].plot(pcc[0])
                axs[max_step][int(newbits + "1", 2)].plot(pcc[1])


        

        next_maxstep = min(max_step+1, n - start_step - 1)

        #Once a bit guess has been made, update all the point triples to accomodate for this choice
        new_points_and_hamming_weights(start_step, start_step, test_directory, KEYS, BITS + bestbits)
        update_points(start_step + 1, int(bestbits[0]))
        spikes = update_spikes(spikes, int(bestbits[0]), next_maxstep)
        start_step += 1
        BITS += bestbits[0]



        if(plotting[1]):
            plt.show()

            # plt.savefig("./plots/newplot.png")        
            plt.close(fig)

            if(start_step < n):
                fig = plt.figure(constrained_layout=True, figsize=(18,9))
                gs = fig.add_gridspec(next_maxstep, 2**(next_maxstep+1))
                axs = update_plots(axs, fig, gs, N, int(bestbits[0]), next_maxstep)


    #EXIT instantly if a single wrong guess is made to save time
    maxcorrect = 0
    while maxcorrect < n:
        if BITS[maxcorrect] == sk[maxcorrect]:
            maxcorrect += 1
        else:
            break

    print("Key: " + str(one_key).zfill(3) + ". " + str(maxcorrect).zfill(3) + " out of 217 bits correct. Time: " + "{0:5.2f}".format(time.time()-t0), "seconds")



    return (BITS[: n - depth_search_steps] == sk[: n - depth_search_steps])





#INPUT
number_of_bits = 217


#corrtype = correlation type - either Pearson correlation or Spearman correlation. Pearson is better.

#num_keys_used = number of Alice public keys / power traces which is used in the attack
#    1 = horizontal/ephemeral attack
#   >1 = vertical attack with a static Bob private key.

#num_test = number of tests to do.
#   If num_keys_used=1 then we do num_test tests with public keys selected consecutively starting from public_key = 0, to key=num_tests - 1
#   If num_keys_used>1 then tests are done with randomly chosen public keys from the pool of available ones

#number_of_bits = the number of bits at which we stop guessing.
#   A single Bob private key has 217 bits.
#   If dir_num is selected to be "0", then the traces of all 217 bits are available.
#   We can stop sooner. If dir_num is selected to be "1", then only 6 bits are available due to the limited capacity of the ChipWhisperer.
    
for corrtype in ["pearson"]:#, "spearman"]:
    for num_keys_used in [1]:
        correct = 0
        for num_tests in range(460):
            correct += main(number_of_bits, num_keys_used, corrtype, num_tests)
        print(corrtype + ",", num_keys_used, "keys used, ", number_of_bits, "bits tested")
        print(str(correct) + "% correctness\n")









#    sk = "10110100 11010001 11101011 10010111 11001110 10000000\
# 00101101 01000111 11011111 01111011 10111110 10111110\
# 00101110 01101010 10001111 00001111 01100000 11001010\
# 01000100 10101110 01101010 11100110 01101111 01111010\
# 10100010 00110000 00000100 00000001"
#    sk = "".join([ a[::-1] for a in sk.split(" ")])[:217]






    ## O L D   S E A R C H
    ## N O  B A C K T R A C K I N G

 

     # if plotting[1]:
    #     fig, axs = plt.subplots(2,  1, sharex=False, sharey=True, figsize=(12,6))
    #     # plt.xticks(np.arange(0, N, N//6))
    #     # plt.setp(axs, xticks = np.arange(0, N, N//6), yticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
    #     plt.ylim(-40, 40)
    #     # fig.suptitle('Number of keys used = ' + str(num_keys_used), fontsize = 16)
    #     plotting[0] = axs


    # for step in range(0, n):

    #     # fig = plt.figure(figsize=(18,9))
    #     # plt.xticks(np.arange(0, N, N//40))
    #     # plt.yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
    #     # plt.ylim(-0.8, 0.8)
    #     # fig.suptitle('Number of keys used = ' + str(num_keys_used), fontsize = 16)

    #     positioning_data = read_positioning_data(step, N)
    #     pcc, trace = correlate(step, step, KEYS, BITS, positioning_data, corrtype, test_directory, N)
    #     # bit, spike, m = key_bit(pcc, step, plotting, N) 
    #     bit, spike, m = spike_check(pcc, N) 

    #     # bit = int(sk[step])
    #     BITS += str(bit)
    #     update_points(step+1, bit)

    #     axs[0].grid(True, linewidth=0.15, color='gray', linestyle='-')
    #     axs[0].plot([x/1000 for x in trace[0]], color="blue", label="bit 0", linewidth=0.60)
    #     axs[0].legend(['Pre-processing'], loc="upper right")
    #     axs[0].set_xlabel('Points [Pt]')
    #     axs[0].set_ylabel('Power consumption [V]')

    #     # plt.savefig("./plots/keys/key_" + str(one_key).zfill(3) + "/bit_" + str(step).zfill(3) + ".png")        
    #     # plt.close(fig)


    # if plotting[1]:
    #     # plt.show()

    #     plt.savefig("./plots/pre-post-proc.png")        
    #     plt.close(fig)







