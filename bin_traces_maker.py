#Guessing the key using the traces of the file capture.py and the C function key_generator

import os
import time
import math
import random
import numpy as np
from scipy import stats
from array import array

#plotting dependencies
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#progress bar dependencies
from tqdm import tqdm







def genkeys(nk):
    L = 0
    out = []
    if nk == 100:
        return range(100)
    if nk == 250:
        return range(250)
    if nk == 460:
        return range(460)

    while (L < nk):
        k = random.randint(0,99)
        if k not in out:
            out.append(k)
            L += 1

    return out




def make_bin_traces(bit, KEYS, INDICES, tot_keys, test_directory, N):
    nk = len(KEYS)
    ni = len(INDICES)
    


    # READ TRACES
    out = []
    for ind in INDICES:
        for k in KEYS:
            filename = test_directory + "Traces_old/traces_" + str(k).zfill(3) + "/bit" + str(bit).zfill(3) + "/trace_" + str(ind).zfill(2) + ".txt"

            f = open(filename, 'r')
            out_ind = f.readlines()
            f.close()

            for i in out_ind:
                out.append([float(j) for j in i.split(" ")[:N]])



    # WRITE TRACES AS BINARY FILES
    output_file = open(test_directory + 'Traces/mp_' + str(bit), 'wb')
    float_array = array('d', [x for y in out for x in y])
    float_array.tofile(output_file)
    output_file.close()



    return out








def test_function(n, num_keys_used):

    dir_num = input("Please insert test directory number: ")

    test_directory = "Test_" + dir_num + "/"
    N = 629 if (dir_num == "0") else 5000 if (dir_num == "3") else 1000

    KEYS = genkeys(num_keys_used)
    INDICES = [3, 7, 14, 15, 16, 17, 18, 19, 20, 21]
    INDICES = range(10)

    for i in tqdm(range(n)):

        traces = make_bin_traces(i, KEYS, INDICES, len(KEYS), test_directory, N)




    return 1




numer_of_bits = 217
numer_of_keys = 460

test_function(numer_of_bits, numer_of_keys)







    # out = []

    # for ind in INDICES:
    #     filename = "Test_all_mp/Traces/mp_" + str(bit) + "_" + str(ind)

    #     f = open(filename, 'r')
    #     out_ind = f.readlines()
    #     f.close()

    #     for i in out_ind:
    #         out.append([float(j) for j in i.split(" ")])

    
    # for k in KEYS:
    #     out[5*nk + k] = out[5*nk + k][4:] + out[5*nk + k][:4]
 


    # out_all

    # filename = "Test_all_mp/Traces_compact/mp_" + str(bit)

    # f = open(filename, 'r')
    # out_all = f.readlines()
    # f.close()

    # out = [out_all[mp_addfast_index * tot_keys + KEYS[key]] for mp_addfast_index in range(ni) for key in range(nk)]
    # out = [[float(j) for j in i.split(" ")] for i in out]


    # for k in KEYS:
    #     test = out[0*nk + k][:400]
    #     output = [0]*len(INDICES)

    #     for j in range(1, len(INDICES)):
    #         compare = out[j*nk + k][(400 - 1)::-1]

    #         comparison = list(np.convolve(test, compare, "full"))
    #         m = max(comparison, key=abs)

    #         shift = (400-1) - comparison.index(m)
    #         output[j] = shift
    #         out[j*nk + k] = out[j*nk + k][shift:] + out[j*nk + k][:shift]

    #     print(output)




    # input_file = open('Test_all_mp/Traces_bin/mp_' + str(bit), 'rb')
    # out2 = array('d')
    # out2.frombytes(input_file.read())

    # out2 = [list(out2[(mp_addfast_index*tot_keys + k)*1000 : (mp_addfast_index*tot_keys + k + 1) * 1000]) for mp_addfast_index in range(ni) for k in KEYS]
    

    # for i in range(len(INDICES)):
    #     print(out2[i*nk] == out[i*nk])
        


    # plt.plot(out[0])
    # plt.plot(out[5*nk])
    # plt.show()

    # exit()
