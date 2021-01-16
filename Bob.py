#This code is used for generations of Bob's keys. 
#The secret key is stored in "bob_keys/bob_sk" and the public key in "bob_keys/bob_pk"
#A binary file containing only the sk of bob is saved in "./bob_keys/only_bob_sk"

import numpy as np
import os

KEYLEN = (218-1+7)//8
C = "./bob_key"




# THIS PART OF CODE IS OBSOLETE.
# IT CAN STILL BE EXECUTED AND IT DOES NOT INTERFERE WITH KEY GENERATION
# IF YOU WISH TO GENERATE A SECRET 

if (0):
    Bob_sk = int.to_bytes(int("aa"*KEYLEN, 16), byteorder="little", length=KEYLEN)
    Bob_sk = os.urandom(KEYLEN)
    #print(Bob_sk)
    key_string = ""
    for i in range(KEYLEN):
        key_string+= hex(Bob_sk[i])[2:].zfill(2)

    # BOB's SECRET KEY HAS 218-1 = 217 BITS (SEE SIDH/SIKE PARAMETERS)
    # WE GENERATE A RANDOM NUMBER OF 224 BITS AND THEN ZERO OUT THE LAST 7 BITS
    key_string = key_string[:54] + hex((int((key_string[-2:]),16)&(0x01)))[2:].zfill(2)
    
    C += " " + str(key_string)
    # print(key_string)



# A random key will be created with the C code in src/side-channel.c which uses original SIDH code.

return_value = os.system(C)