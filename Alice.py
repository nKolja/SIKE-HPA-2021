#This code is used for generations of Alice's keys in function of Bob's public key. 

#The code reads Bob's public key (from the file bob_keys/bob_pk) and
#saves multiple Alice public keys in alice_keys/alice_pk_***** the *s representing the index of said public key of Alice.

#INPUT NEEDED FOR NUMBER OF PUBLIC KEYS WE WISH TO CREATE

import os
x = input('How many Alice keys do you wish to create?\n')
NB_Alice_key = int(x)
# print(NB_Alice_key)
C = "./alice_key "+ str(NB_Alice_key)
return_value = os.system(C)
