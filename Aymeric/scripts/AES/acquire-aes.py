import random
import time
import datetime
import numpy as np

SEED = "Extra sugar, extra salt, extra oil and MSG"
KEYLEN = 128//8
PTLEN = 128//8
N_traces = 100

# Generate fixed key
key = int.to_bytes(int("A3112233445566778899AABBCCDDEEFF", 16), byteorder="big", length=KEYLEN)

# Write key to card
simpleserial(target, 'k', key)

# Keep log
filename = datetime.datetime.now().strftime("data/%Y-%m-%d_%H-%M-%S_AES-pts.txt")
f = open(filename, 'w')
print(f"Opened {filename} !")

# Start capture
try:
    plaintexts = []
    random.seed(SEED)
    for i in range(N_traces):
        # Generate random plaintext
        pt = int.to_bytes(random.randint(0, 2**(PTLEN*8)-1), byteorder="big", length=PTLEN)
        f.write(''.join([hex(x)[2:].zfill(2) for x in pt]) + '\n')

        # Capture
        target.simpleserial_write('p', pt)
        time.sleep(0.2)

        # Clear buffer (if any)
        num_char = target.in_waiting()
        rd = ""
        while num_char > 0:
            rd += target.read(timeout=10)
            time.sleep(0.01)
            num_char = target.in_waiting()

        # Save text + trace
        plaintexts += [pt]

        # Print progress
        if i % (N_traces//10) == 0:
            print(f"Captured {i}/{N_traces}...")

    # Finished
    print(f"Captured {N_traces}/{N_traces}...")

finally:
    print(f"Closing {filename}...")
    f.close()
