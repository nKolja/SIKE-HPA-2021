import os
import time
import datetime
import numpy as np
import random

### Constants ###
# Random seed
SEED = "Extra sugar, extra salt, extra oil and MSG"

# Add of two field numbers (in p434) => closest is 56*8 = 440, hence 2x56
PTLEN = 112

# Total number of traces to collect.
N_traces = 1000

# Keep log
filename = datetime.datetime.now().strftime("data/%Y-%m-%d_%H-%M-%S_mpaddfast-pts.txt")
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
        target.simpleserial_write('q', pt)
        time.sleep(0.5)

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
