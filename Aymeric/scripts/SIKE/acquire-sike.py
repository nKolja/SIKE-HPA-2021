import os
import time
import datetime
import numpy as np
import random

# Total number of bits to be recorded (max = 218).
NBITS = 218

# Total number of traces to collect.
N_traces = 250

# Read Bob's secret key from precomputed data
with open("data/SIKE/bob_sk.bin", "rb") as f:
    bob_sk = f.read()

key = bob_sk[16:16+28]

assert len(key) == KEYLEN # sanity check
simpleserial(target, 'k', key)
simpleserial(target, 'n', int.to_bytes(NBITS, length=1, byteorder="little", signed=False))

# Read plaintexts from precomputed data
alice_rootfiles = "data/SIKE/alice_pk_"

# Start capture
plaintexts = []
for i in range(N_traces):
    alice_file = alice_rootfiles + str(i).zfill(5) + ".bin"
    if not os.path.isfile(alice_file):
        print(f"ERROR: {alice_file} does not exist!")
        break

    # Read entire plaintext from file
    with open(alice_file, "rb") as f:
        pt = f.read()

    assert len(pt) == PTLEN # Sanity check

    # Capture
    now = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    print(f"{now}\tsending... {pt.hex()}")
    target.simpleserial_write('p', pt)
    time.sleep(4*60 + 15) # Extra long time to make sure that procedure has finished

    # Clear buffer (if any)
    num_char = target.in_waiting()
    rd = ""
    while num_char > 0:
        rd += target.read(timeout=1)
        time.sleep(0.01)
        num_char = target.in_waiting()

    # Save text + trace
    plaintexts += [pt]

    # Print progress
    if i != 0 and i % (N_traces//10) == 0:
        print(f"Captured {i}/{N_traces}...")

# Finished
print(f"Captured {N_traces}/{N_traces}...")
