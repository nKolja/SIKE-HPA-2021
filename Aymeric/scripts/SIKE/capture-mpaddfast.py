import os
import time
import random
import datetime

### Constants ###
# Random seed
SEED = "Extra sugar, extra salt, extra oil and MSG"

# Add of two field numbers (in p434) => closest is 56*8 = 440, hence 2x56
PTLEN = 112

# Total number of traces to collect.
N_traces = 1000

# Total number of samples (20 us)
scope.adc.samples = int(29.54*20)

# Keep log
filename = datetime.datetime.now().strftime("data/%Y-%m-%d_%H-%M-%S_mpaddfast-pts.txt")
f = open(filename, 'w')
print(f"Opened {filename} !")

# Start capture
try:
    plaintexts = []
    traces = []
    random.seed(SEED)
    for i in range(N_traces):
        # Generate random plaintext
        pt = int.to_bytes(random.randint(0, 2**(PTLEN*8)-1), byteorder="big", length=PTLEN)
        f.write(''.join([hex(x)[2:].zfill(2) for x in pt]) + '\n')

        # Arm oscilloscope
        scope.arm()
        target.flush()

        # Capture
        target.simpleserial_write('q', pt)
        ret = scope.capture()
        if ret:
            print('Timeout happened during acquisition')
            continue

        # Clear buffer (if any)
        num_char = target.in_waiting()
        rd = ""
        while num_char > 0:
            rd += target.read(timeout=10)
            time.sleep(0.01)
            num_char = target.in_waiting()

        # Save text + trace
        plaintexts += [pt]
        traces += [scope.get_last_trace()]

        # Print progress
        if i % (N_traces//10) == 0:
            print(f"Captured {i}/{N_traces}...")

    # Finished
    print(f"Captured {N_traces}/{N_traces}...")

finally:
    print(f"Closing {filename}...")
    f.close()
