import os
import re
import scipy.signal
import numpy as np
from utils import load_trace

def single_invert_trace(tr):
    return list(map(lambda x: -x, tr))

def single_wavelet_denoise(tr, N=1, wavelet='db1'):
    dn_tr = pywt.dwt(tr, wavelet)[0]
    if N > 1:
        return single_wavelet_denoise(dn_tr, N=N-1, wavelet=wavelet)
    else:
        return dn_tr

def list_files(traces_dir="/home/lemonade/Data/lecroy/SIKE/218bits/"):
    p = re.compile(".*--([0-9]+)\\.trc")
    listfiles = []
    for file in os.listdir(traces_dir):
        m = p.match(file)
        if m:
            listfiles += [file]
    return sorted(listfiles)

def convert_all_traces(traces_dir="/home/lemonade/Data/lecroy/SIKE/218bits/2020-11-27_250x218x22/", out_dir="./data/SIKE/",
                       length=5000, Nseq=22, Nbits=218, height=20000, Nwv=3, wavelet='db3'):
    # 1. List all files
    listfiles = list_files(traces_dir=traces_dir)
    Ntraces = len(listfiles)
    print(f"Found: {Ntraces}\n")

    # 2. Find peaks in first subtrace
    file = listfiles[0]
    print(f"Looking for peak in {file}...")
    fulltrace = load_trace(os.path.join(traces_dir, file))
    tr = fulltrace[:length]
    pk = scipy.signal.find_peaks(tr, height=height)
    if len(pk[0]) < 1:
        print(f"ERROR: no peak found!! Displaying trace now:")
        print(tr)
        return (fulltrace, [])
    x_ref = pk[0][0]
    print(f"Reference peak: tr[{x_ref}] = {tr[x_ref]}\n")

    # 3. Process all files and align on peak found in first subtrace
    bit = 0 # from 0 to Nbits
    idx = 0 # from 0 to Ntraces
    for file in listfiles:
        print(f"Loading {file}...")
        fulltrace = load_trace(os.path.join(traces_dir, file))
        seqlength = len(fulltrace)//Nseq
        for i in range(Nseq):
            tr = fulltrace[seqlength*i:seqlength*i+length]
            #print(f"\tfulltrace[{seqlength*i}:{seqlength*i+length}] ({len(tr)}, ", end="")
            assert len(tr) == length, "ERROR: slicing of fulltrace failed!"
            pk = scipy.signal.find_peaks(tr, height=height)
            if len(pk[0]) < 1:
                print(f"ERROR: no peak found!! Displaying trace now:")
                print(tr)
                return (fulltrace, [])
            x = pk[0][0]
            delta = x_ref - x
            if delta == 0:
                print(f"delta={delta}, no shifting... tr[{x_ref}+0={x}] = {tr[x]}")
            elif delta > 0:
                print(f"delta={delta}, shifting right... tr[{x_ref}-{delta}={x}] = {tr[x]}")
                tr = delta*[0] + tr[:-delta]
            else:
                delta = abs(delta)
                print(f"delta={delta}, shifting left... tr[{x_ref}+{delta}={x}] = {tr[x]}")
                tr = tr[delta:] + delta*[0]
            tr = single_wavelet_denoise(tr, N=Nwv, wavelet=wavelet)
            #print(f"\tres = {len(tr)}")
            assert len(tr) == (length + 4 + 8 + 16 + 8)//(2**3), "ERROR: wavelet transform did not shrink correctly"
            dest_dir = f"{out_dir}traces_{str(idx).zfill(3)}{os.path.sep}bit{str(bit).zfill(3)}{os.path.sep}"
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            np.savetxt(f"{dest_dir}trace_{str(i).zfill(2)}.txt", tr, fmt="%d", delimiter=" ", newline=" ", header="", footer="", comments="", encoding="latin1")
        print(f"done!\n")

        # Update bit and trace index
        bit += 1
        if bit >= Nbits:
            bit = 0
            idx += 1