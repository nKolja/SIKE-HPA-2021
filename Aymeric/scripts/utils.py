import lecroy2sigrok
import os
import re
import scipy.signal
import numpy as np
import random
import pywt

def prinths(s):
    print(''.join([hex(ord(c))[2:].zfill(2) for c in s]))

def printba(b):
    print(''.join([hex(c)[2:].zfill(2) for c in b]))

def read_pts(filename):
    with open(filename, 'r') as f:
        plaintexts = f.read().split('\n')
    return [int.to_bytes(int(pt, 16), byteorder="big", length=len(pt)//2) for pt in plaintexts if pt]

def write_pts(plaintexts, dir="./data/"):
    assert os.path.isdir(dir), f"Error: {dir} is not an existing directory"
    np.savetxt(f"{dir}{os.path.sep}plaintexts.txt", [t.hex() for t in plaintexts], fmt="%s", delimiter="\n", header="", footer="", comments="", encoding="latin1")

def load_trace(file, start=0, end=None):
    data, dataFormat, horizInterval = lecroy2sigrok.parseFile(file)
    end = min(end, len(data)//2) if end else len(data)//2
    tr = []
    for i in range(start, end):
        tr += [int.from_bytes(data[2*i:2*(i+1)], byteorder="little", signed=True)]
    return tr

def load_all_traces(traces_dir, start=0, end=None):
    p = re.compile(".*--([0-9]+)\\.trc")
    traces = []
    listfiles = []
    for file in os.listdir(traces_dir):
        m = p.match(file)
        if m:
            listfiles += [file]
    for file in sorted(listfiles):
        print(f"{file}...", end="")
        traces += [load_trace(os.path.join(traces_dir, file), start=start, end=end)]
        print(" done!")
    return traces

def read_all_traces(traces_file):
    assert os.path.isfile(traces_file), f"Error: {traces_file} is not an existing file"
    return np.loadtxt(traces_file, dtype=int)

def write_all_traces(all_traces, dir="./data/", start=0, end=None):
    assert os.path.isdir(dir), f"Error: {dir} is not an existing directory"
    if end:
        assert start < end, f"Error: invalid range: [{start}:{end}]"
    end = min(end, len(all_traces[0])) if end else len(all_traces[0])
    np.savetxt(f"{dir}{os.path.sep}traces_{len(all_traces)}x{end-start}.txt", [t[start:end] for t in all_traces], fmt="%d", delimiter=" ", header="", footer="", comments="", encoding="latin1")

def invert_traces(traces):
    return [list(map(lambda x: -x, t)) for t in traces]

def align_peaks(traces, height):
    # Reference peak
    pk = scipy.signal.find_peaks(traces[0], height=height)
    ref_x = pk[0][0] # peaks will be aligned on this coordinate

    aligned_traces = [traces[0]]
    print(f"x_ref = {ref_x}, peaks will be aligned on this point");
    for t in traces[1:]:
        pk = scipy.signal.find_peaks(t, height=height)
        if len(pk[0]) < 1:
            print(f"ERROR: no peak found!!")
            continue
        x = pk[0][0]
        delta = ref_x - x
        if delta == 0:
            print(f"delta={delta}, no shifting... {t[x]} <> {traces[0][ref_x]}")
            aligned_traces += [t]
        elif delta > 0:
            print(f"delta={delta}, shifting right... {t[x]} <> {traces[0][ref_x]}")
            aligned_traces += [delta*[0] + t[:-delta]]
        else:
            print(f"delta={delta}, shifting left... {t[x]} <> {traces[0][ref_x]}")
            delta = abs(delta)
            aligned_traces += [t[delta:] + delta*[0]]
    return aligned_traces

def fft_traces(traces, sample_rate=1E6):
    traces_f = [list(map(lambda x: abs(x), t[:len(traces[0])//2])) for t in np.fft.rfft(traces)]
    freqs = np.fft.rfftfreq(len(traces[0]), d=1./sample_rate)
    assert len(freqs) - 1 == len(traces_f[0])
    return (freqs[1:], traces_f)

# Hamming weights of all bytes (from 0x00 to 0xff)
HW = [bin(n).count("1") for n in range(256)]

def wavelet_denoise(traces, N=1, wavelet='db1'):
    dn_tr = list(map(lambda t: pywt.dwt(t, wavelet)[0], traces))
    if N > 1:
        return wavelet_denoise(dn_tr, N=N-1, wavelet=wavelet)
    else:
        return dn_tr

def corr_HW(hw,tr):
    """
    @input hw List of Hamming Weight to correlate with trace
    @input tr List of Power Traces
    """
    n = len(tr) # Number of traces
    m = len(tr[0]) # Length of traces

    pcc = []
    # For each sample point
    for p in range(m):
        T = [tr[i][p] for i in range(n)] # Gets each sample from all traces at said point
        pcc.append(np.corrcoef(T, hw)[0,1])

    return pcc

def expand_seed(seed, N=100, length=16):
    random.seed(seed)
    return [int.to_bytes(random.randint(0, 2**(length*8)-1), byteorder="big", length=length) for _ in range(N)]
