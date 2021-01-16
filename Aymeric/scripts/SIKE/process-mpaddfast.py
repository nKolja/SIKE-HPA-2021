import numpy as np

PTLEN = 112 # 2x56 bytes (a: 56 bytes, b: 56 bytes)

def process_HW_chunks(a, b, n=4):
    c_a = [int.from_bytes(a[i:i+n], byteorder="little") for i in range(0,len(a),n)]
    c_b = [int.from_bytes(b[i:i+n], byteorder="little") for i in range(0,len(b),n)]

    res = [p[0] + p[1] for p in zip(c_a, c_b)]
    return [bin(r).count("1") for r in res]

n_chunks = 4 # one chunk = 4 bytes
res_hw = {}
for i in range(PTLEN//(2*n_chunks)):
    res_hw[i] = []

for p in plaintexts:
    res = process_HW_chunks(p[:PTLEN//2], p[PTLEN//2:], n=n_chunks)
    for i in range(PTLEN//(2*n_chunks)):
        res_hw[i] += [res[i]]

all_pcc = []
for i in res_hw:
    all_pcc += [corr_HW(res_hw[i], traces)]
    print(f"Chunk {i}")