
sys.path.append('../')
from setup import simpleserial, chipwhisperersetup

CRYPTO_TARGET='SIKE'
SCOPETYPE='OPENADC'
PLATFORM='CWLITEARM'

# Path to SIKE compiled code
fw_path = "/home/lemonade/Code/chipwhisperer-5.2.1/hardware/victims/firmware/simpleserial-sike/simpleserial-sike-{}.hex".format(PLATFORM)

# Connects to chipwhisperer
(target, scope) = chipwhisperersetup(fw_path, CRYPTO_TARGET, SCOPETYPE, PLATFORM)

# Full SIKE trace
scope.adc.samples = 5000

### Constants ###
# Bob's private key is 218-bit long, the closest number of bytes is 28
KEYLEN = (218-1+7)//8

"""
Alice sends 6x434 bits that corresponds to two elliptic curves points:
    (R[0]->x, R[0]->z, R[1]->x, R[1]->z, R[2]->x, R[2]->z)

The closest multiple of 8 is (434 + 7)/8 = 55     (NBYTES)
The closest multiple of 32 is (434 + 31)/32 = 14  (NWORDS)

Hence, depending on memory alignment, we send either:
    (14*32)/8 = 448/8 = 56 bytes per coordinate (i.e., elements in Fp343)
    (55*8)/8 = 440/8 = 55 bytes per coordinate (i.e., elements in Fp343)

!!! Note that 346 = 330 (= 2*3*55, the three points of 2 coordinates) +
                     16 (= MSG_BYTES, the hashed message with the two ciphertexts).
"""
PTLEN = 346
