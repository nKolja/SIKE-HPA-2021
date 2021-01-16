
sys.path.append('../')
from setup import simpleserial, chipwhisperersetup

CRYPTO_TARGET='AES'
SCOPETYPE='OPENADC'
PLATFORM='CWLITEARM'

# Path to AES compiled code
fw_path = "/home/lemonade/Code/chipwhisperer-5.2.1/hardware/victims/firmware/simpleserial-aes/simpleserial-aes-{}.hex".format(PLATFORM)

# Connects to chipwhisperer
(target, scope) = chipwhisperersetup(fw_path, CRYPTO_TARGET, SCOPETYPE, PLATFORM)

# AES 1st round
scope.adc.samples = 5000
