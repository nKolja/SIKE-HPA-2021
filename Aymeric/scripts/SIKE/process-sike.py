PTLEN = 346
KEYLEN = (218-1+7)//8

def print_key_bitstring(key):
    """
    Prints the key bits in same order they are procesesd by simpleserial-sike.
    The value of 2 below coresponds to the '0b' skip, -1 reverse order.
    """
    print(bin(int.from_bytes(key, byteorder="little", signed=False))[:2:-1])