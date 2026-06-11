from math import gcd
from random import randint
import numpy as np
from sympy import isprime, nextprime

def total_modulus(moduli):
    modulus = 1 # Compute the big modulus
    for m in moduli:
        assert(m != 0)
        modulus *= m
    return modulus

def rns_precompute(moduli):
    modulus = total_modulus(moduli)
    precomputed = []
    for m in moduli:
        rest = modulus // m # 0 mod all the other moduli
        inverse = pow(rest % m, -1, m) # factor to make 1 mod this moduli
        icrt_val = (rest * inverse) % modulus # combine
        precomputed.append(icrt_val)
    return precomputed

def rns_reconstruct(residues, moduli, precomputed):
    assert(len(residues) == len(moduli))
    assert(len(moduli) == len(precomputed))
    output = 0
    for (i, r) in enumerate(residues):
        output += precomputed[i] * r
    return output % total_modulus(moduli)

def rns_reconstruct_full(residues, moduli):
    return rns_reconstruct(residues, moduli, rns_precompute(moduli))

def to_rns(x, moduli, check=True):
    if check:
        assert(x < total_modulus(moduli))
    return [x % m for m in moduli]

def ceil_div(x, y):
    return (x + y - 1) // y


def gen_extend(target, state=-1, prev=1, modbits=52):
    t = state
    moduli = []
    modulus = 1
    total = modulus * prev
    while (modulus < target):
        t += 2
        m_candidate = (2**modbits) - t
        if gcd(m_candidate, total) == 1: # relatively prime
            modulus *= m_candidate
            total *= m_candidate
            moduli.append(m_candidate)
    return moduli, state, modulus

def word_len(x, word_len):
    return (x.bit_length() + word_len - 1) // word_len

def byte_len(x):
    return word_len(x, 8)

def word_reinterpret(x, word_length, length=-1):
    if length == -1:
        length = word_len(x, word_length)
    assert(x < (2**word_length)**length)
    return [(x >> (word_length*i)) & (2**word_length - 1) for i in range(length)]

def byte_reinterpret(x, length=-1):
    return word_reinterpret(x, 8, length)

def int_reconstruct(x, word_length):
    return sum([int(d) << (word_length * i) for (i, d) in enumerate(x)])

def concat(x, y):
    out = []
    out.extend(x)
    out.extend(y)
    return out

def redundance(val, mod):
    return (val // mod).bit_length()

def set_get_modulus(args):
    modulus = int(args.modulus)
    if modulus == 0:
        modbits = int(args.modbits)
        modulus = 4 # not prime
        while(modulus.bit_length() != modbits):
            modulus = nextprime(randint(2**(modbits - 1), (2**modbits) - 1))
        print(f"Picked {modulus.bit_length()} bit  modulus ", modulus)
        args.modulus = modulus
    return modulus

def rns_array(py_array, moduli):
    a_mat = np.zeros((len(py_array), len(moduli)), dtype=np.uint32)
    for i, v in enumerate(py_array):
        a_mat[i, :] = to_rns(v, moduli)
    return a_mat

def random_array(length, modulus, mod_words=None):
    if mod_words == None:
        mod_words = word_len(modulus, 32)
    np_arr = np.zeros((length, mod_words), dtype=np.uint32)
    py_arr = []
    for i in range(length):
        v = randint(0, modulus)
        py_arr.append(v)
        np_arr[i, :] = word_reinterpret(v, 32, mod_words)
    return np_arr, py_arr
