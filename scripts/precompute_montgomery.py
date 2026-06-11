
import numpy as np
from rns_helpers import random_array, word_len

# A array of size 1 that will be allocated needlessly
UNUSED = np.array([0], dtype=np.uint32) 

def get_mont_factor(modulus, r=2**32):
    return pow(-modulus, -1, r)

def get_mont_factor2(modulus, r=2**32):
    return pow(modulus, -1, r)

def to_mont(x, modulus, limbs, limb_size=32):
    return (x << (limb_size*limbs)) % modulus

def to_montt(x, modulus, mont_factor):
    return (x * mont_factor) % modulus

def to_limbs(val, limb_count):
    assert(2**(32 * limb_count) > val)
    mask = 2**32 - 1
    return [((val >> (32 * i)) & mask) for i in range(limb_count)]

def from_limbs(limbs):
    return sum([(int(limbs[i]) << (32*i)) for i in range(len(limbs))])

def restride_baseline(arr, args):
    elements = args.total
    limbs = args.limbs
    out = np.zeros((elements // args.blocksize, limbs, args.blocksize), dtype=np.uint32)
    for i in range(elements):
        tid = i % args.blocksize
        bid = i // args.blocksize
        for k in range(limbs):
            out[bid, k, tid] = arr[i, k]
    return out

def unstride_baseline(out, args):
    elements = args.total
    limbs = args.limbs
    arr = np.zeros((elements, limbs), dtype=np.uint32)
    for i in range(elements):
        tid = i % args.blocksize
        bid = i // args.blocksize
        for k in range(limbs):
            arr[i, k] = out[bid, k, tid]
    return arr

def gen_baselines(args, modulus):
    if not args.skip_correct:
        assert(modulus.bit_length() % 32 != 0) # Redundant form requires one bit
    args.limbs = word_len(modulus, 32)
    args.moduli_in = 4 * args.limbs # Hacky but it's not a fragment so multiply by 4
    # Unused variables required for compilation
    args.moduli1 = 0
    args.moduli2 = 0
    args.shared_size = 0
    args.precision = 0
    a_np, a_py = random_array(args.total, modulus, None)
    b_np, b_py = random_array(args.total, modulus, None)
    mont_factor = get_mont_factor(modulus)
    mont_factor_np = np.array([mont_factor], dtype=np.uint32)
    mod_np = np.array(to_limbs(modulus, args.limbs), dtype=np.uint32)
    return [restride_baseline(a_np, args), restride_baseline(b_np, args), mont_factor_np, mod_np, UNUSED, UNUSED, UNUSED, UNUSED], [a_py, b_py, None, None]

def result_func_baseline(args):
    length = args.total
    limbs = args.limbs
    return restride_baseline(np.zeros((length, limbs), dtype=np.uint32), args), 0, args.blocksize * (args.m // 8) 

def correct_func_baseline(args, correct, output):
    modulus = args.modulus
    elements = args.total
    limbs = args.limbs
    arr = unstride_baseline(output, args)
    failed = 0
    for i in range(args.total):
        v = from_limbs(arr[i, :])
        if (to_mont(v, modulus, limbs) != correct[i]):
            if failed < 10:
                print(i, to_mont(v, modulus, limbs), correct[i])
            failed += 1
    return failed == 0