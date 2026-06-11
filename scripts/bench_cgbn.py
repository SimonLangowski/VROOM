from random import randint, seed
import numpy as np
from basic_benchmarker import run_benchmark
from rns_helpers import set_get_modulus, word_len, word_reinterpret
from gpu import MODBITS_ARG, MODULUS_ARG, MATMUL_M_ARG, TPI_ARG, POW_ARG

def random_array(length, modulus, mod_words=None, cache_key=None):
    if mod_words == None:
        mod_words = word_len(modulus, 32)
    np_arr = np.zeros((length, mod_words), dtype=np.uint32)
    py_arr = []
    for i in range(length):
        v = randint(0, modulus)
        py_arr.append(v)
        np_arr[i, :] = word_reinterpret(v, 32, mod_words)  
    return np_arr, py_arr

def gen_func_randints(args):
    length = int(args.total)
    modulus = set_get_modulus(args)
    # For PYREPLACE_LIMBS
    args.limbs = word_len(modulus, 32)
    a_np, a_py = random_array(length, modulus)
    b_np, b_py = random_array(length, modulus)
    # print(a_py[736] * b_py[736])
    return [a_np, b_np], [a_py, b_py]

def int_reconstruct(x, word_length):
    return sum([int(d) << (word_length * i) for (i, d) in enumerate(x)])

def to_py_array(np_array, words, word_length=32):
    return [int_reconstruct(np_array[i*words:(i+1)*words], word_length) for i in range(len(np_array) // words)]

# Not certain which version of Montgomery CGBN uses
def compute_mont_factor(modulus):
    return pow(-modulus, -1, 2**32)

def to_mont(x, modulus, modbits):
    return (x << modbits) % modulus

def gen_func_cgbn(args):
    seed(5)
    args.tpi = args.btpi
    [a_np, b_np], [a_py, b_py] = gen_func_randints(args)
    modulus = args.modulus
    args.mont_factor = compute_mont_factor(modulus)
    print(args.mont_factor)
    mod_np = np.zeros(args.limbs, dtype=np.uint32)
    mod_np[:] = word_reinterpret(modulus, 32, args.limbs)
    modbits = 32 * args.limbs
    for (i, b) in enumerate(b_py):
        b_np[i, :] = word_reinterpret(to_mont(b, modulus, modbits), 32, args.limbs)
    args.mult_class = "CGBN"
    args.variation_name = "CGBN"
    args.reduction_method_name = "montmul"
    return [a_np, b_np, mod_np], [a_py, b_py]

def result_func_cgbn(args):
    modulus = int(args.modulus)
    mod_words = word_len(modulus, 32)
    length = int(args.total)
    return np.zeros(length * mod_words, dtype=np.uint32), 0, args.blocksize // args.tpi

def correct_func_cgbn(args, input, output):
    a_py, b_py = input
    modulus = args.modulus
    c_t = to_py_array(output, args.limbs)
    assert(len(c_t) == len(a_py))
    for i in range(len(a_py)):
        if ((a_py[i] * pow(b_py[i], args.pow, modulus)) % modulus != c_t[i] % modulus):
            print(i, (a_py[i] * b_py[i]) % modulus, c_t[i])
            print(word_reinterpret(a_py[i], 32, args.limbs), word_reinterpret(b_py[i], 32, args.limbs), word_reinterpret(modulus, 32, args.limbs))
            return False
    return True

def run_cgbn_benchmark(args=None, _=None):
    return run_benchmark("cgbn_baseline", "../gpu/algorithms/baselines/cgbn_mul.cuh", "MULT_CGBN", [MODULUS_ARG, MODBITS_ARG, TPI_ARG, POW_ARG, MATMUL_M_ARG], gen_func_cgbn, result_func_cgbn, correct_func_cgbn, args)

if __name__ == "__main__":
    run_cgbn_benchmark()