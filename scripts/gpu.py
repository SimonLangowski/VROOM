from copy import deepcopy
from random import seed
from precompute_montgomery import UNUSED, gen_baselines, get_mont_factor, get_mont_factor2, result_func_baseline, correct_func_baseline
from basic_benchmarker import run_benchmark
from rns_helpers import *
from rearrangement import *
from precompute_matrix import *
from python_reduction import *
from debug import *

# Arguments determining problem size parameters
MODULUS_ARG = ("--modulus", int, 0, "Modulus for arithmetic")
MODBITS_ARG = ("--modbits", int, 256, "Choose random k-bit modulus (if no modulus provided)")
BOUND_ARG = ("--bound", int, 0, "Maximum number of additions performed between multiplications")

# Arguments controlling variations and parameters of algorithm
MATMUL_M_ARG = ("--m", int, 16, "Batch size, must be multiple of 16")
Q_CORRECT_ARG = ("--qcorrect", bool, False, "Apply second correction step, completing modular reduction")

variation_choices = {
    -1: "MontRedundantBaselineReduce", # Baseline using redundant form with Montgomery (no tensor cores)
    0: "BasicRNSReduce", # Basic RNS based algorithm
    1: "DualRNSReduce", # RNS reduction using two matrices
    2: "MontRNSReduce", # big integer Montgomery algorithm emulated with RNS reduction
    3: "MontRNSReduceGeneral", # Version where matrix precomputation is independent of modulus
    4: "MontRNSReduceVec", # Version with vector cores
}
VARIATION_ARG = ("--variation", int, 2, f"Variation choice: {variation_choices}")

EPILOGUE_ARG = ("--epilogue", int, 1, f"Epilogue: 1 - Correct, 2 - Full reduction, 0 - Implicit (when applicable)")

reduction_method_choices = {
    -1: "None", # No reduction
    0: "Psuedomersenne", # Moduli of the form 2**32 - s for small s
    1: "MontRedundant30", # Montgomery reduction with redundant form, with 30 bit moduli
    2: "Mont32", # Montgomery reduction with 32 bit moduli
}
REDUCTION_METHOD_ARG = ("--reduction-method", int, 2, f"Choice of elementwise reduction method: {reduction_method_choices}")

# Arguments for testing
application_choices = {
    0: "bench_modmult2"
}
APP_ARG = ("--app", int, 0, f"Application choice: {application_choices}")
POW_ARG = ("--pow", int, 1, "Iterations for repeated multiplication microbenchmark")
DET_ARG = ("--det", bool, False, "Use a determinstic seed") # I couldn't consistently spell deterministic on the command line so this got shortened

TPI_ARG = ("--btpi", int, 4, "Base threads used to store each number (power of 2)")

STRIDE_ARG = ("--stride", int, 1, "Parallelization factor (power of 2)")

VECTOR = 4

# Populate names from numbers for argument options
def translate_arg_dicts(args):
    args.reduction_method_name = reduction_method_choices[args.reduction_method]
    args.variation_name = variation_choices[args.variation]
    args.app_name = application_choices[args.app]
    # preprocessor requires int not bool
    args.debug = 1 if args.debug else 0
    if args.debug == 1:
        args.iterations = 1
        args.det = True
    if args.det:
        seed(5)

# Fill in the template parameters as required by the variation
def construct_class(args):
    assert((args.epilogue >= 0) and (args.epilogue <= 2))
    if args.variation == 0:
        args.mult_class = f"{args.variation_name}<{args.moduli2}, RegisterFragment, {args.reduction_method_name}, {args.epilogue}>"
    elif args.variation == -1:
        args.mult_class = f"{args.variation_name}<{args.limbs}>"
    elif args.variation == VECTOR:
        args.mult_class = f"{args.variation_name}<{args.moduli1}, {args.moduli2}, {args.mpt1}, {args.mpt2}, {args.tpi}, {args.stride}, {args.reduction_method_name}>"
    elif args.variation >= 1:
        args.mult_class = f"{args.variation_name}<{args.moduli1}, {args.moduli2}, RegisterFragment, {args.reduction_method_name}, {args.epilogue}>"

# Create a precomputed matrix following the specification
def gen_matrix_pre(precompute, args, vector=False, q_correct=False):
    # precompute_ints = gen_precompute_wide(moduli_in_info, target, moduli_out_info, premult, postmult, 4, 8, 1, 32, 32)
    # print(str(precompute_ints).replace(" ", "\n").replace("[", "").replace("]", ""), file=open("debug.txt", "w"))
    if precompute == None:
        return np.zeros(1, dtype=np.uint32)
    precompute_ints = precompute.precompute
    precompute_mats = precompute_to_matrix_i(precompute_ints)
    if vector:
        return basic_pack(precompute_mats, args)
    # print(precompute_mats, file=open("debug.txt", "w"))
    # return pack_b(precompute_mats, q_correct)

def to_numpy_arr(args, int_list):
    num_entries = len(int_list)
    mx, nx = int_list[0]
    tpi = args.btpi
    skip_cols = 2 if args.variation == VECTOR else 0
    m, n = len(mx) + skip_cols, len(nx) + skip_cols
    # Is MPT a power of 2, or just TPI?
    args.mpt1 = ceil_div(m, tpi)
    args.mpt2 = ceil_div(n, tpi)
    padded_m, padded_n = tpi*args.mpt1, tpi*args.mpt2
    am = np.zeros((num_entries, padded_m), dtype=np.uint32)
    an = np.zeros((num_entries, padded_n), dtype=np.uint32)
    if args.variation == 4:
        args.limbs = padded_n
    for i, e in enumerate(int_list):
        ml, nl = e
        am[i, skip_cols:m] = ml.store()
        an[i, skip_cols:n] = nl.store()
    return am, an

# Convert the input into RNS form
def input_to_rns(args, a_py, b_py, alg):
    a_conv = [alg.to_mont_avx(a) for a in a_py]
    b_conv = [alg.to_mont_avx(b) for b in b_py]
    a1_np, a2_np = to_numpy_arr(args, a_conv)
    b1_np, b2_np = to_numpy_arr(args, b_conv)
    a1_f = flatten_a32(args, a1_np)
    b1_f = flatten_a32(args, b1_np)
    a2_f = flatten_a32(args, a2_np)
    b2_f = flatten_a32(args, b2_np)
    print(f"{(a1_f.nbytes + a2_f.nbytes) // 4} flattened u32s per input")
    return a1_f, b1_f, a2_f, b2_f

# Pack the moduli info into an array
def prepare_moduli(args, moduli_info):
    if args.reduction_method == 0:
        mt = [2**32 - x for x in moduli_info]
        while len(mt) % 4 != 0:
            mt.append(0)
        return np.array(mt, dtype=np.uint32)
    else:
        moduli = [m for m in moduli_info]
        mont_factors = [pow(m, -1, 2**32) for m in moduli]
        # np.pad((0, 0), (2, (len(moduli) + 2) % 4)))
        if args.variation == 4:
            for i in range(2):
                moduli.insert(0, 0)
                mont_factors.insert(0, 0)
        # * args.mpt?
        while len(moduli) % args.btpi != 0:
            moduli.append(0)
            mont_factors.append(0)
        return strided(np.array([moduli, mont_factors], dtype=np.uint32), args)
    
def padded(a, b):
    return b * ceil_div(a, b)

def gen_func_rns_red(args):
    translate_arg_dicts(args)
    
    modulus = set_get_modulus(args)
    
    length = args.total
    if args.variation < 0:
        return gen_baselines(args, modulus)
    
    args.tpi = args.stride * args.btpi
    
    _, a_py = random_array(length, modulus, None)
    _, b_py = random_array(length, modulus, None)

    # Maybe more complicated in 2t^2 case
    reducers = [PythonReduce, MontgomeryReduce, MontgomeryReduce]
    algorithms = [BasicRNS, None, IntRNS2, None, IntRNS2]

    target = modulus
    if args.variation == VECTOR:
        params = Params.gpuv(reducers[args.reduction_method])
    else:
        params = Params.gput(reducers[args.reduction_method])
        assert(args.tpi == 4)
    if args.reduction_method == 1:
        params.set30()
    alg = algorithms[args.variation](target, params)

    moduli1_info = alg.m1
    moduli2_info = alg.m2
    precompute1 = alg.r1
    precompute2 = alg.r2
    args.precision = params.precision
    b_pre1 = gen_matrix_pre(precompute1, args, args.variation == VECTOR)
    b_pre2 = gen_matrix_pre(precompute2, args, args.variation == VECTOR)

    args.moduli1 = len(moduli1_info)
    args.moduli2 = len(moduli2_info)
    moduli1_arr = prepare_moduli(args, moduli1_info)
    moduli2_arr = prepare_moduli(args, moduli2_info)

    args.limbs = args.moduli1 if args.moduli2 == 0 else padded(args.moduli1, args.tpi) + args.moduli2
    args.moduli_in = args.limbs

    # Option to convert on GPU instead of here?
    a1_f, b1_f, a2_f, b2_f = input_to_rns(args, a_py, b_py, alg)

    print(f"m1 {args.moduli1}, m2 {args.moduli2}, m in {args.moduli_in}, bpre1 {b_pre1.nbytes}, bpre2 {b_pre2.nbytes}, a1 {a1_f.nbytes}, a2 {a2_f.nbytes}, b1 {b1_f.nbytes}, b2 {b2_f.nbytes} mpt1 {args.mpt1} mpt2 {args.mpt2} tpi {args.tpi}, bpt {args.stride}")

    # print(b_pre1.tolist(), file=open("debugn.txt", "w"))

    args.shared_size = args.tpi * 4
    return [a1_f, b1_f, b_pre1, moduli1_arr, a2_f, b2_f, b_pre2, moduli2_arr], [a_py, b_py, alg]

def result_func_rns_red(args):
    construct_class(args)
    if args.variation < 0:
        return result_func_baseline(args)
    length = args.total
    limbs = args.btpi * args.mpt2 #args.limbs # stride added by flatten function
    args.resultshape = (length, limbs)
    return flatten_a32(args, np.zeros((length, limbs), dtype=np.uint32)), args.shared_size, (args.blocksize * args.m) // args.tpi

def correct_func_rns_red(args, input, output):
    a_py, b_py, alg = input
    modulus = args.modulus
    correct = [(a_py[i] * pow(b_py[i], args.pow, modulus)) % modulus for i in range(len(a_py))]
    if (args.variation < 0):
        return correct_func_baseline(args, correct, output)
    # print(np.array2string(a_np))
    u = unflatten_32(args, output, args.resultshape)
    rows, cols = u.shape
    assert(rows == len(correct))
    failed = 0
    last_fail = -1
    skip_cols = 2 if args.variation == VECTOR else 0
    for i in range(rows):
        # Are we returning m2 only? Or both m1 and m2?
        c = alg.from_mont_avx(AVXVector([int(x) for x in u[i, skip_cols:skip_cols+len(alg.m2)]]))
        if(c % modulus != correct[i]):
            failed += 1
            if last_fail == -1 or failed < 10:
                print(i, c % modulus, correct[i])
                last_fail = i
        else:
            if last_fail >= 0:
                print(i, c % modulus, correct[i])
                last_fail = -1
        if ((i == 0) and (args.debug)):
            # generate trace
            aavx, bavx = alg.to_mont_avx(a_py[i]), alg.to_mont_avx(b_py[i])
            DEBUG[0] = True
            for _ in range(args.pow):
                cavx = alg.mulreduce(*aavx, *bavx)
                aavx = deepcopy(cavx)
            DEBUG[0] = False
            ans = alg.from_mont_avx(cavx[1])
            assert(ans == correct[i])
        if (args.variation == VECTOR):
            # Only checking one entry for latency versions
            break
    if failed > 0:
        if (args.variation != VECTOR):
            print(f"{rows - failed}/{rows}")
        return False
    else:
        print("Passed!")
        return True
    

def run_rns_latency_benchmark(args=None, _=None):
    return run_benchmark("latencygpu", "../gpu/benching/loadingvec.cuh", "MULT_APP",
        [MODULUS_ARG, MODBITS_ARG, POW_ARG, REDUCTION_METHOD_ARG, APP_ARG, DET_ARG, EPILOGUE_ARG, TPI_ARG, MATMUL_M_ARG, VARIATION_ARG, STRIDE_ARG], gen_func_rns_red, result_func_rns_red, correct_func_rns_red, args, False)
