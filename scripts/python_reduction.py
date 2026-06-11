import os
from random import seed
from precompute_io import PrecomputeIO
from rns_helpers import *
from precompute_matrix import *
from precompute_montgomery import get_mont_factor, get_mont_factor2, to_montt
from debug import *


# MODBITS = 52
# MULBITS = 52
# PRECISION = 64
# INBITS = MULBITS
# INDIGITS = 1

# MODBITS = 32
# MULBITS = 32
# PRECISION = 64
# INBITS = 32
# INDIGITS = 1

# MODBITS = 32
# MULBITS = 32
# PRECISION = 64 # seems to require w = 32 not w = 8 like basic did
# INBITS = 8
# INDIGITS = 4

# MODBITS = 32
# MULBITS = 32
# PRECISION = 32 # It's just the conversion step that requires higher precision.  In Basic, the conversion k is 0 or 1 which makes sense, but in INTRNS2 it's handling the montgomery conversion.  Needs to take in 8 bit numbers to compensate.
# INBITS = 8
# INDIGITS = 4

class Params:

    def __init__(self, modbits, mulbits, precision, inbits, indigits, reducer):
        self.modbits = modbits
        self.mulbits = mulbits
        self.precision = precision
        self.inbits = inbits
        self.indigits = indigits
        self.reducer = reducer
        self.outbits = mulbits
        self.r = 1
        if reducer.IS_MONT:
            self.r = 2**self.mulbits
        AVXVector.set_params(self)

    def refresh(self):
        AVXVector.set_params(self)

    @staticmethod
    def avx(reducer):
        return Params(modbits=52, mulbits=52, precision=64, inbits=52, indigits=1, reducer=reducer)

    @staticmethod
    def gput(reducer):
        return Params(modbits=32, mulbits=32, precision=32, inbits=8, indigits=4, reducer=reducer)
    
    @staticmethod
    def gpuv(reducer):
        return Params(modbits=32, mulbits=32, precision=64, inbits=32, indigits=1, reducer=reducer)
    
    def set30(self):
        self.modbits = 30
        self.refresh()

def mask(bits):
    return (1 << bits) - 1

# For python version, slice are unneeded, because when M1+M2 is the input it is read in scalars
# And when M1+M2 is the output, we have M2=0
class AVXVector:


    def __init__(self, values, verify=True):
        self.data = values
        self.length = len(values)
        self.mulbits = AVXVector.mulbits
        self.mask = (1 << self.mulbits) - 1
        # Verify no overflow
        if verify:
            for i in range(self.length):
                assert(self.data[i] < 2**64)

    @staticmethod
    def set_params(params):
        AVXVector.mulbits = params.mulbits

    @staticmethod
    def zero(length):
        return AVXVector([0 for _ in range(length)])
    
    @staticmethod
    def scalar(val, length):
        return AVXVector([val for _ in range(length)])

    def v(self, values):
        return AVXVector(values)

    def mulhi(self, a, b):
        return self.v([self.data[i] + (((a.data[i] & self.mask) * (b.data[i] & self.mask)) >> self.mulbits) for i in range(self.length)])
    
    def mullo(self, a, b):
        return self.v([self.data[i] + (((a.data[i] & self.mask) * (b.data[i] & self.mask)) & self.mask) for i in range(self.length)])
    
    def mulhi_scalar(self, a, b):
        return self.mulhi(a, self.v([b for _ in range(self.length)]))
    
    def mullo_scalar(self, a, b):
        return self.mullo(a, self.v([b for _ in range(self.length)]))
    
    def srli(self, bits):
        return self.v([(x >> bits) for x in self.data])
    
    def slli(self, bits):
        return self.v([(x << bits) for x in self.data])
    
    def add(self, other):
        return self.v([self.data[i] + other.data[i] for i in range(self.length)])
    
    def mask_cmp_ge(self, other):
        return self.v([self.data[i] >= other.data[i] for i in range(self.length)])
    
    def mask_sub(self, mask, other):
        return self.v([self.data[i] - (mask.data[i] * other.data[i]) for i in range(self.length)])
    
    def sub(self,other):
        return self.v([self.data[i] - other.data[i] for i in range(self.length)])
    
    def mask_add(self, mask, other):
        return self.v([self.data[i] + (mask.data[i] * other.data[i]) for i in range(self.length)])
    
    def mask_cmp_gt(self, other):
        return self.v([self.data[i] > other.data[i] for i in range(self.length)])
    
    def cmp_sub(self, other):
        mask = self.mask_cmp_ge(self, other)
        return self.mask_sub(mask, other)
    
    def extract0(self):
        return self.data[0]
    
    def store(self):
        return [x for x in self.data]
    
    def __and__(self, other):
        return self.v([self.data[i] & other.data[i] for i in range(self.length)])
    
    def __str__(self):
        return str(self.data)
    
    def __len__(self):
        return len(self.data)
    
    def bl(self):
        return [x.bit_length() for x in self.data]
    
    def cmp_eq(self, bound):
        # Extraction of first lane is fast
        l = self.data[0]
        if l > bound:
            return False
        # Compare_eq_mask
        # Check if mask8 is all 1's mask
        raise NotImplementedError("todo")
    
    def __eq__(self, other):
        return self.data == other.data

class PythonReduce:

    IS_MONT = False

    def __init__(self, moduli, params, io=None):
        self.moduli = moduli
        self.mulbits = params.mulbits
    
    def reduce_full(self, hi, low):
        return AVXVector([((hi.data[i] << self.mulbits) + low.data[i]) % self.moduli[i] for i in range(len(self.moduli))])
    
    def reduce_small(self, hi, low):
        return self.reduce_full(hi, low)
    
    def reduce32(self, wide):
        return AVXVector([wide.data[i] % self.moduli[i] for i in range(self.length)])
    
# 50 bit in 52 bit montgomery should not require correction and uses R=2**52
class PythonMontgomeryReduce:

    IS_MONT = True

    def __init__(self, moduli, params, io=None):
        self.moduli = moduli
        self.r = 2**params.mulbits
        assert(self.r == params.r)
        self.r_inv = AVXVector([pow(self.r, -1, m) for m in moduli])
        self.full = PythonReduce(moduli, params=params)
        self.length = len(moduli)

    def reduce_small(self, hi, low):
        return self.reduce_full(hi, low)
    
    def reduce_full(self, hi, low):
        x = self.full.reduce_full(hi, low)
        ex_hi = AVXVector.zero(self.length).mulhi(x, self.r_inv)
        ex_lo = AVXVector.zero(self.length).mullo(x, self.r_inv)
        return self.full.reduce_full(ex_hi, ex_lo)
    
# Much trickier with AVX
# Use case is on GPU to allow bitshifts of non-32
class MontgomeryReduce:

    IS_MONT = True

    def __init__(self, moduli, params, io=None):
        self.test = PythonMontgomeryReduce(moduli, params)
        self.r = 2**params.mulbits
        assert(self.r == params.r)
        self.mont_factor = AVXVector([pow(m, -1, self.r) for m in moduli])
        self.moduli = moduli
        self.length = len(moduli)
        self.modulia = AVXVector(moduli)
        self.neg_modulia = AVXVector([2**64 - m for m in moduli]) # 64 bit accumulator negation
        self.mulbits = params.mulbits
        self.bits = params.modbits
        self.hi_bit_shift = self.bits + self.bits - self.mulbits
        self.hi_mask = AVXVector.scalar((1 << self.hi_bit_shift) - 1, len(moduli))
        self.lo_mask = AVXVector.scalar((1 << self.mulbits) - 1, len(moduli))
        self.t = AVXVector([2**(self.hi_bit_shift) % m for m in moduli])
        self.t_squared = AVXVector([(2**(2*self.bits)) % m for m in moduli])
        if io is not None:
            io.vector(self.moduli)
            io.vector(self.mont_factor.data)
            # io.vector(self.t_squared.data)
            io.vector(self.t.data)

    # 52 bit moduli
    def reduce_small(self, hi, low):
        q = AVXVector.zero(self.length).mullo(low, self.mont_factor)
        h = AVXVector.zero(self.length).mulhi(q, self.modulia)
        uh = hi.sub(h)
        if self.mulbits == self.bits:
            mask = h.mask_cmp_gt(hi)
            return uh.mask_add(mask, self.modulia)
        else:
            return uh.add(self.modulia)

    # # need hi < m, lo < 2^52
    def reduce_full_old(self, hi, low):

        hi_hi = hi.srli(self.hi_bit_shift) 
        hi_rem = hi & self.hi_mask
        # require hi_hi * t_squared < 2^52
        # required low < 2^64 - 2^52
        lo_acc = low.mullo(hi_hi, self.t_squared)

        lo_hi = lo_acc.srli(self.mulbits)
        lo_rem = lo_acc & self.lo_mask

        # require hi_rem < 2(2^52 - t)
        hi_acc = hi_rem.add(lo_hi)
        mask = hi_acc.mask_cmp_ge(self.modulia)
        hi_acc.mask_sub(mask, self.modulia)

        # Certifies that hi < 2^52 - t so quotient will fit
        return self.reduce_small(hi_acc, lo_rem)
    
    def reduce_full(self, hi, low):
        lo_hi = low.srli(self.mulbits)
        lo_rem = low & self.lo_mask
        hi_acc = hi.add(lo_hi)

        # v4 = self.test.reduce_full(hi_acc, lo_rem)
        
        hi_hi = hi_acc.srli(self.hi_bit_shift)
        hi_lo = hi_acc & self.hi_mask
        h = hi_lo.mullo(hi_hi, self.t)

        mask = h.mask_cmp_ge(self.modulia)
        h.mask_sub(mask, self.modulia)

        v = self.reduce_small(h, lo_rem)
        # v2 = self.reduce_full_old(hi, low)
        # v3 = self.test.reduce_full(hi, low)
        # assert(v2 == v3)
        # assert(v4 == v3)
        # assert(v == v3)
        return v

def ele_mult(a, b, reducer):
    debug_printvec("a", a)
    debug_printvec("b", b)
    assert(len(a) == len(b))
    ab_lo = AVXVector.zero(len(a)).mullo(a, b)
    ab_hi = AVXVector.zero(len(a)).mulhi(a, b)
    ab = reducer.reduce_small(ab_hi, ab_lo)
    debug_printvec("product", ab)
    return ab

class RNSMatrix:
    
    def __init__(self, moduli_in, moduli_out, params, target_modulus=None, premult=1, postmult=1, io=None, inbits=None, indigits=None, outdigits=1):
        if inbits == None:
            inbits = params.inbits
        if indigits == None:
            indigits = params.indigits
        precompute = gen_precompute_wide(moduli_in, target_modulus, moduli_out, premult, postmult, indigits, inbits, outdigits, params.mulbits, params.precision)
        self.precompute = precompute # For extraction via GPU stuff
        rns_values, correction, shifted_quotient_estimations, target_shifted_quotient_estimations, tqc, tp, precision_requirement, _, _, _ = precompute
        self.rns_matrix = [AVXVector(r) for r in rns_values]
        self.correction = AVXVector(correction)
        self.shifted_quotient_estimations = shifted_quotient_estimations
        self.target_shifted_quotient_estimations = target_shifted_quotient_estimations
        self.tqc = tqc
        self.tp = AVXVector(tp)
        self.precision = precision_requirement
        self.moduli_out = moduli_out
        self.inbits = inbits
        self.indigits = indigits
        self.mulbits = params.mulbits
        if io is not None:
            # print("before", io.write_counter)
            [io.vector(rv.data) for rv in self.rns_matrix]
            # print("rns mat", io.write_counter, len(self.rns_matrix), len(self.rns_matrix[0].data))
            io.vector(self.shifted_quotient_estimations)
            io.vector(self.correction.data)
            # io.vector(self.target_shifted_quotient_estimations)
            # io.scalar(self.tqc)
            # io.vector(self.tp.data)

    def rns_reduce(self, residues1, reducer, q_correct=False, acc=None):
        scalarsr = residues1.store()
        debug_printvec("a-input32", scalarsr)
        scalars = flatten([word_reinterpret(s, self.inbits, self.indigits) for s in scalarsr])
        assert(len(scalars) == len(self.rns_matrix))
        assert(len(scalars) == len(self.shifted_quotient_estimations))
        out_hi = AVXVector.zero(len(self.correction))
        out_lo = AVXVector.zero(len(self.correction))
        if acc is not None:
            a, b = acc
            out_hi = out_hi.mulhi(a, b)
            out_lo = out_lo.mullo(a, b)
        k_raw = 0
        q_raw = 0
        # debug_printvec("scalars", scalars)
        debug_printvec("sqe", self.shifted_quotient_estimations)
        sqe_hi = [x >> self.mulbits for x in self.shifted_quotient_estimations]
        debug_printvec("sqe hi", sqe_hi)
        sqe_lo = [x & mask(self.mulbits) for x in self.shifted_quotient_estimations]
        debug_printvec("sqe lo", sqe_lo)
        k_raw_lo = 0
        k_raw_hi = 0
        for (i, s) in enumerate(scalars):
            # The scalar should be combined from low and high (so there is one set1 instruction) in c code
            debug_printvec("bpre", self.rns_matrix[i])
            out_hi = out_hi.mulhi_scalar(self.rns_matrix[i], s)
            out_lo = out_lo.mullo_scalar(self.rns_matrix[i], s)
            k_raw += s * self.shifted_quotient_estimations[i]
            if q_correct:
                q_raw += s * self.target_shifted_quotient_estimations[i]
            # Only computed to compare with gpu version
            k_raw_hi += s * sqe_hi[i]
            k_raw_lo += s * sqe_lo[i]
        # k
        # debug_printvec("sqe", self.shifted_quotient_estimations)
        debug_printvec("kraw", k_raw)
        debug_printvec("k parts", [k_raw_hi >> self.mulbits, k_raw_hi & mask(self.mulbits), k_raw_lo >> self.mulbits, k_raw_lo & mask(self.mulbits)])
        k = k_raw >> self.precision
        # Maybe, given that it's vector, we just accept a premultiplication by invicrt
        # Technically only needed on second mult but could just do on both.
        # Alternatively have 2 threads for fixed point, and slightly more complicated syncronization.
        # Just pad to a whole warp and don't worry about batching
        # If 3~wide instruction for premult
        # How much is the extra syncronization?  bitshift, 32 bit sync?  64 bit sync add?
        # print(k_raw.bit_length(), k.bit_length())
        debug_printvec("k", k)
        # print(k_raw.bit_length(), k.bit_length())
        # Broadcast once k
        
        # mask = AVXVector(2**52 - 1, len(self.correction))
        # k_vec = AVXVector(k, len(self.correction))
        # mask_k = k_vec & mask
        # k_hi = k_vec.srli(52)
        # reduced_k = mask_k.mullo(k_hi, self.t_vec)
        # reduced_k.cmp_sub(self.moduli_out)

        # Then shift in AVX form is fast
        out_hi = out_hi.mulhi_scalar(self.correction, k)
        out_lo = out_lo.mullo_scalar(self.correction, k)
        out_hi = out_hi.mullo_scalar(self.correction, k >> self.mulbits)
        lp = AVXVector.zero(len(out_hi)).mulhi_scalar(self.correction, k >> self.mulbits)
        out_hi = out_hi.add(lp.slli(self.mulbits))
        # l
        if (q_correct):
            q_raw += k * self.tqc
            q = q_raw >> self.precision
            out_hi = out_hi.mulhi_scalar(self.tp, q)
            out_lo = out_lo.mullo_scalar(self.tp, q)
            out_hi = out_hi.mullo_scalar(self.tp, q >> self.mulbits)
            lp = AVXVector.zero(len(out_hi)).mulhi_scalar(self.tp, q >> self.mulbits)
            out_hi = out_hi.add(lp.slli(self.mulbits))
        if reducer == None:
            return out_hi, out_lo
        return reducer.reduce_full(out_hi, out_lo)

    
class ConvertToRNS:

    def __init__(self, premult, postmult, target, moduli_out, params, inbase, io=None):
        self.inbase = inbase
        self.indigits = word_len(target, inbase)
        self.conversion_matrix = RNSMatrix([target], moduli_out, params, None, premult, postmult, io, self.inbase, self.indigits)

    def convert_to_rns(self, integer, reducer):
        # digits = word_reinterpret(integer, self.inbase, self.indigits)
        # Hacky: clearly no reason to convert to avx and back in c++
        return self.conversion_matrix.rns_reduce(AVXVector([integer], False), reducer)

class ConvertFromRNS:

    def __init__(self, premult, postmult, target, moduli_in, params, outbase, io=None, q_correct=False):
        self.outbase = outbase
        self.outdigits = word_len(2*target, outbase) # Conditional correction will remove final redundance and therefore requires 2p output

        self.conversion_matrix = RNSMatrix(moduli_in, [target], params, None, premult, postmult, io, params.inbits, params.indigits, self.outdigits)
        self.target = target
        self.q_correct = q_correct

    def convert_from_rns(self, residues):
        # Not quite certain how to get the merged barrett to pop out
        digits_hi, digits_lo = self.conversion_matrix.rns_reduce(residues, None, False)
        # For now let's do it dumb
        i = int_reconstruct(digits_lo.store(), self.outbase)
        i += (int_reconstruct(digits_hi.store(), self.outbase) << self.outbase)
        # Better way
        # vperm cyclic shift
        # extract 0 element as carry
        # input 0 element from previous carry
        # add hi and lo
        # shift 52
        # cyclic shift
        # add
        # Then do 52 -> 64 bit extraction
        # print(i.bit_length(), i.bit_length()-target.bit_length())
        return i % self.target


class BasicRNS:

    def __init__(self, target, params, io=None, max_bound=256*100, q_correct=False):
        modbits = params.modbits
        mulbits = params.mulbits
        reducer = params.reducer
        self.r = params.r
        moduli, _, _ = gen_extend((target * max_bound)**2, -1, 1, modbits)
        self.r1 = RNSMatrix(moduli, moduli, params, target, self.r, self.r, io)  
        self.r2 = None
        self.m1 = []
        self.m2 = moduli
        self.reducer = reducer(self.m2, params, io)
        self.q_correct = q_correct
        self.r = params.r
        self.factor2 = total_modulus(self.m2)
        self.convert_to = ConvertToRNS(1, self.r, target, self.m2, params, params.inbits, io)
        self.convert_from = ConvertFromRNS(1, 1, target, self.m2, params, params.outbits, io)

    def mulreduce(self, a_m1, a_m2, b_m1, b_m2):
        a = a_m2
        b = b_m2
        assert(len(a) == len(b))
        ab = ele_mult(a, b, self.reducer)
        return AVXVector([]), self.r1.rns_reduce(ab, self.reducer, self.q_correct)

    def to_mont_avx(self, a):
        c_m2 = self.convert_to.convert_to_rns(a, self.reducer)
        return AVXVector([]), c_m2
    
    def from_mont_avx(self, a_m2):
        return self.convert_from.convert_from_rns(a_m2)

class IntRNS2:

    bound = 3

    def __init__(self, target, params, io=None, max_bound=9):
        modbits = params.modbits
        reducer = params.reducer
        moduli1, state, m1 = gen_extend(target * max_bound, -1, 1, modbits)
        moduli2, _, _ = gen_extend(target * 6, state, m1, modbits)
        
        self.limbs1 = len(moduli1)
        self.limbs2 = len(moduli2)
        self.limbs = self.limbs1 + self.limbs2
        
        self.m1 = moduli1
        self.m2 = moduli2
        m1 = total_modulus(self.m1)
        m2 = total_modulus(self.m2)
        assert(m2 > 6 * target)
        self.mont_factor = get_mont_factor(target, m1)
        self.inv_factor = pow(m1, -1, m2)
        self.factor = m1
        self.factor2 = m2

        self.r = params.r
       
        self.r1 = RNSMatrix(moduli1, moduli2, params, None, self.mont_factor * pow(self.r, -1, m1), target * self.inv_factor * self.inv_factor * self.r * self.r, io)  
        self.r2 = RNSMatrix(moduli2, moduli1, params, None, self.factor * pow(self.r, -1, m2), self.r * self.r, io)
        self.reducer1 = reducer(self.m1, params, io) 
        self.reducer2 = reducer(self.m2, params, io)

        io2 = None # could write conversion to file as well
        self.convert_to = ConvertToRNS(self.factor, self.inv_factor * self.r * self.r, target, moduli2, params, params.inbits, io2)
        self.convert_from = ConvertFromRNS(self.factor * pow(self.r, -1, self.factor2), pow(self.factor, -1, target), target, self.m2, params, params.outbits, io2, False) # TODO: set to true and reduce big integer arithmetic requirement
    
    def to_mont_avx(self, a):
        # It would be slightly more efficient to use 1 matrix, but the slicing is really annoying with AVX
        # Because to have fast slicing I need to divide on a vector boundary which means there is padding
        # in the middle of the vector (rather than just at the end)
        c_m2 = self.convert_to.convert_to_rns(a, self.reducer2)
        c_m1 = self.r2.rns_reduce(c_m2, self.reducer1)
        return c_m1, c_m2

    def from_mont(self, a_m1, a_m2):
        a_m1 = [(a_m1[i] * pow(self.r, -1, m)) % m for i,m in enumerate(self.m1)]
        a_m2 = [(a_m2[i] * self.factor * pow(self.r, -1, m)) % m for i,m in enumerate(self.m2)]
        return a_m1, a_m2
    
    def from_mont_avx(self, a_m2):
        return self.convert_from.convert_from_rns(a_m2)

    def mulreduce(self, a_m1, a_m2, b_m1, b_m2):
        ab_m1 = ele_mult(a_m1, b_m1, self.reducer1)
        c_m2 = self.r1.rns_reduce(ab_m1, self.reducer2, False, (a_m2, b_m2))
        c_m1 = self.r2.rns_reduce(c_m2, self.reducer1)
        return c_m1, c_m2

def test_mulreduction(method, params, target, power=1, io=None, seeded=None):
    if seeded:
        seed(seeded)
    r = method(target, params, io)
    a = randint(0, target)
    b = randint(0, target)
   
    a_m1, a_m2 = r.to_mont_avx(a)
    b_m1, b_m2 = r.to_mont_avx(b)
   
    for _ in range(power):
        a_m1, a_m2 = r.mulreduce(a_m1, a_m2, b_m1, b_m2)

    c = r.from_mont_avx(a_m2)

    correct = ((a * pow(b, power, target)) % target)
    if seeded is not None:
        return c % target == correct
    assert(c % target == correct)
    print(f"pass {target.bit_length()} {method.__name__} {power}")

def test_a_lot():
    cpu_params_mont = Params.avx(MontgomeryReduce)
    for i in range(10):
        seed(i)
        modbits = randint(32, 4096)
        modulus = nextprime(randint(1, 2**modbits)) # probably needs to be prime or montgomery may fail?
        test_mulreduction(IntRNS2, cpu_params_mont, modulus, 10)

if __name__ == "__main__":
    target = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
    s = randint(0, 2**64)
    print("seed", s)
    seed(s)

    cpu_params_py = Params.avx(PythonReduce)
    cpu_params_mont = Params.avx(MontgomeryReduce)
    gpu_params_mont = Params.gput(MontgomeryReduce)

    gpu_params_30 = Params.gput(MontgomeryReduce)
    gpu_params_30.set30()

    gpu_params_vec = Params.gpuv(MontgomeryReduce)

    # gpu_params_mont30 = Params.gpu30(MontgomeryReduce)
    params = [cpu_params_py, cpu_params_mont, gpu_params_mont, gpu_params_30, gpu_params_vec]
    for p in params:
        p.refresh()
        test_mulreduction(BasicRNS, p, target, 10)
        test_mulreduction(IntRNS2, p, target, 10)
    test_a_lot()
    
    