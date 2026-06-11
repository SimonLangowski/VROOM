from math import floor, sqrt
from precompute_montgomery import to_mont
from rns_helpers import *

# From google: flatten array
def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]

def transpose_flatten(nested_list):
    length = len(nested_list[0])
    assert(all([len(sublist) == length for sublist in nested_list]))
    return flatten([[sublist[i] for sublist in nested_list] for i in range(length)])

def to_radix_rns(value, moduli, base, digits):
    residues = to_rns(value, moduli)
    return flatten([word_reinterpret(r, base, digits) for r in residues])

# Outputs a matrix that takes x in RNS in moduli_info
# The product with this matrix will be y = [[(x*premult)%M1] *tmult %target] * postmult in the RNS moduli2_info
def gen_precompute_wide(moduli_in, target_modulus, moduli_out=None, premult=1, postmult=1, indigits=1, inbits=8, outdigits=1, outbits=8, precision=64):
    if moduli_out == None:
        moduli_out = moduli_in
    modulus = total_modulus(moduli_in)
    modulus_out = total_modulus(moduli_out)
    if target_modulus == None:
        target_modulus = modulus_out
    icrt_factors = rns_precompute(moduli_in)
    inbase = 2**inbits
    shifted_icrt_factors = [[(premult * (inbase**j) * i) % modulus for i in icrt_factors] for j in range(indigits)]
    target_shifted_icrt_factors = [[((f % target_modulus) * postmult) % modulus_out for f in r] for r in shifted_icrt_factors]

    # Precision required to correct quotient
    precision_requirement = precision #max_k.bit_length() + 2
    fixed_point = 2**precision_requirement

    shifted_quotient_estimations = [[ceil_div(p * fixed_point, modulus) for p in f] for f in shifted_icrt_factors]
    target_shifted_quotient_estimations = [[ceil_div((p % target_modulus) * fixed_point, target_modulus) for p in f] for f in shifted_icrt_factors]

    rns_values = [[to_radix_rns(f, moduli_out, outbits, outdigits) for f in r] for r in target_shifted_icrt_factors]
    correction = to_radix_rns((((-modulus) % target_modulus) * postmult) % modulus_out, moduli_out, outbits, outdigits)
    tqc = ceil_div((((-modulus) % target_modulus)) * fixed_point, target_modulus)
    tp = to_radix_rns((-target_modulus * postmult) % modulus_out, moduli_out, outbits, outdigits)
    return transpose_flatten(rns_values), correction, transpose_flatten(shifted_quotient_estimations), transpose_flatten(target_shifted_quotient_estimations), tqc, tp, precision_requirement, moduli_in, moduli_out, (shifted_icrt_factors, target_shifted_icrt_factors, modulus, modulus_out, target_modulus, premult, postmult)