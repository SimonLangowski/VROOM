import numpy as np
from rns_helpers import ceil_div, int_reconstruct

# Convert python arrays to numpy matrices
def precompute_to_matrix_i(precompute):
    rns_values, correction, shifted_quotient_estimations, target_shifted_quotient_estimations, tqc, tp, precision_requirement, moduli, moduli_out, _ = precompute
    flatdim = len(rns_values)
    residuedim = len(rns_values[0])
    rns_mat = np.zeros((flatdim, residuedim), dtype=np.uint32)
    sqe_mat = np.zeros((flatdim, 1), dtype=np.uint64)
    tsqe_mat = np.zeros((flatdim, 1), dtype=np.uint64)
    for flat_idx in range(flatdim):
        for k in range(residuedim):
            rns_mat[flat_idx, k] = rns_values[flat_idx][k]
        sqe_mat[flat_idx, 0] = shifted_quotient_estimations[flat_idx]
        tsqe_mat[flat_idx, 0] = target_shifted_quotient_estimations[flat_idx]
    cor_mat = np.array(correction, dtype=np.uint32)
    tp_mat = np.array(tp, dtype=np.uint32)
    return rns_mat, cor_mat, sqe_mat, tsqe_mat, tqc, tp_mat, precision_requirement, moduli, moduli_out

# Split u32 into 4 u8
def get_parts(mat):
    u32mat = mat.astype(np.uint32)
    assert(u32mat.dtype == np.uint32)
    u32bytes = u32mat.view(np.uint8)
    u32bytesByByte = [u32bytes[:,::4], u32bytes[:,1::4], u32bytes[:,2::4], u32bytes[:,3::4]]
    return u32bytesByByte

# Will have to call on each part of b_pre (idk why the correction and matrix are still together)
def strided(mat, args, zero_fill=True):
    stride = args.stride
    rows, cols = mat.shape
    if zero_fill:
        # Space out array at striding and fill with zeros
        assert(cols == args.btpi) # true if mpt = 1
        out = np.zeros((rows, cols*stride), dtype=np.uint32)
        out[:,::stride] = mat[:,:]
        # for i in range(rows):
        #     for j in range(stride):
        #         out[i, j * stride] = mat[i, j]
        return out
    else:

        # Maybe copy the cpp function and then just match the columns here
        # Because cpp function will get more complicated with bank conflicts... (LDS.128 with no conflicts is ideal, but not always divisible...)
        # Restride matrix by splitting columns in new entries
        assert(cols == args.btpi)
        assert(cols * stride == args.tpi)
        out = np.zeros((ceil_div(rows, stride), cols*stride), dtype=np.uint32)
        i = 0
        per = ceil_div(rows, stride)
        for i in range(rows):
            nr = i % per
            nc = i // per
            out[nr, nc::stride] = mat[i, :]
        return out

def unstride(mat, args):
    stride = args.stride
    return mat[:,::stride]

def basic_pack(precompute_mats, args):
    rns_mat, cor_mat, sqe_mat, tsqe_mat, tqc, tp_mat, precision_requirement, moduli, moduli_out = precompute_mats
    flatdim, residuedim = rns_mat.shape
    # Decided to put sqe in first two threads
    sqe_threads = 2
    tpi = args.btpi
    mpt = ceil_div(residuedim + sqe_threads, tpi)
    # for now t will just be set to num moduli
    t = flatdim
    packed_b = np.zeros((tpi * mpt, t), dtype=np.uint32)
    packed_b[0, :] = sqe_mat[:, 0] >> 32
    packed_b[1, :] = sqe_mat[:, 0] & (2**32 - 1)
    for i in range(residuedim):
        packed_b[i + 2, :] = rns_mat[:, i]
    # First two entries are empty (compute fixed point only)
    # It could be one thread or two depending on mpt
    np_pad_cor_mat = np.zeros(tpi*mpt, dtype=np.uint32)
    for i in range(residuedim):
        np_pad_cor_mat[i + 2] = cor_mat[i]
    
    np_pad_cor_mat = strided(np_pad_cor_mat[:, np.newaxis].T, args)
    packed_b = strided(packed_b.T, args, False)
    # Stack
    return np.vstack((np_pad_cor_mat, packed_b))
    #stacked = np.hstack((np_pad_cor_mat[:, np.newaxis], packed_b))
    # Transpose?
    #return np.transpose(stacked).copy()

# Arrange the elements of the rns to align with the batching shape for IMMA
def flatten_a32(args, a_mat):
    m, k = a_mat.shape
    if args.variation == 4:
        return strided(a_mat, args)

# Invert flattena32 so we can get outputs in python
def unflatten_32(args, c_mat, shape):
    if args.variation == 4:
        return unstride(c_mat, args)
