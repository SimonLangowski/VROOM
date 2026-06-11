import itertools
from random import seed
from rns_helpers import *
from precompute_matrix import *
from precompute_montgomery import get_mont_factor
from debug import *
from sympy.ntheory import sqrt_mod

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

    # Do you want this for faster elementwise reduction?
    # modbits = mulbits - 2 should skip a round of reduction (at least in Montgomery, also in psuedomersenne?)
    # Or skipping the correction after addition for the first mult?
    @staticmethod
    def avx_redundant(reducer):
        return Params(modbits=50, mulbits=52, precision=64, inbits=52, indigits=1, reducer=reducer)

    @staticmethod
    def avx_very_redundant(reducer):
        return Params(modbits=44, mulbits=52, precision=52, inbits=44, indigits=1, reducer=reducer)

    @staticmethod
    def gput(reducer):
        return Params(modbits=32, mulbits=32, precision=32, inbits=8, indigits=4, reducer=reducer)
    
    @staticmethod
    def gpuv(reducer):
        return Params(modbits=32, mulbits=32, precision=64, inbits=32, indigits=1, reducer=reducer)
    
    def set30(self):
        self.modbits = 30
        self.refresh()

    @staticmethod
    def gput30(reducer):
        p = Params.gput(reducer)
        p.set30()
        return p

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

    def perm(self, other):
        if isinstance(other, AVXVector):
            indices = other.data
        else:
            indices = other
        for index in indices:
            if index < 0 or index >= self.length:
                raise ValueError(f"Index {index} out of bounds {indices}")
        return self.v([self.data[indices[i]] for i in range(len(indices))])

    def broadcast_lane(self, lane):
        return self.perm([lane] * self.length)

    def maskz_madd52lo(self, mask, b, c):
        mask52 = (1 << self.mulbits) - 1
        out = []
        for i in range(self.length):
            if (mask >> i) & 1:
                prod = ((b.data[i] & mask52) * (c.data[i] & mask52)) & mask52
                out.append(self.data[i] + prod)
            else:
                out.append(0)
        return self.v(out)

    def maskz_madd52hi(self, mask, b, c):
        mask52 = (1 << self.mulbits) - 1
        out = []
        for i in range(self.length):
            if (mask >> i) & 1:
                prod = ((b.data[i] & mask52) * (c.data[i] & mask52)) >> self.mulbits
                out.append(self.data[i] + prod)
            else:
                out.append(0)
        return self.v(out)

    def madd52lo(self, b, c):
        mask52 = (1 << self.mulbits) - 1
        return self.v([
            self.data[i] + ((b.data[i] & mask52) * (c.data[i] & mask52) & mask52)
            for i in range(self.length)
        ])

    def madd52hi(self, b, c):
        mask52 = (1 << self.mulbits) - 1
        return self.v([
            self.data[i] + (((b.data[i] & mask52) * (c.data[i] & mask52)) >> self.mulbits)
            for i in range(self.length)
        ])


def matvec_rows_scalar(matrix_rows, vec):
    """out[j] += sum_i matrix_rows[i][j] * vec[i]  (RNS-style row accumulation)."""
    n = len(vec)
    out = [0] * n
    for i, row in enumerate(matrix_rows):
        s = vec[i]
        for j in range(n):
            out[j] += row[j] * s
    return out


DOT_LANE = 0
AVX8_WIDTH = 8


def mod_positive(a, n):
    return (a % n + n) % n


def count_true_moduli(moduli):
    """Logical moduli count (exclude mod-1 padding entries)."""
    return sum(1 for m in moduli if m > 1)


def pad_perm_indices(indices, width=AVX8_WIDTH):
    """Pad gather indices to AVX width; extra lanes map to themselves."""
    p = list(indices)
    while len(p) < width:
        p.append(len(p))
    if len(p) > width:
        raise ValueError(f"perm indices length {len(p)} exceeds AVX width {width}")
    return p


def advance_perm_mt_c_avx8(n_moduli, width=AVX8_WIDTH):
    """advance_perm_mt_c on lanes 0..n_moduli, identity on padding lanes."""
    assert n_moduli > 0 and n_moduli < width
    return pad_perm_indices(advance_perm_mt_c(n_moduli), width)


def rns_matrix_perm_core(matrix_rows, n_moduli):
    """n_moduli×n_moduli active block: drop mod-1 input row/col and trailing zero rows."""
    rows = strip_mod1_row(list(matrix_rows))
    rows = strip_trailing_zero_rows(rows)
    assert len(rows) == n_moduli, f"expected {n_moduli} active rows, got {len(rows)}"
    ncol = len(rows[0])
    mod1_input_col = ncol > n_moduli and all(row[0] == 0 for row in rows)
    if mod1_input_col:
        return [row[1:n_moduli + 1] for row in rows]
    assert ncol == n_moduli, f"expected {n_moduli} cols (or mod-1 + padding), got {ncol}"
    return [row[0:n_moduli] for row in rows]


def matvec_mt_c_scalar(matrix_rows, c, vec):
    """M^T @ v in lanes 1..n and c·v in lane 0 (packed layout)."""
    n = len(vec)
    mt = [0] * n
    dot = 0
    for i in range(n):
        dot += c[i] * vec[i]
        for j in range(n):
            mt[j] += matrix_rows[i][j] * vec[i]
    return [dot] + mt


def matvec_mt_c_fixed(matrix_rows, c, vec):
    """6×7 fixed layout: skip mod-1 input row; vec is [v_mod1, v_m0, …, v_m5]."""
    nrows = len(matrix_rows)
    ncol = len(vec)
    assert ncol == nrows + 1, f"expected {nrows + 1} lanes, got {ncol}"
    assert all(len(row) == ncol for row in matrix_rows)
    dot = sum(c[j] * vec[j] for j in range(ncol))
    mt = [0] * ncol
    for i in range(nrows):
        s = vec[i + 1]
        for j in range(ncol):
            mt[j] += matrix_rows[i][j] * s
    return [dot] + mt[1:]


def pack_vec_with_dot_slot(vec):
    """Pack v as [0, v[0], ..., v[n-1]] so lane 0 can accumulate c·v."""
    return AVXVector([0] + list(vec))


def pack_vec_with_dot_slot_last(vec):
    """Pack v as [v[0], ..., v[n-1], 0] so the last lane accumulates c·v."""
    return AVXVector(list(vec) + [0])

def advance_perm_mt(n):
    # Simple cyclic shift
    return [((j - 1) % n) for j in range(n)]

def advance_perm_mt_c(n):
    """step_perms[0] for packed [0, v0, ..., v_{n-1}]: used for every cumulative advance."""
    return [1] + [((j - 1) % n) + 1 for j in range(n)]


def advance_perm_mt_c_last(n):
    """Cumulative advance for packed [v0, ..., v_{n-1}, 0] (dot in last lane)."""
    old = advance_perm_mt_c(n)
    def to_new(o):
        return n if o == 0 else o - 1
    return [to_new(old[l + 1]) for l in range(n)] + [to_new(old[0])]


def _lane_v_indices_after_advances(n, n_advances, advance_perm, dot_packed_index=0):
    """After n_advances steps, v-index at each lane (None = dot slot / zero)."""
    labels = list(range(n + 1))
    for _ in range(n_advances):
        labels = [labels[p] for p in advance_perm]
    return [None if x == dot_packed_index else x - 1 for x in labels]


def _lane_v_indices_first_dot(n, n_advances):
    return _lane_v_indices_after_advances(n, n_advances, advance_perm_mt_c(n), dot_packed_index=0)


def _lane_v_indices_last_dot(n, n_advances):
    return _lane_v_indices_after_advances(n, n_advances, advance_perm_mt_c_last(n), dot_packed_index=0)


def matvec_mt_c_scalar_last(matrix_rows, c, vec):
    """M^T @ v in lanes 0..n-1 and c·v in lane n (last column layout)."""
    n = len(vec)
    mt = [0] * n
    dot = 0
    for i in range(n):
        dot += c[i] * vec[i]
        for j in range(n):
            mt[j] += matrix_rows[i][j] * vec[i]
    return mt + [dot]


def pre_permute(matrix_rows):
    """Pre-shift M^T for diagonal matvec; returns coeffs and per-step gather indices.

    matrix_rows[i][j] is M[i][j].  step_perms[d][j] = (j-d)%n gathers from the original v.
    """
    n = len(matrix_rows)
    assert all(len(row) == n for row in matrix_rows)
    perm_mat = []
    step_perms = []
    for d in range(n):
        step_perms.append([(j - d) % n for j in range(n)])
        coeffs = [matrix_rows[(j - d) % n][j] for j in range(n)]
        perm_mat.append(AVXVector(coeffs))
    return perm_mat, step_perms


def pre_permute_fixed_mt_c(matrix_rows, c):
    """Fixed [k | residues] layout: 6 rows × 7 cols, no vector rotation.

    Row d uses input scalar vec[d+1]; lane 0 accumulates c[d+1]*vec[d+1] for k.
    matrix_rows[d] is the d-th row after dropping the all-zero mod-1 input row.
    c has one entry per input modulus (len = ncols = nrows + 1).
    """
    nrows = len(matrix_rows)
    ncols = len(matrix_rows[0])
    assert len(c) == ncols
    assert ncols == nrows + 1
    perm_mat = []
    for d in range(nrows):
        coeffs = list(matrix_rows[d])
        coeffs[DOT_LANE] = c[d + 1]
        perm_mat.append(AVXVector(coeffs))
    return perm_mat


def perm_mat_mul_fixed(perm_mat, vec):
    """Diagonal matvec with residues fixed in lanes 1..n; lane 0 holds k."""
    nlanes = len(vec.data)
    assert nlanes == len(perm_mat[0])
    out_hi = AVXVector.zero(nlanes)
    out_lo = AVXVector.zero(nlanes)
    for d in range(len(perm_mat)):
        s = vec.data[d + 1]
        out_hi = out_hi.mulhi_scalar(perm_mat[d], s)
        out_lo = out_lo.mullo_scalar(perm_mat[d], s)
    return out_hi, out_lo


def pre_permute_mt_c_avx8(matrix_rows, c, width=AVX8_WIDTH):
    """pre_permute_mt_c with perm rows and advance_perm padded to width AVX lanes."""
    perm_mat, advance_perm = pre_permute_mt_c(matrix_rows, c)
    perm_mat_avx = [AVXVector(pad_row_to_avx8(row.store())) for row in perm_mat]
    advance_perm_avx = pad_perm_indices(advance_perm, width)
    return perm_mat_avx, advance_perm_avx


def pre_permute_mt(matrix_rows):
    """Pre-shift M^T for cumulative diagonal matvec via advance_perm_mt.

    perm_mat_mul_cumulative(..., skip_first_advance=True) advances before steps
    d>=1 only.  After d advances, lane j holds v[(j-d)%n]; coeffs[d][j] = M[(j-d)%n][j].
    """
    advance_perm = advance_perm_mt(len(matrix_rows))
    perm_mat = []
    n = len(matrix_rows)
    for d in range(len(matrix_rows)):
        coeffs = [matrix_rows[(j - d) % n][j] for j in range(n)]
        perm_mat.append(AVXVector(coeffs))
    return perm_mat, advance_perm


def pre_permute_mt_avx8(matrix_rows, width=AVX8_WIDTH):
    """pre_permute_mt with perm rows and advance_perm padded to width AVX lanes."""
    perm_mat, advance_perm = pre_permute_mt(matrix_rows)
    perm_mat_avx = [AVXVector(pad_row_to_avx8(row.store())) for row in perm_mat]
    advance_perm_avx = pad_perm_indices(advance_perm, width)
    return perm_mat_avx, advance_perm_avx


def pre_permute_avx8(matrix_rows, width=AVX8_WIDTH):
    """pre_permute with per-step gather indices padded to width AVX lanes."""
    perm_mat, step_perms = pre_permute(matrix_rows)
    perm_mat_avx = [AVXVector(pad_row_to_avx8(row.store())) for row in perm_mat]
    step_perms_avx = [pad_perm_indices(sp, width) for sp in step_perms]
    return perm_mat_avx, step_perms_avx


def pre_permute_mt_c(matrix_rows, c):
    """Pre-shift [c | M^T] for cumulative rotation via advance_perm_mt_c.

    Online: apply advance_perm_mt_c before each step (not identity at d=0).
    Lane 0 cycles through all v[i]; c[d] column is rearranged to c[(-d)%n].
    Matrix coeffs use the (j-d-1) diagonal aligned to the rotated lanes.
    Returns (perm_mat, advance_perm).
    """
    n = len(matrix_rows)
    assert len(c) == n
    assert all(len(row) == n for row in matrix_rows)
    advance_perm = advance_perm_mt_c(n)
    perm_mat = []
    dot_lane = 0
    for d in range(n):
        lanes = _lane_v_indices_first_dot(n, d + 1)
        coeffs = [0] * (n + 1)
        coeffs[dot_lane] = c[(-d) % n]
        for l in range(1, n + 1):
            k = (l - d - 2) % n
            coeffs[l] = matrix_rows[k][l - 1]
        perm_mat.append(AVXVector(coeffs))
    return perm_mat, advance_perm


def pre_permute_mt_c_last(matrix_rows, c):
    """Like pre_permute_mt_c but dot product accumulates in the last lane.

    Packed [v0, ..., v_{n-1}, 0]; advance_perm_mt_c_last before each step.
    Lane n cycles v[(-d)%n] with c[(-d)%n]; lanes 0..n-1 carry the M^T diagonals.
    """
    n = len(matrix_rows)
    assert len(c) == n
    assert all(len(row) == n + 1 for row in matrix_rows)
    advance_perm = advance_perm_mt_c_last(n)
    perm_mat = []
    dot_lane = n
    for d in range(n):
        lanes = _lane_v_indices_last_dot(n, d + 1)
        coeffs = [0] * (n + 1)
        coeffs[dot_lane] = c[(-d) % n]
        for l in range(n):
            k = (l - d - 1) % n
            coeffs[l] = matrix_rows[k][l]
        perm_mat.append(AVXVector(coeffs))
    return perm_mat, advance_perm


def perm_mat_mul(perm_mat, vec, step_perms):
    """Diagonal matvec: gather from original v via step_perms[d], then mulhi/mullo."""
    nsteps = len(perm_mat)
    nlanes = len(vec.data)
    assert len(step_perms) == nsteps
    assert all(len(p) == nlanes for p in step_perms)
    assert all(len(perm_mat[d]) == nlanes for d in range(nsteps))
    out_hi = AVXVector.zero(nlanes)
    out_lo = AVXVector.zero(nlanes)
    for d in range(nsteps):
        v_step = vec.perm(step_perms[d])
        out_hi = out_hi.mulhi(perm_mat[d], v_step)
        out_lo = out_lo.mullo(perm_mat[d], v_step)
    return out_hi, out_lo


def perm_mat_mul_cumulative(perm_mat, vec, advance_perm, out_hi=None, out_lo=None, skip_first_advance=False):
    """Diagonal matvec: advance_perm applied to running vector before each multiply.

    skip_first_advance=True (pre_permute_mt / 8×8 r1): step 0 uses v as-is; no
    leading rotate.  mt_c paths pass skip_first_advance=False (advance at d=0).
    """
    nsteps = len(perm_mat)
    nlanes = len(vec.data)
    assert len(advance_perm) == nlanes
    if out_hi is None:
        out_hi = AVXVector.zero(nlanes)
        out_lo = AVXVector.zero(nlanes)
    assert len(out_hi) == nlanes
    assert len(out_lo) == nlanes
    v_cur = vec
    for d in range(nsteps):
        if not (skip_first_advance and d == 0):
            v_cur = v_cur.perm(advance_perm)
        assert len(v_cur) == len(perm_mat[d])
        out_hi = out_hi.mulhi(perm_mat[d], v_cur)
        out_lo = out_lo.mullo(perm_mat[d], v_cur)
    return out_hi, out_lo


def extract_k_broadcast_m128_lane0(out_hi, out_lo, mulbits=52):
    """Lane-0 k: hi[0] + (lo[0] >> mulbits), broadcast to all lanes."""
    k = out_hi.data[0] + (out_lo.data[0] >> mulbits)
    return AVXVector([k] * out_hi.length)


def add_correction_perm(out_hi, out_lo, correction, mulbits=52):
    """Full-width madd52 correction (no lane mask)."""
    corr = correction if isinstance(correction, AVXVector) else AVXVector(correction)
    k_broadcast = extract_k_broadcast_m128_lane0(out_hi, out_lo, mulbits)
    out_hi = out_hi.madd52hi(corr, k_broadcast)
    out_lo = out_lo.madd52lo(corr, k_broadcast)
    return out_hi, out_lo


def correction_madd_mask(n_lanes, dot_lane):
    """Bitmask for maskz_madd52: apply correction on all lanes except the dot/k lane."""
    return ((1 << n_lanes) - 1) & ~(1 << dot_lane)


def extract_k_broadcast(out_hi, out_lo, dot_lane, mulbits=52):
    """k from dot lane: hi + (lo >> mulbits); upper (64-mulbits) bits of lo complete the carry."""
    k_vec = out_hi.add(out_lo.srli(mulbits))
    return k_vec.broadcast_lane(dot_lane)


def add_correction_perm_masked(out_hi, out_lo, correction, mulbits, dot_lane=DOT_LANE):
    """Vector correction: maskz_madd52 on moduli lanes; dot lane zeroed by maskz."""
    n_lanes = len(out_hi)
    mask = correction_madd_mask(n_lanes, dot_lane)
    correction_padded = pad_correction_for_dot(correction, n_lanes, dot_lane)
    k_broadcast = extract_k_broadcast(out_hi, out_lo, dot_lane, mulbits)
    out_hi = out_hi.maskz_madd52hi(mask, correction_padded, k_broadcast)
    out_lo = out_lo.maskz_madd52lo(mask, correction_padded, k_broadcast)
    return out_hi, out_lo


def pad_correction_for_dot(correction, n_lanes, dot_lane=DOT_LANE):
    """Correction aligned to n_lanes; already n-wide or embed with 0 on dot lane."""
    data = correction.store() if isinstance(correction, AVXVector) else list(correction)
    if len(data) == n_lanes:
        return AVXVector(data)
    assert n_lanes == len(data) + 1
    if dot_lane == 0:
        return AVXVector([0] + list(data))
    if dot_lane == len(data):
        return AVXVector(list(data) + [0])
    raise ValueError(f"unsupported dot_lane={dot_lane} for n={len(data)}")


def slice_moduli_lanes(v, dot_lane, n_moduli):
    """Drop dot lane; return n_moduli residue lanes for reduction."""
    if dot_lane == 0:
        return AVXVector(v.data[1:n_moduli + 1])
    if dot_lane == n_moduli:
        return AVXVector(v.data[:n_moduli])
    raise ValueError(f"unsupported dot_lane={dot_lane}")


def wide_vec_to_ints(hi, lo, mulbits=None):
    """Reconstruct per-lane wide integers from hi/lo accumulators."""
    if mulbits is None:
        mulbits = hi.mulbits
    return [(hi.data[i] << mulbits) + lo.data[i] for i in range(hi.length)]


def matvec_rows_avx_scalar(matrix_rows, vec):
    """Reference AVX path: broadcast each scalar across a row (mulhi_scalar)."""
    n = len(vec.data)
    out_hi = AVXVector.zero(n)
    out_lo = AVXVector.zero(n)
    for i, row in enumerate(matrix_rows):
        row_v = AVXVector(row)
        out_hi = out_hi.mulhi_scalar(row_v, vec.data[i])
        out_lo = out_lo.mullo_scalar(row_v, vec.data[i])
    return out_hi, out_lo


def test_diagonal_matmul(n=6, trials=200, seed_value=0):
    """Compare diagonal perm_mat_mul against scalar and broadcast AVX references."""
    from random import seed as rand_seed, randint
    rand_seed(seed_value)
    mulbits = AVXVector.mulbits if hasattr(AVXVector, "mulbits") else 52
    max_val = (1 << 20)  # keep products small enough for exact integer check
    for t in range(trials):
        matrix_rows = [[randint(0, max_val) for _ in range(n)] for _ in range(n)]
        vec_vals = [randint(0, max_val) for _ in range(n)]
        vec = AVXVector(vec_vals)
        expect = matvec_rows_scalar(matrix_rows, vec_vals)
        perm_mat, step_perms = pre_permute(matrix_rows)
        hi_d, lo_d = perm_mat_mul(perm_mat, vec, step_perms)
        got_d = wide_vec_to_ints(hi_d, lo_d, mulbits)
        hi_s, lo_s = matvec_rows_avx_scalar(matrix_rows, vec)
        got_s = wide_vec_to_ints(hi_s, lo_s, mulbits)
        if got_d != expect or got_s != expect:
            raise AssertionError(
                f"trial {t}: expect={expect}\n  diagonal={got_d}\n  scalar_avx={got_s}"
            )
    print(f"diagonal matmul: {trials} trials x {n} lanes passed (mulbits={mulbits})")


def test_diagonal_matmul_mt_c_avx8(n=6, width=AVX8_WIDTH, trials=200, seed_value=0):
    """M^T+c·v in lane 0, full width AVX8 simulation."""
    from random import seed as rand_seed, randint
    rand_seed(seed_value)
    mulbits = AVXVector.mulbits if hasattr(AVXVector, "mulbits") else 52
    max_val = (1 << 20)
    for t in range(trials):
        matrix_rows = [[randint(0, max_val) for _ in range(n)] for _ in range(n)]
        c = [randint(0, max_val) for _ in range(n)]
        vec_vals = [randint(0, max_val) for _ in range(n)]
        expect_core = matvec_mt_c_scalar(matrix_rows, c, vec_vals)
        expect = pad_row_to_avx8(expect_core)
        perm_mat, advance_perm = pre_permute_mt_c_avx8(matrix_rows, c, width)
        packed = AVXVector(pad_row_to_avx8([0] + vec_vals))
        hi_d, lo_d = perm_mat_mul_cumulative(perm_mat, packed, advance_perm)
        got = wide_vec_to_ints(hi_d, lo_d, mulbits)
        if got != expect:
            raise AssertionError(
                f"trial {t}: expect={expect}\n  diagonal={got}"
            )
    print(f"diagonal M^T+c·v AVX{width}: {trials} trials x {n} moduli passed (mulbits={mulbits})")


def test_diagonal_matmul_mt_c(n=6, trials=200, seed_value=0):
    """M^T@v plus c·v in lane 0 via duplicated gather indices."""
    from random import seed as rand_seed, randint
    rand_seed(seed_value)
    mulbits = AVXVector.mulbits if hasattr(AVXVector, "mulbits") else 52
    max_val = (1 << 20)
    for t in range(trials):
        matrix_rows = [[randint(0, max_val) for _ in range(n)] for _ in range(n)]
        c = [randint(0, max_val) for _ in range(n)]
        vec_vals = [randint(0, max_val) for _ in range(n)]
        expect = matvec_mt_c_scalar(matrix_rows, c, vec_vals)
        perm_mat, advance_perm = pre_permute_mt_c(matrix_rows, c)
        packed = pack_vec_with_dot_slot(vec_vals)
        hi_d, lo_d = perm_mat_mul_cumulative(perm_mat, packed, advance_perm)
        got = wide_vec_to_ints(hi_d, lo_d, mulbits)
        if got != expect:
            raise AssertionError(
                f"trial {t}: expect={expect}\n  diagonal={got}"
            )
    print(f"diagonal M^T+c·v: {trials} trials x {n}+1 lanes passed (mulbits={mulbits})")


def test_diagonal_matmul_mt_c_last(n=6, trials=200, seed_value=0):
    """M^T@v in lanes 0..n-1, c·v in last lane (cumulative rotation)."""
    from random import seed as rand_seed, randint
    rand_seed(seed_value)
    mulbits = AVXVector.mulbits if hasattr(AVXVector, "mulbits") else 52
    max_val = (1 << 20)
    for t in range(trials):
        matrix_rows = [[randint(0, max_val) for _ in range(n)] for _ in range(n)]
        c = [randint(0, max_val) for _ in range(n)]
        vec_vals = [randint(0, max_val) for _ in range(n)]
        expect = matvec_mt_c_scalar_last(matrix_rows, c, vec_vals)
        perm_mat, advance_perm = pre_permute_mt_c_last(matrix_rows, c)
        packed = pack_vec_with_dot_slot_last(vec_vals)
        hi_d, lo_d = perm_mat_mul_cumulative(perm_mat, packed, advance_perm)
        got = wide_vec_to_ints(hi_d, lo_d, mulbits)
        if got != expect:
            raise AssertionError(
                f"trial {t}: expect={expect}\n  diagonal={got}"
            )
    print(f"diagonal M^T+c·v (last dot): {trials} trials x {n}+1 lanes passed (mulbits={mulbits})")


def test_diagonal_matmul_mt(n=8, width=AVX8_WIDTH, trials=200, seed_value=0):
    """M^T @ v via pre_permute_mt + cumulative advance_perm_mt (8×8 BLS12-381 r1 path)."""
    from random import seed as rand_seed, randint
    rand_seed(seed_value)
    mulbits = AVXVector.mulbits if hasattr(AVXVector, "mulbits") else 52
    max_val = (1 << 20)
    for t in range(trials):
        matrix_rows = [[randint(0, max_val) for _ in range(n)] for _ in range(n)]
        vec_vals = [randint(0, max_val) for _ in range(n)]
        expect = matvec_rows_scalar(matrix_rows, vec_vals)
        if width > n:
            expect = pad_row_to_avx8(expect)
        perm_mat, advance_perm = pre_permute_mt_avx8(matrix_rows, width)
        vec = AVXVector(pad_row_to_avx8(vec_vals) if width > n else vec_vals)
        hi_d, lo_d = perm_mat_mul_cumulative(
            perm_mat, vec, advance_perm, skip_first_advance=True)
        got = wide_vec_to_ints(hi_d, lo_d, mulbits)
        if got != expect:
            raise AssertionError(
                f"trial {t}: expect={expect}\n  diagonal_mt={got}"
            )
    print(f"diagonal M^T (pre_permute_mt): {trials} trials x {n} moduli AVX{width} passed (mulbits={mulbits})")


def random_residues(moduli, seed_lane0_zero=True):
    """Full-width residue vector: [0, m0..m_{n-1}, 0] for padded mod-1 layout."""
    from random import randint
    lanes = []
    for i, m in enumerate(moduli):
        if m == 1:
            lanes.append(0)
        elif i == 0 and seed_lane0_zero:
            lanes.append(0)
        else:
            lanes.append(randint(0, m - 1))
    return lanes


def perm_reduced_residue_lanes(vec, n_moduli):
    """Productive moduli lanes 1..n_moduli (exclude k lane 0 and trailing mod-1 pad)."""
    return vec.store()[1:n_moduli + 1]


def test_rns_reduce_perm(intrns, trials=200, seed_value=42):
    """rns_reduce_perm (8-wide advance_perm_mt_c) must match scalar on residue lanes."""
    from random import seed as rand_seed
    rand_seed(seed_value)
    r2 = intrns.r2
    if r2.perm_mat is None:
        raise AssertionError("perm tables not built")
    assert r2.perm_n_lanes == AVX8_WIDTH
    assert len(r2.perm_mat) == r2.perm_n_moduli
    assert len(r2.advance_perm) == AVX8_WIDTH
    assert len(r2.rns_matrix) == len(intrns.m2)
    n = r2.perm_n_moduli
    reducer = intrns.reducer1

    def check_residues(ref, got, msg):
        if perm_reduced_residue_lanes(ref, n) != perm_reduced_residue_lanes(got, n):
            raise AssertionError(f"{msg}: ref={ref.store()} got={got.store()}")

    acc_a = AVXVector(random_residues(intrns.m2))
    acc_b = AVXVector(random_residues(intrns.m2))
    for _ in range(trials):
        rv = AVXVector(random_residues(intrns.m2))
        ref = r2._rns_reduce_scalar(rv, reducer)
        got = r2.rns_reduce_perm(rv, reducer)
        check_residues(ref, got, "perm mismatch")
        ref_acc = r2._rns_reduce_scalar(rv, reducer, acc=(acc_a, acc_b))
        got_acc = r2.rns_reduce_perm(rv, reducer, acc=(acc_a, acc_b))
        check_residues(ref_acc, got_acc, "perm+acc mismatch")
        ref_dispatch = r2.rns_reduce(rv, reducer)
        check_residues(ref, ref_dispatch, "rns_reduce dispatch mismatch scalar reference")
        ref_acc_dispatch = r2.rns_reduce(rv, reducer, acc=(acc_a, acc_b))
        check_residues(ref_acc, ref_acc_dispatch, "rns_reduce dispatch mismatch scalar reference with acc")
    print(f"rns_reduce_perm: {trials} trials passed (8-wide advance_perm_mt_c, dispatch)")


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
        self.mont_factor = AVXVector([pow(m, -1, self.r) for m in moduli])

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
        # print(f"  [DEBUG Python reduce_small] input hi={hi.store()}, low={low.store()}")
        q = AVXVector.zero(self.length).mullo(low, self.mont_factor)
        # print(f"  [DEBUG Python reduce_small] q={q.store()}, mont_factor={self.mont_factor.store()}")
        h = AVXVector.zero(self.length).mulhi(q, self.modulia)
        # print(f"  [DEBUG Python reduce_small] h={h.store()}, modulia={self.modulia.store()}")
        uh = hi.sub(h)
        # print(f"  [DEBUG Python reduce_small] uh={uh.store()}")
        if self.mulbits == self.bits:
            mask = h.mask_cmp_gt(hi)
            # print(f"  [DEBUG Python reduce_small] mask={mask.store()}, mulbits={self.mulbits}, bits={self.bits}")
            result = uh.mask_add(mask, self.modulia)
            # print(f"  [DEBUG Python reduce_small] final result={result.store()}")
            return result
        else:
            result = uh.add(self.modulia)
            # print(f"  [DEBUG Python reduce_small] final result (else)={result.store()}")
            return result

    # # need hi < m, lo < 2^52
    def reduce_full_old(self, hi, low):
        hi_hi = hi.srli(self.hi_bit_shift)
        hi_rem = hi & self.hi_mask
        # print(f"  [DEBUG Python reduce_full] hi_rem={hi_rem.store()}, hi_mask={self.hi_mask.store()}")
        # require hi_hi * t_squared < 2^52
        # required low < 2^64 - 2^52
        lo_acc = low.mullo(hi_hi, self.t_squared)
        # print(f"  [DEBUG Python reduce_full] lo_acc={lo_acc.store()}, t_squared={self.t_squared.store()}")

        lo_hi = lo_acc.srli(self.mulbits)
        lo_rem = lo_acc & self.lo_mask
        # print(f"  [DEBUG Python reduce_full] lo_hi={lo_hi.store()}, lo_rem={lo_rem.store()}, lo_mask={self.lo_mask.store()}")

        # require hi_rem < 2(2^52 - t)
        hi_acc = hi_rem.add(lo_hi)
        # print(f"  [DEBUG Python reduce_full] hi_acc={hi_acc.store()}")
        mask = hi_acc.mask_cmp_ge(self.modulia)
        hi_acc.mask_sub(mask, self.modulia)
        # print(f"  [DEBUG Python reduce_full] hi_acc after mask_sub={hi_acc.store()}")

        # Certifies that hi < 2^52 - t so quotient will fit
        result = self.reduce_small(hi_acc, lo_rem)
        # print(f"  [DEBUG Python reduce_full] final result={result.store()}")
        return result
    
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

def skip_mod1_input_row(matrix_rows):
    """Drop row 0 (mod-1 input); residues stay in lanes 1..n."""
    if matrix_rows and all(x == 0 for x in matrix_rows[0]):
        return matrix_rows[1:]
    return matrix_rows


class RNSMatrix:
    
    def __init__(self, moduli_in, moduli_out, params, target_modulus=None, premult=1, postmult=1, permute=False, skip_correction=False, inbits=None, indigits=None, outdigits=1, debug=False):
        if inbits == None:
            inbits = params.inbits
        if indigits == None:
            indigits = params.indigits
        self.premult = premult
        self.postmult = postmult
        precompute = gen_precompute_wide(moduli_in, target_modulus, moduli_out, premult, postmult, indigits, inbits, outdigits, params.mulbits, params.precision, debug=debug)
        self.precompute = precompute # For extraction via GPU stuff
        rns_values, correction, shifted_quotient_estimations, _, _, _, precision_requirement, _, _, _ = precompute
        self.rns_matrix = [AVXVector(r) for r in rns_values]
        self.skip_correction = skip_correction
        self.correction = AVXVector(correction)
        self.shifted_quotient_estimations = shifted_quotient_estimations
        self.precision = precision_requirement
        self.moduli_out = moduli_out
        self.inbits = inbits
        self.indigits = indigits
        self.mulbits = params.mulbits
        self.perm_n_lanes = len(self.correction)
        self.perm_n_moduli = count_true_moduli(moduli_in)
        if permute:
            wide_rows = [row.store() for row in self.rns_matrix]
            core = rns_matrix_perm_core(wide_rows, self.perm_n_moduli)
            if self.perm_n_moduli < 8:
                sqe = self.shifted_quotient_estimations[1:self.perm_n_moduli + 1]
                self.perm_mat, self.advance_perm = pre_permute_mt_c_avx8(core, sqe, self.perm_n_lanes)
            elif self.perm_n_moduli == 8:
                # 8×8 cyclic advance_perm_mt; no k column; skip leading advance at d=0
                self.perm_mat, self.advance_perm = pre_permute_mt_avx8(core, self.perm_n_lanes)
                self.perm_skip_first_advance = True
            else:
                raise ValueError(f"Unsupported number of moduli: {self.perm_n_moduli}")
            assert self.perm_n_lanes == AVX8_WIDTH
            if self.perm_n_moduli < 8:
                self.perm_skip_first_advance = False
        else:
            self.perm_mat = None
            self.advance_perm = None
            self.perm_skip_first_advance = False

    def _perm_accumulate_acc(self, out_hi, out_lo, acc):
        if acc is None:
            return out_hi, out_lo
        a, b = acc
        assert len(a) == self.perm_n_lanes and len(b) == self.perm_n_lanes
        ah = AVXVector.zero(self.perm_n_lanes).mulhi(a, b)
        al = AVXVector.zero(self.perm_n_lanes).mullo(a, b)
        return out_hi.add(ah), out_lo.add(al)

    def rns_reduce(self, residues1, reducer, acc=None):
        """Reduce via fixed-layout perm when tables exist; else scalar loop."""
        if self.perm_mat is not None:
            return self.rns_reduce_perm(residues1, reducer, acc=acc, skip_correction=self.skip_correction)
        return self._rns_reduce_scalar(residues1, reducer, acc)

    def _rns_reduce_scalar(self, residues1, reducer, acc=None):
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
        # debug_printvec("scalars", scalars)
        debug_printvec("sqe", self.shifted_quotient_estimations)
        sqe_hi = [x >> self.mulbits for x in self.shifted_quotient_estimations]
        debug_printvec("sqe hi", sqe_hi)
        sqe_lo = [x & mask(self.mulbits) for x in self.shifted_quotient_estimations]
        debug_printvec("sqe lo", sqe_lo)
        k_raw_lo = 0
        k_raw_hi = 0
        for (i, s) in enumerate(scalars):
            debug_printvec("bpre", self.rns_matrix[i])
            out_hi = out_hi.mulhi_scalar(self.rns_matrix[i], s)
            out_lo = out_lo.mullo_scalar(self.rns_matrix[i], s)
            k_raw += s * self.shifted_quotient_estimations[i]
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
        # print(f"  [DEBUG Python rns_reduce] k_raw={k_raw}, k={k}")
        # print(f"  [DEBUG Python rns_reduce] out_hi before correction={out_hi.store()}")
        # print(f"  [DEBUG Python rns_reduce] out_lo before correction={out_lo.store()}")
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
        # print(f"  [DEBUG Python rns_reduce] out_hi after correction={out_hi.store()}")
        # print(f"  [DEBUG Python rns_reduce] out_lo after correction={out_lo.store()}")
        if reducer == None:
            return out_hi, out_lo
        reduced = reducer.reduce_full(out_hi, out_lo)
        # print(f"  [DEBUG Python rns_reduce] final reduced={reduced.store()}")
        return reduced

    def rns_reduce_perm(self, residues1, reducer, acc=None, skip_correction=False):
        if self.perm_mat is None or len(residues1) != self.perm_n_lanes:
            return self._rns_reduce_scalar(residues1, reducer, acc)
        out_hi = AVXVector.zero(self.perm_n_lanes)
        out_lo = AVXVector.zero(self.perm_n_lanes)
        out_hi, out_lo = self._perm_accumulate_acc(out_hi, out_lo, acc)
        vec = residues1 if isinstance(residues1, AVXVector) else AVXVector(residues1)
        out_hi, out_lo = perm_mat_mul_cumulative(
            self.perm_mat, vec, self.advance_perm, out_hi, out_lo,
            skip_first_advance=self.perm_skip_first_advance)
        if not skip_correction:
            out_hi, out_lo = add_correction_perm(
                out_hi, out_lo, self.correction, self.mulbits)
        if reducer is None:
            return out_hi, out_lo
        return reducer.reduce_full(out_hi, out_lo)
    
class ConvertToRNS:

    def __init__(self, premult, postmult, target, moduli_out, params, inbase, io=None):
        self.inbase = inbase
        self.indigits = word_len(target, inbase)
        self.premult = premult
        self.postmult = postmult
        self.target = target
        self.moduli_out = moduli_out
        self.conversion_matrix = RNSMatrix(
            [target], moduli_out, params,
            target_modulus=None,
            premult=premult,
            postmult=postmult,
            inbits=self.inbase,
            indigits=self.indigits,
        )

    def convert_to_rns(self, integer, reducer):
        # digits = word_reinterpret(integer, self.inbase, self.indigits)
        # Hacky: clearly no reason to convert to avx and back in c++
        return self.conversion_matrix.rns_reduce(AVXVector([integer], False), reducer)

class ConvertFromRNS:

    def __init__(self, premult, postmult, target, moduli_in, params, outbase, io=None):
        self.outbase = outbase
        self.outdigits = word_len(2*target, outbase)

        self.conversion_matrix = RNSMatrix(
            moduli_in, [target], params,
            target_modulus=None,
            premult=premult,
            postmult=postmult,
            inbits=params.inbits,
            indigits=params.indigits,
            outdigits=self.outdigits,
        )
        self.target = target

    def convert_from_rns(self, residues):
        # Not quite certain how to get the merged barrett to pop out
        digits_hi, digits_lo = self.conversion_matrix.rns_reduce(residues, None)
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

def mod_sqrt(a, moduli):
    square_roots = []
    for m in moduli:
        s = sqrt_mod(a % m, m)
        if s is None:
            return None
        square_roots.append(s)
    return rns_reconstruct_full(square_roots, moduli)

def compute_small_factor(moduli):
    M = total_modulus(moduli)
    icrt_factors = [M//m for m in moduli]
    return rns_reconstruct_full(icrt_factors, moduli)

def generate_moduli_candidates(num_moduli, modulus_bits):
    modulus = 1
    mod_cand = 2**modulus_bits - 1
    moduli = []
    while len(moduli) < num_moduli:
        if gcd(modulus, mod_cand) == 1:
            moduli.append(mod_cand)
            modulus *= mod_cand
        mod_cand -= 2 # only odd numbers are coprime to Montgomery
    return moduli

SECP256K1_MODULUS = 115792089237316195423570985008687907853269984665640564039457584007908834671663
BN254_MODULUS = 21888242871839275222246405745257275088696311157297823662689037894645226208583
BLS12_381_MODULUS = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab


class IntRNS4:

    def __init__(self, target, params, num_moduli=6, bound=9):

        # multiply M by sqrt(ICRT * current pre mult factor)
        # Set premult factor to ICRT^-1
        # Find square root here, or let be provided?

        num_cands = 4 * num_moduli
        candidates = generate_moduli_candidates(num_cands, params.modbits)
        self.r = params.r
        # Iterate over candidate sets
        # Check if the number has a square root
        
        for combo in itertools.combinations(range(num_cands), num_moduli):
            if combo == (0, 1, 2, 3, 4, 5):
                if target == SECP256K1_MODULUS:
                    combo = (4, 5, 7, 11, 14, 17)  # secp256k1 fast path
                elif target == BN254_MODULUS:
                    combo = (1, 8, 9, 14, 15, 16)  # bn254 fast path
            if combo == (0, 1, 2, 3, 4, 5, 6, 7):
                combo = (3, 5, 9, 15, 16, 17, 19, 24)  # bls12-381 fast path
            moduli1 = [candidates[i] for i in combo]
            m1 = total_modulus(moduli1)
            self.mont_factor = get_mont_factor(target, m1)
            r1_premult = self.mont_factor * pow(self.r, -1, m1)
            m1_opt = compute_small_factor(moduli1)
            qr_premult = mod_sqrt(r1_premult * pow(m1_opt, -1, m1), moduli1)
            if qr_premult is not None:
                print(combo)
                break
        if qr_premult is None:
            raise ValueError("No square root found")
        self.qr_premult = qr_premult
        self.m1_opt = m1_opt
        
        j = 0
        moduli2 = []
        while len(moduli2) < num_moduli:
            if j not in combo:
                moduli2.append(candidates[j])
            j += 1

        # perpend 1 to each moduli set
        if num_moduli < 8:
            moduli1 = [1] + moduli1
            moduli2 = [1] + moduli2
        while len(moduli1) % 8 != 0:
            moduli1.append(1)
        while len(moduli2) % 8 != 0:
            moduli2.append(1)
        
        self.limbs1 = len(moduli1)
        self.limbs2 = len(moduli2)
        self.limbs = self.limbs1 + self.limbs2
        
        self.m1 = moduli1
        self.m2 = moduli2
        
        m2 = total_modulus(self.m2)
        assert(m2 > 6 * target)
        
        self.inv_factor = pow(m1, -1, m2)
        self.factor = m1
        self.factor2 = m2
        self.target = target

        
       
        print("Creating r1 matrix (Python):")
       
        r1_postmult = target * self.inv_factor * self.inv_factor * self.r * self.r
        r2_premult = self.factor * pow(self.r, -1, m2)
        r2_postmult = self.r * self.r
        convert_to_premult = self.factor
        convert_to_postmult = self.inv_factor * self.r * self.r
        convert_from_premult = self.factor * pow(self.r, -1, self.factor2)
        convert_from_postmult = pow(self.factor, -1, target)

        rot_factor = 1
        # Apply the premult sqrt on r2
        r1_premult = m1_opt
        r2_postmult *= qr_premult
        rot_factor *= qr_premult

        reducer = params.reducer
        # Base matrix matches IntRNS2 (target_modulus=None); perm_mat is the diagonal
        # M^T layout with shifted_quotient_estimations in lane 0 for k accumulation.
        self.r1 = RNSMatrix(
            moduli1, moduli2, params,
            target_modulus=None,
            premult=r1_premult,
            postmult=r1_postmult,
            permute=True,
            skip_correction=True,
        )
        # 8-moduli (BLS12-381): r1 uses 8×8 perm without k column; r2 stays scalar.
        r2_permute = num_moduli < 8
        self.r2 = RNSMatrix(
            moduli2, moduli1, params,
            target_modulus=None,
            premult=r2_premult,
            postmult=r2_postmult,
            permute=r2_permute,
            skip_correction=False,
        )
        self.reducer1 = reducer(self.m1, params)
        self.reducer2 = reducer(self.m2, params)
        self.convert_to = ConvertToRNS(convert_to_premult, convert_to_postmult, target, moduli2, params, params.modbits)
        self.convert_from = ConvertFromRNS(convert_from_premult, convert_from_postmult, target, self.m2, params, params.outbits)

        self.total_modulus = self.factor * self.factor2
        icrt1 = self.factor2 * pow(self.factor2, -1, self.factor)
        icrt2 = self.factor * pow(self.factor, -1, self.factor2)
        rotation_factor = ((rot_factor * icrt1 + self.inv_factor * icrt2)) % self.total_modulus
        self.total_rotation = (self.reducer2.r * rotation_factor) % self.total_modulus
        self.total_inv_rotation = pow(self.total_rotation, -1, self.total_modulus)

        # print(self.r1.shifted_quotient_estimations)

    def to_mont_avx(self, a):
        # It would be slightly more efficient to use 1 matrix, but the slicing is really annoying with AVX
        # Because to have fast slicing I need to divide on a vector boundary which means there is padding
        # in the middle of the vector (rather than just at the end)
        c_m2 = self.convert_to.convert_to_rns(a, self.reducer2)
        c_m1 = self.r2.rns_reduce(c_m2, self.reducer1)
        return c_m1, c_m2

    def to_mont_exact(self, a, wide=False):
        quo_a = a // self.target
        mont_a = (a * self.factor) % self.target
        if wide:
            mont_a = (mont_a * self.factor) % self.target
        mont_a += quo_a * self.target
        mont_a_r = (mont_a * self.total_rotation) % self.total_modulus
        if wide:
            mont_a_r = (mont_a * self.total_rotation) % self.total_modulus
        c_m1 = [mont_a_r % m for m in self.m1]
        c_m2 = [mont_a_r % m for m in self.m2]
        return AVXVector(c_m1), AVXVector(c_m2)

    def to_mont_py(self, a):
        return self.to_mont_exact(a)

    def from_mont_py(self, a):
        return self.from_mont_exact(a)

    def from_mont_exact(self, a, wide=False):
        a_m1, a_m2 = a
        mont_a_r = rns_reconstruct_full(concat(a_m1.store(), a_m2.store()), concat(self.m1, self.m2))
        mont_a = (mont_a_r * self.total_inv_rotation) % self.total_modulus
        if wide:
            mont_a = (mont_a * self.total_inv_rotation) % self.total_modulus
        # if mont_a > self.total_modulus // 2:
        #     mont_a -= self.total_modulus
        value = (mont_a * pow(self.factor, -1, self.target)) % self.target
        if wide:
            value = (value * pow(self.factor, -1, self.target)) % self.target
        # mont_a is the integer encoded in the RNS representation, which is used to understand the redundancy
        redundance = mont_a // self.target
        assert(redundance * self.target + mont_a % self.target == mont_a)
        return value, redundance, mont_a

    def from_mont_avx(self, a_m2):
        return self.convert_from.convert_from_rns(a_m2)

    def from_mont_avx2(self, a):
        a_m1, a_m2 = a
        return self.from_mont_avx(a_m2)

    def mulreduce(self, a_m1, a_m2, b_m1, b_m2):
        ab_m1 = ele_mult(a_m1, b_m1, self.reducer1)
        c_m2 = self.r1.rns_reduce(ab_m1, self.reducer2, acc=(a_m2, b_m2))
        c_m1 = self.r2.rns_reduce(c_m2, self.reducer1)
        return c_m1, c_m2

def test_mulreduction(method, params, target, power=1, io=None, seeded=None, num_moduli=6):
    if seeded:
        seed(seeded)
    if method is IntRNS4:
        r = method(target, params, num_moduli)
    else:
        r = method(target, params, io)
    a = randint(0, target)
    b = randint(0, target)

    assert(r.from_mont_avx(r.to_mont_py(a)[1])==a)
    assert(r.from_mont_avx(r.to_mont_py(b)[1])==b)
   
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

def format_array(arr):
    """Format a Python list as a C++ std::array initializer."""
    return "{" + ", ".join(str(x) for x in arr) + "}"


def format_matrix_rows(rows):
    """Format rows as C++ nested std::array initializer {{ {...}, {...} }}."""
    return "{{" + ", ".join(format_array(row) for row in rows) + "}}"

def strip_mod1_row(matrix_rows):
    """Drop row 0 when it is the all-zero mod-1 input row."""
    if matrix_rows and all(x == 0 for x in matrix_rows[0]):
        return matrix_rows[1:]
    return matrix_rows


def strip_trailing_zero_rows(matrix_rows):
    """Drop trailing all-zero rows (mod-1 output padding)."""
    rows = list(matrix_rows)
    while rows and all(x == 0 for x in rows[-1]):
        rows.pop()
    return rows


def pad_row_to_avx8(row):
    """Pad any row to 8 AVX lanes; keep column 0 and append zeros at the end."""
    data = list(row)
    while len(data) < 8:
        data.append(0)
    if len(data) > 8:
        raise ValueError(f"row has {len(data)} lanes, cannot pad to 8")
    return data


def perm_mat_rows_to_8(perm_mat):
    """Export perm_mat AVX rows as 6×8 (sqe already in lane 0)."""
    rows = [pad_row_to_avx8(row.store()) for row in perm_mat]
    return strip_trailing_zero_rows(rows)


def print_rns_matrix_avx8(f, rns_matrix, prefix, skip_zero_input_row=False):
    """Print matrix rows as N×8 AVX tables; optionally drop all-zero mod-1 input row."""
    rows = [vec.store() for vec in rns_matrix.rns_matrix]
    if skip_zero_input_row:
        rows = strip_mod1_row(rows)
        rows = strip_trailing_zero_rows(rows)
    rows_8 = [pad_row_to_avx8(row) for row in rows]
    f.write(f'constexpr std::array<std::array<uint64_t, 8>, {len(rows_8)}> {prefix}_rns_mat = {format_matrix_rows(rows_8)};\n')
    correction_8 = pad_row_to_avx8(rns_matrix.correction.store())
    f.write(f'constexpr std::array<uint64_t, 8> {prefix}_correction = {format_array(correction_8)};\n')
    sqe = rns_matrix.shifted_quotient_estimations
    f.write(f'constexpr std::array<uint64_t, {len(sqe)}> {prefix}_shifted_quotient_estimations = {format_array(sqe)};\n')


def print_rns_matrix(f, rns_matrix, prefix):
    """Print RNS matrix fields (rns_mat, correction, shifted_quotient_estimations) to file.
    
    Args:
        f: File handle to write to
        rns_matrix: RNSMatrix object
        prefix: String prefix for variable names (e.g., "r1")
    """
    # Convert rns_matrix (list of AVXVector) to 2D list
    rns_mat_2d = [vec.store() for vec in rns_matrix.rns_matrix]
    
    if len(rns_mat_2d) > 0:
        inner_size = len(rns_mat_2d[0])
        outer_size = len(rns_mat_2d)
        # Format 2D array - need 3 levels of braces for nested std::array
        inner_arrays = [format_array(row) for row in rns_mat_2d]
        rns_mat_str = "{" + ", ".join(inner_arrays) + "}"
        f.write(f'const std::array<std::array<uint64_t, {inner_size}>, {outer_size}> {prefix}_rns_mat = {{{rns_mat_str}}};\n')
    
    # Print correction (AVXVector)
    correction_list = rns_matrix.correction.store()
    f.write(f'const std::array<uint64_t, {len(correction_list)}> {prefix}_correction = {format_array(correction_list)};\n')
    
    # Print shifted_quotient_estimations (list)
    sqe_list = rns_matrix.shifted_quotient_estimations
    f.write(f'const std::array<uint64_t, {len(sqe_list)}> {prefix}_shifted_quotient_estimations = {format_array(sqe_list)};\n')

def secp256k1_short_weierstrass_rhs(x, p, b=7):
    """Right-hand side of Y² = X³ + b (mod p) for secp256k1 G1."""
    return (pow(x, 3, p) + b) % p


def secp256k1_affine_on_curve(x, y, p, b=7):
    """All EC test inputs/outputs must satisfy this (projective formulas assume curve points)."""
    return (y * y - secp256k1_short_weierstrass_rhs(x, p, b)) % p == 0


def bn254_short_weierstrass_rhs(x, p, b=3):
    """Right-hand side of Y² = X³ + b (mod p) for BN254 G1."""
    return (pow(x, 3, p) + b) % p


def bn254_affine_on_curve(x, y, p, b=3):
    """All EC test inputs/outputs must satisfy this (projective formulas assume curve points)."""
    return (y * y - bn254_short_weierstrass_rhs(x, p, b)) % p == 0


def print_bn254_ec_tests(path, num_cases=50, seed_value=42):
    """Generate G1 EC test vectors for BN254 (py_ecc bn128 reference, b=3)."""
    from py_ecc.bn128 import G1, add, double, multiply, field_modulus, curve_order, b as curve_b

    seed(seed_value)
    p = int(field_modulus)
    order = int(curve_order)
    if int(curve_b) != 3:
        raise ValueError(f"expected BN254 b=3, got b={curve_b}")
    curve_b_int = 3

    def rand_scalar():
        return randint(1, order - 1)

    def aff(pt):
        return int(pt[0]), int(pt[1])

    def check_xy(label, x, y):
        if not bn254_affine_on_curve(x, y, p, curve_b_int):
            raise ValueError(
                f"{label} off curve Y²=X³+{curve_b_int}: x={x}, y={y}, "
                f"residue={(y * y - bn254_short_weierstrass_rhs(x, p, curve_b_int)) % p}"
            )

    cases = []
    for case_idx in range(num_cases):
        k1, k2, k3 = rand_scalar(), rand_scalar(), rand_scalar()
        P = multiply(G1, k1)
        Q = multiply(G1, k2)
        p_x, p_y = aff(P)
        q_x, q_y = aff(Q)

        R_add = add(P, Q)
        add_x, add_y = aff(R_add)

        # py_ecc is affine-only; mixed-add (Z=1) matches full add at these coords.
        mix_x, mix_y = add_x, add_y

        R_dbl = double(P)
        dbl_x, dbl_y = aff(R_dbl)

        Pz = multiply(G1, k3)
        pz_x, pz_y = aff(Pz)
        Q2 = multiply(G1, rand_scalar())
        q2_x, q2_y = aff(Q2)
        R_mz = add(Pz, Q2)
        mz_x, mz_y = aff(R_mz)

        case = {
            "p_x": p_x, "p_y": p_y, "q_x": q_x, "q_y": q_y,
            "add_x": add_x, "add_y": add_y,
            "mix_x": mix_x, "mix_y": mix_y,
            "dbl_x": dbl_x, "dbl_y": dbl_y,
            "pz_x": pz_x, "pz_y": pz_y, "q2_x": q2_x, "q2_y": q2_y,
            "mz_x": mz_x, "mz_y": mz_y,
        }
        for key, (x, y) in [
            ("P", (p_x, p_y)),
            ("Q", (q_x, q_y)),
            ("add", (add_x, add_y)),
            ("mix", (mix_x, mix_y)),
            ("dbl", (dbl_x, dbl_y)),
            ("Pz", (pz_x, pz_y)),
            ("Q2", (q2_x, q2_y)),
            ("mixed_z", (mz_x, mz_y)),
        ]:
            check_xy(f"case {case_idx} {key}", x, y)
        cases.append(case)

    with open(path, "w") as f:
        f.write("#pragma once\n")
        f.write("// BN254 G1 EC vectors (py_ecc bn128); Y²=X³+3; all coords on-curve\n")
        f.write(f"constexpr int EC_NUM_CASES = {num_cases};\n\n")
        f.write("struct Bn254EcCase {\n")
        f.write("    const char* p_x; const char* p_y;\n")
        f.write("    const char* q_x; const char* q_y;\n")
        f.write("    const char* add_x; const char* add_y;\n")
        f.write("    const char* mix_x; const char* mix_y;\n")
        f.write("    const char* dbl_x; const char* dbl_y;\n")
        f.write("    const char* pz_x; const char* pz_y;\n")
        f.write("    const char* q2_x; const char* q2_y;\n")
        f.write("    const char* mz_x; const char* mz_y;\n")
        f.write("};\n\n")
        f.write("constexpr Bn254EcCase EC_CASES[] = {\n")
        for c in cases:
            f.write("    {")
            keys = ["p_x", "p_y", "q_x", "q_y", "add_x", "add_y", "mix_x", "mix_y",
                    "dbl_x", "dbl_y", "pz_x", "pz_y", "q2_x", "q2_y", "mz_x", "mz_y"]
            f.write(", ".join(f'"{c[k]}"' for k in keys))
            f.write("},\n")
        f.write("};\n")
    print(f"wrote {path} ({num_cases} cases, seed={seed_value})")


def print_ec_tests(path, num_cases=50, seed_value=42):
    """Generate G1 EC test vectors for secp256k1 (ecdsa reference, b=7)."""
    from ecdsa import SECP256k1
    from ecdsa.ellipticcurve import PointJacobi

    seed(seed_value)
    curve = SECP256k1.curve
    order = SECP256k1.order
    G = SECP256k1.generator
    p = curve.p()
    curve_b = curve.b()
    if curve_b != 7:
        raise ValueError(f"expected secp256k1 b=7, got b={curve_b}")

    def rand_scalar():
        return randint(1, order - 1)

    def aff(pt):
        a = pt.to_affine()
        return int(a.x()), int(a.y())

    def check_xy(label, x, y):
        if not secp256k1_affine_on_curve(x, y, p, curve_b):
            raise ValueError(
                f"{label} off curve Y²=X³+{curve_b}: x={x}, y={y}, "
                f"residue={(y * y - secp256k1_short_weierstrass_rhs(x, p, curve_b)) % p}"
            )

    cases = []
    for case_idx in range(num_cases):
        k1, k2, k3 = rand_scalar(), rand_scalar(), rand_scalar()
        P = k1 * G
        Q = k2 * G
        p_x, p_y = aff(P)
        q_x, q_y = aff(Q)

        R_add = P + Q
        add_x, add_y = aff(R_add)

        Q_aff = PointJacobi.from_affine(Q.to_affine(), generator=G)
        R_mixed = P + Q_aff
        mix_x, mix_y = aff(R_mixed)

        R_dbl = P.double()
        dbl_x, dbl_y = aff(R_dbl)

        Pz = k3 * G
        pz_x, pz_y = aff(Pz)
        Q2 = rand_scalar() * G
        q2_x, q2_y = aff(Q2)
        Q2_aff = PointJacobi.from_affine(Q2.to_affine(), generator=G)
        R_mz = Pz + Q2_aff
        mz_x, mz_y = aff(R_mz)

        case = {
            "p_x": p_x, "p_y": p_y, "q_x": q_x, "q_y": q_y,
            "add_x": add_x, "add_y": add_y,
            "mix_x": mix_x, "mix_y": mix_y,
            "dbl_x": dbl_x, "dbl_y": dbl_y,
            "pz_x": pz_x, "pz_y": pz_y, "q2_x": q2_x, "q2_y": q2_y,
            "mz_x": mz_x, "mz_y": mz_y,
        }
        for key, (x, y) in [
            ("P", (p_x, p_y)),
            ("Q", (q_x, q_y)),
            ("add", (add_x, add_y)),
            ("mix", (mix_x, mix_y)),
            ("dbl", (dbl_x, dbl_y)),
            ("Pz", (pz_x, pz_y)),
            ("Q2", (q2_x, q2_y)),
            ("mixed_z", (mz_x, mz_y)),
        ]:
            check_xy(f"case {case_idx} {key}", x, y)
        cases.append(case)

    with open(path, "w") as f:
        f.write("#pragma once\n")
        f.write("// secp256k1 G1 EC vectors (ecdsa); Y²=X³+7; 6×45 RNS; all coords on-curve\n")
        f.write(f"constexpr int EC_NUM_CASES = {num_cases};\n\n")
        f.write("struct SecpEcCase {\n")
        f.write("    const char* p_x; const char* p_y;\n")
        f.write("    const char* q_x; const char* q_y;\n")
        f.write("    const char* add_x; const char* add_y;\n")
        f.write("    const char* mix_x; const char* mix_y;\n")
        f.write("    const char* dbl_x; const char* dbl_y;\n")
        f.write("    const char* pz_x; const char* pz_y;\n")
        f.write("    const char* q2_x; const char* q2_y;\n")
        f.write("    const char* mz_x; const char* mz_y;\n")
        f.write("};\n\n")
        f.write("constexpr SecpEcCase EC_CASES[] = {\n")
        for c in cases:
            f.write("    {")
            keys = ["p_x", "p_y", "q_x", "q_y", "add_x", "add_y", "mix_x", "mix_y",
                    "dbl_x", "dbl_y", "pz_x", "pz_y", "q2_x", "q2_y", "mz_x", "mz_y"]
            f.write(", ".join(f'"{c[k]}"' for k in keys))
            f.write("},\n")
        f.write("};\n")
    print(f"wrote {path} ({num_cases} cases, seed={seed_value})")


def intrns4_true_moduli(moduli, num_moduli):
    """Logical moduli (exclude mod-1 padding used when num_moduli < 8)."""
    if num_moduli < 8:
        return moduli[1:1 + num_moduli]
    return moduli[:num_moduli]


def print_intrns4_params(f, intrns, target, num_moduli, namespace, dim_prefix):
    """Write C++ constants for IntRNS4 (perm + scalar r2 matrices, rotation, moduli)."""
    n = num_moduli
    padded_dim = len(intrns.m1)
    mod1 = intrns4_true_moduli(intrns.m1, n)
    mod2 = intrns4_true_moduli(intrns.m2, n)
    skip_mod1_layout = n < 8
    f.write('#pragma once\n')
    f.write('// Generated by scripts/rns_secp256k1.py — do not edit\n\n')
    f.write('#include <array>\n\n')
    f.write(f'namespace {namespace} {{\n\n')
    f.write(f'constexpr int {dim_prefix}_RNS_DIM = {n};\n')
    f.write(f'constexpr int {dim_prefix}_PADDED_DIM = {padded_dim};\n')
    f.write(f'constexpr int {dim_prefix}_CONVERT_TO_DIGITS = {intrns.convert_to.indigits};\n')
    f.write(f'constexpr int {dim_prefix}_CONVERT_TO_RADIX_BITS = {intrns.convert_to.inbase};\n\n')

    f.write(f'constexpr std::array<uint64_t, {n}> moduli1 = {format_array(mod1)};\n')
    f.write(f'constexpr std::array<uint64_t, {n}> moduli2 = {format_array(mod2)};\n')
    f.write(f'constexpr std::array<uint64_t, {padded_dim}> moduli1_padded = {format_array(intrns.m1)};\n')
    f.write(f'constexpr std::array<uint64_t, {padded_dim}> moduli2_padded = {format_array(intrns.m2)};\n\n')

    f.write(f'const char* target_modulus = "{target}";\n')
    f.write(f'const char* M1 = "{intrns.factor}";\n')
    f.write(f'const char* M2 = "{intrns.factor2}";\n')
    f.write(f'const char* M1M2 = "{intrns.total_modulus}";\n')
    f.write(f'const char* total_rotation = "{intrns.total_rotation}";\n')
    f.write(f'const char* total_inv_rotation = "{intrns.total_inv_rotation}";\n')
    f.write(f'const char* qr_premult = "{intrns.qr_premult}";\n')
    f.write(f'const char* inv_factor_target = "{pow(intrns.factor, -1, target)}";\n\n')

    f.write('// Premult/postmult for IntRNS2System factory (matches Python IntRNS4)\n')
    f.write(f'const char* r1_premult = "{intrns.r1.premult}";\n')
    f.write(f'const char* r1_postmult = "{intrns.r1.postmult}";\n')
    f.write(f'const char* r2_premult = "{intrns.r2.premult}";\n')
    f.write(f'const char* r2_postmult = "{intrns.r2.postmult}";\n')
    f.write(f'const char* convert_to_premult = "{intrns.convert_to.conversion_matrix.premult}";\n')
    f.write(f'const char* convert_to_postmult = "{intrns.convert_to.conversion_matrix.postmult}";\n')
    f.write(f'const char* convert_from_premult = "{intrns.convert_from.conversion_matrix.premult}";\n')
    f.write(f'const char* convert_from_postmult = "{intrns.convert_from.conversion_matrix.postmult}";\n\n')

    if intrns.r1.perm_mat is not None:
        r1_perm = perm_mat_rows_to_8(intrns.r1.perm_mat)
        if n == 8:
            f.write('// r1 perm rows (pre_permute_mt): 8×8 cyclic advance_perm_mt\n')
        else:
            f.write('// r1 perm rows (pre_permute_mt_c): N×8 AVX; shift_perm from advance_perm_mt_c\n')
        f.write(f'constexpr std::array<std::array<uint64_t, 8>, {len(r1_perm)}> r1_perm_mat_fixed = {format_matrix_rows(r1_perm)};\n')
        f.write(f'constexpr bool r1_skip_correction = {"true" if intrns.r1.skip_correction else "false"};\n\n')

    if intrns.r2.perm_mat is not None:
        r2_perm = perm_mat_rows_to_8(intrns.r2.perm_mat)
        f.write(f'constexpr std::array<std::array<uint64_t, 8>, {len(r2_perm)}> r2_perm_mat_fixed = {format_matrix_rows(r2_perm)};\n\n')

    f.write(f'// r1/r2 base matrices ({n} rows × {padded_dim} AVX lanes)\n')
    print_rns_matrix_avx8(f, intrns.r1, "r1", skip_zero_input_row=skip_mod1_layout)
    print_rns_matrix_avx8(f, intrns.r2, "r2", skip_zero_input_row=skip_mod1_layout)

    f.write('\n// convert_to/from matrices\n')
    print_rns_matrix_avx8(f, intrns.convert_to.conversion_matrix, "convert_to")
    print_rns_matrix_avx8(f, intrns.convert_from.conversion_matrix, "convert_from")

    reducer1 = intrns.reducer1
    reducer2 = intrns.reducer2
    f.write('\n// MontgomeryReduce mont_factor (reducer1=m1, reducer2=m2)\n')
    f.write(f'constexpr std::array<uint64_t, {padded_dim}> reducer1_mont_factor = {format_array(reducer1.mont_factor.store())};\n')
    f.write(f'constexpr std::array<uint64_t, {padded_dim}> reducer2_mont_factor = {format_array(reducer2.mont_factor.store())};\n')

    f.write(f'\n}} // namespace {namespace}\n')


def print_bls12_381_intrns4_params(path="bls12_381_intrns4_params.hpp", test_path="../test_data/test_values_bls12_381.hpp"):
    """Generate BLS12-381 G1 field IntRNS4 export (8×50-bit moduli, r1 perm, scalar r2)."""
    target = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
    params = Params.avx_redundant(MontgomeryReduce)
    seed(5)
    a = randint(0, target)
    b = randint(0, target)
    intrns = IntRNS4(target, params, 8)
    a_m1, a_m2 = intrns.to_mont_avx(a)
    b_m1, b_m2 = intrns.to_mont_avx(b)
    c_m1, c_m2 = intrns.mulreduce(a_m1, a_m2, b_m1, b_m2)
    c = intrns.from_mont_avx(c_m2)
    correct = (a * b) % target
    with open(path, "w") as f:
        print_intrns4_params(
            f, intrns, target, 8, "bls12_381_intrns4", "BLS12_381")
    with open(test_path, "w") as f:
        f.write(f'const char* modulus = "{target}";\n')
        f.write(f'const char* a = "{a}";\n')
        f.write(f'const char* b = "{b}";\n')
        f.write(f'const char* c_correct = "{correct}";\n')
        f.write(f'const std::array<uint64_t, {len(a_m1.store())}> a_m1 = {format_array(a_m1.store())};\n')
        f.write(f'const std::array<uint64_t, {len(a_m2.store())}> a_m2 = {format_array(a_m2.store())};\n')
        f.write(f'const std::array<uint64_t, {len(b_m1.store())}> b_m1 = {format_array(b_m1.store())};\n')
        f.write(f'const std::array<uint64_t, {len(b_m2.store())}> b_m2 = {format_array(b_m2.store())};\n')
        f.write(f'const std::array<uint64_t, {len(c_m1.store())}> c_m1 = {format_array(c_m1.store())};\n')
        f.write(f'const std::array<uint64_t, {len(c_m2.store())}> c_m2 = {format_array(c_m2.store())};\n')
    print(f"wrote {path}")
    print(f"wrote {test_path}")
    print("moduli1:", intrns.m1)
    print("moduli2:", intrns.m2)
    print("qr_premult:", intrns.qr_premult)
    print("total_rotation:", intrns.total_rotation)
    if intrns.r1.perm_mat is not None:
        print("r1_perm_mat_fixed:")
        for row in perm_mat_rows_to_8(intrns.r1.perm_mat):
            print(" ", row)
    return intrns


def print_curve_intrns4_exports(params, target, namespace, dim_prefix, params_path, test_path, num_moduli=6, seed_value=5):
    """Generate IntRNS4 C++ params + mul test vectors for a 6-moduli curve."""
    seed(seed_value)
    a = randint(0, target)
    b = randint(0, target)
    intrns = IntRNS4(target, params, num_moduli)
    a_m1, a_m2 = intrns.to_mont_avx(a)
    b_m1, b_m2 = intrns.to_mont_avx(b)
    c_m1, c_m2 = intrns.mulreduce(a_m1, a_m2, b_m1, b_m2)
    c = intrns.from_mont_avx(c_m2)
    correct = (a * b) % target

    a_m1_list = a_m1.store()
    a_m2_list = a_m2.store()
    b_m1_list = b_m1.store()
    b_m2_list = b_m2.store()
    c_m1_list = c_m1.store()
    c_m2_list = c_m2.store()

    with open(params_path, "w") as f:
        print_intrns4_params(f, intrns, target, num_moduli, namespace, dim_prefix)
    print(f"wrote {params_path}")

    with open(test_path, "w") as f:
        f.write(f'const char* modulus = "{target}";\n')
        f.write(f'const char* a = "{a}";\n')
        f.write(f'const char* b = "{b}";\n')
        f.write(f'const char* c = "{c}";\n')
        f.write(f'const std::array<uint64_t, {len(a_m1_list)}> a_m1 = {format_array(a_m1_list)};\n')
        f.write(f'const std::array<uint64_t, {len(a_m2_list)}> a_m2 = {format_array(a_m2_list)};\n')
        f.write(f'const std::array<uint64_t, {len(b_m1_list)}> b_m1 = {format_array(b_m1_list)};\n')
        f.write(f'const std::array<uint64_t, {len(b_m2_list)}> b_m2 = {format_array(b_m2_list)};\n')
        f.write(f'const std::array<uint64_t, {len(c_m1_list)}> c_m1 = {format_array(c_m1_list)};\n')
        f.write(f'const std::array<uint64_t, {len(c_m2_list)}> c_m2 = {format_array(c_m2_list)};\n')
        f.write(f'const char* c_correct = "{correct}";\n')
        f.write(f'const char* total_rotation_py = "{intrns.total_rotation}";\n')
        f.write(f'const char* total_inv_rotation_py = "{intrns.total_inv_rotation}";\n')
        f.write(f'const std::array<uint64_t, {len(intrns.m1)}> moduli1 = {format_array(intrns.m1)};\n')
        f.write(f'const std::array<uint64_t, {len(intrns.m2)}> moduli2 = {format_array(intrns.m2)};\n')
        
        # Print RNS matrix fields
        f.write('\n')
        print_rns_matrix(f, intrns.r1, "r1")
        print_rns_matrix(f, intrns.r2, "r2")
        print_rns_matrix(f, intrns.convert_to.conversion_matrix, "convert_to")
        print_rns_matrix(f, intrns.convert_from.conversion_matrix, "convert_from")
        
        # Print premult/postmult values for comparison
        f.write('\n')
        f.write(f'// Premult/postmult values from Python\n')
        f.write(f'const char* r1_premult_py = "{intrns.r1.premult}";\n')
        f.write(f'const char* r1_postmult_py = "{intrns.r1.postmult}";\n')
        f.write(f'const char* r2_premult_py = "{intrns.r2.premult}";\n')
        f.write(f'const char* r2_postmult_py = "{intrns.r2.postmult}";\n')
        f.write(f'const char* convert_to_premult_py = "{intrns.convert_to.conversion_matrix.premult}";\n')
        f.write(f'const char* convert_to_postmult_py = "{intrns.convert_to.conversion_matrix.postmult}";\n')
        f.write(f'const char* convert_from_premult_py = "{intrns.convert_from.conversion_matrix.premult}";\n')
        f.write(f'const char* convert_from_postmult_py = "{intrns.convert_from.conversion_matrix.postmult}";\n')
        
        # Print MontgomeryReduce constants for reducer1 (moduli1)
        f.write('\n')
        f.write('// MontgomeryReduce constants for reducer1 (m1/moduli1)\n')
        reducer1 = intrns.reducer1
        reducer1_mont_factor = reducer1.mont_factor.store()
        # reducer1_t_squared = reducer1.t_squared.store()
        f.write(f'const std::array<uint64_t, {len(reducer1_mont_factor)}> reducer1_mont_factor = {format_array(reducer1_mont_factor)};\n')
        #f.write(f'const std::array<uint64_t, {len(reducer1_t_squared)}> reducer1_t_squared = {format_array(reducer1_t_squared)};\n')
        # f.write(f'const uint64_t reducer1_hi_mask = {reducer1.hi_mask.data[0]};\n')
        # f.write(f'const uint64_t reducer1_lo_mask = {reducer1.lo_mask.data[0]};\n')
        
        # Print MontgomeryReduce constants for reducer2 (moduli2)
        f.write('\n')
        f.write('// MontgomeryReduce constants for reducer2 (m2/moduli2)\n')
        reducer2 = intrns.reducer2
        reducer2_mont_factor = reducer2.mont_factor.store()
        # reducer2_t_squared = reducer2.t_squared.store()
        f.write(f'const std::array<uint64_t, {len(reducer2_mont_factor)}> reducer2_mont_factor = {format_array(reducer2_mont_factor)};\n')
        # f.write(f'const std::array<uint64_t, {len(reducer2_t_squared)}> reducer2_t_squared = {format_array(reducer2_t_squared)};\n')
        # f.write(f'const uint64_t reducer2_hi_mask = {reducer2.hi_mask.data[0]};\n')
        # f.write(f'const uint64_t reducer2_lo_mask = {reducer2.lo_mask.data[0]};\n')

    return intrns


def print_secp256k1_intrns4_exports(params, target=SECP256K1_MODULUS):
    return print_curve_intrns4_exports(
        params, target, "secp256k1_intrns4", "SECP256K1",
        "secp256k1_intrns4_params.hpp", "../test_data/test_values_secp256k1.hpp")


def print_bn254_intrns4_exports(params, target=BN254_MODULUS):
    return print_curve_intrns4_exports(
        params, target, "bn254_intrns4", "BN254",
        "bn254_intrns4_params.hpp", "../test_data/test_values_bn254.hpp")


if __name__ == "__main__":
    # print_tests()

    # Diagonal matvec (independent of curve); needs mulbits set on AVXVector
    # _p = Params(45, 52, 52, 32, 2, MontgomeryReduce)
    # test_diagonal_matmul(n=6, trials=500, seed_value=42)
    # test_diagonal_matmul(n=8, trials=200, seed_value=43)
    # test_diagonal_matmul_mt(n=8, trials=500, seed_value=48)
    # test_diagonal_matmul_mt(n=6, width=8, trials=500, seed_value=49)
    # test_diagonal_matmul_mt_c(n=6, trials=500, seed_value=44)
    # test_diagonal_matmul_mt_c(n=8, trials=200, seed_value=45)
    # test_diagonal_matmul_mt_c_last(n=6, trials=500, seed_value=46)
    # test_diagonal_matmul_mt_c_last(n=8, trials=200, seed_value=47)

    p = Params(46, 52, 52, 52, 1, MontgomeryReduce)

    print_secp256k1_intrns4_exports(p)
    print_bn254_intrns4_exports(p)
    print_bn254_ec_tests("../test_data/test_values_bn254_ec.hpp", num_cases=50, seed_value=42)
    # print_ec_tests("../test_data/test_values_secp256k1_ec.hpp", num_cases=50, seed_value=42)
    test_mulreduction(IntRNS4, p, SECP256K1_MODULUS, 10, num_moduli=6)
    test_mulreduction(IntRNS4, p, BN254_MODULUS, 10, num_moduli=6)

    target = BLS12_381_MODULUS
    cpu_params_50 = Params.avx_redundant(MontgomeryReduce)
    test_mulreduction(IntRNS4, cpu_params_50, target, 10, num_moduli=8)
    print_bls12_381_intrns4_params()
    
    