# src

Main C++ templates: RNS Montgomery multiplication, change-base (CRNS), pairing, and EC arithmetic. Build with `make`; see root **README.md** and **ARTIFACT.md**. Integer fallback (debugging): `make -f Makefile.fallback`.

## Precomputed ring parameters

Production curves ship with IntRNS4 exports already wired in:

| Curve | Typedefs | System / perm headers |
|-------|----------|------------------------|
| BLS12-381 (50-bit, MatrixNoK) | `bls12_381_ring_params.hpp` → `BlsRingProduction` | `bls12_381_intrns4_system.hpp`, `bls12_381_perm*.hpp` |
| BN254 (46-bit, FixedPerm) | `bn254_ring_params.hpp` → `Bn254RingProduction` | `bn254_intrns4_system.hpp`, `bn254_perm*.hpp` |

Generated constants originate in `scripts/rns_secp256k1.py` → `scripts/*_intrns4_params.hpp`. Regeneration: **ARTIFACT.md** § Parameter generation.

## Module map (paper terminology)

| Paper / concept | Primary files | Notes |
|-----------------|---------------|-------|
| **IntRNS2 / IntRNS4** | `bounded_ring.hpp`, `*_intrns4_system.hpp` | Two moduli sets m1/m2; Montgomery RNS multiply |
| **CRNS / change-base** | `change_base.hpp`, `cpu/vector/changebase.hpp` | r1/r2 base conversion between residue sets |
| **MatrixNoK (“QR no k”)** | `bls12_381_perm*.hpp`, `matrix_nok_batch.hpp` | Advanced params; 50-bit convert radix |
| **FixedPerm** | `bn254_perm*.hpp`, `bn254_perm_batch.hpp` | Precomputed r1 permutation (BN254) |
| **Delayed reduction** | `ring_element.hpp`, `elementwise_reduce.hpp` | `DelayedProduct`, `reduce_auto` |
| **Integer ↔ RNS** | `cpu/vector/conversion.hpp`, `bounded_ring.hpp` | `from_bigint`, `to_bigint`, `convert_to_rns` |
| **FP2 / FP12 tower** | `fp2.hpp`, `fp12.hpp` | Extension-field arithmetic in RNS |
| **Miller loop / pairing** | `miller.hpp`, `pairing.hpp`, `final_exponentiation.hpp` | BLS12-381 pairing stack |
| **EC / scalar mult** | `ec.hpp`, `scalar_mult.hpp` | G1/G2 with GLV/GLS |

## Layout

- **`*.hpp`** — library headers (included by tests/benches, not a separate `.so`).
- **`test_*.cpp`** — correctness tests (`make test-*`).
- **`bench_*.cpp`** — Google Benchmark harnesses (`make bench-*`).
- **`gmp_wrapper.hpp`** — GMP `BigInt` + `test_data/test_values.hpp` for reference arithmetic.

## Build targets

See **ARTIFACT.md** § Full test and benchmark matrix. Quick checks:

```bash
make test                          # reduce + inversion
make test-bls12-381-field-mul      # single field multiply
make test-bls12-381-pairing        # full pairing
make bench-bls12-381              # paper CPU benchmarks (BLS12-381)
```

Examples: `make -C ../examples` (see `examples/README.md`).
