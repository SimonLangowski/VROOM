# Examples

Minimal programs demonstrating core VROOM APIs. These are the artifact “minimal working examples.”

| Program | What it shows |
|---------|----------------|
| `01_parameter_setup` | BLS12-381 and BN254 ring configurations: limb counts, residue bit widths, moduli m1/m2 |
| `02_singular_modmul` | One field multiply: `from_bigint` → RNS, `modmul`, `to_bigint` |
| `03_sum_of_products` | Delayed multiply (`operator*`), accumulate, `batch_reduce_expand` |
| `04_difference_of_products` | `(a*b) - (c*d)` via `ring.negate` on a multiplicand (no delayed subtract) |
| `05_product_of_sums` | `(a+b)*(c+d)` via `ring.prep` / `ring.prep_left` before delayed multiply |

## Build & run

From the repo root (after `cd blst && make` if you need BLST elsewhere — examples do not require BLST):

```bash
make -C examples
./examples/01_parameter_setup
./examples/02_singular_modmul
./examples/03_sum_of_products
./examples/04_difference_of_products
./examples/05_product_of_sums
```

Requires **AVX-512 IFMA** (same as the main `src/` build). On a machine without IFMA, use the integer fallback:

```bash
make -C examples FALLBACK=1
```

Example 2 links against GMP for reference integers.

## Related code

- Ring typedefs: `src/bls12_381_ring_params.hpp`, `src/bn254_ring_params.hpp`
- Parameter generation: `scripts/rns_secp256k1.py`
- Full tests: `src/test_elementwise_reduce.cpp`, `src/test_bls12_381_field_mul.cpp`
