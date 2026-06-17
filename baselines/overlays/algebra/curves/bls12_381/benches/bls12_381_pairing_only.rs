use ark_algebra_bench_templates::*;
use ark_bls12_381::{Bls12_381};

// Only generate pairing benchmarks
pairing_bench!(Bls12_381);

criterion_main!(pairing::benches);
