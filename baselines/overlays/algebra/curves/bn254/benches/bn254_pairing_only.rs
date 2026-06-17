use ark_algebra_bench_templates::*;
use ark_bn254::Bn254;

// Only generate pairing benchmarks
pairing_bench!(Bn254);

criterion_main!(pairing::benches);
