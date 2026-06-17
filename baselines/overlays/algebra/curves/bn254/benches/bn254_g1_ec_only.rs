use ark_algebra_bench_templates::{criterion, criterion_group, criterion_main};
use ark_bn254::G1Projective as G1;
use ark_ec::CurveGroup;
use ark_ff::AdditiveGroup;
use ark_std::UniformRand;

fn g1_arithmetic(c: &mut criterion::Criterion) {
    const SAMPLES: usize = 1000;
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group("Arithmetic for BN254::G1Projective");

    let left = (0..SAMPLES)
        .map(|_| G1::rand(&mut rng))
        .collect::<Vec<_>>();
    let right = (0..SAMPLES)
        .map(|_| G1::rand(&mut rng))
        .collect::<Vec<_>>();
    let right_affine = G1::normalize_batch(&right);

    group.bench_function("Addition", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            left[i] + right[i]
        })
    });
    group.bench_function("Mixed Addition", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            left[i] + right_affine[i]
        })
    });
    group.bench_function("Double", |b| {
        let mut i = 0;
        b.iter(|| {
            i = (i + 1) % SAMPLES;
            left[i].double()
        })
    });
    group.finish();
}

criterion_group!(benches, g1_arithmetic);
criterion_main!(benches);
