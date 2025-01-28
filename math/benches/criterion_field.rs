use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{Output, PProfProfiler};

mod fields;
use fields::mersenne31::{mersenne31_extension_ops_benchmarks, mersenne31_ops_benchmarks};
use fields::mersenne31_montgomery::mersenne31_mont_ops_benchmarks;
use fields::{
    baby_bear::{
        babybear_extension_ops_benchmarks_p3, babybear_p3_ops_benchmarks,
        babybear_u32_ops_benchmarks, babybear_u64_extension_ops_benchmarks,
        babybear_u64_ops_benchmarks,
    },
    stark252::starkfield_ops_benchmarks,
    u64_goldilocks::u64_goldilocks_ops_benchmarks,
    u64_goldilocks_montgomery::u64_goldilocks_montgomery_ops_benchmarks,
};

criterion_group!(
    name = field_benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));

    targets =babybear_u32_ops_benchmarks,babybear_u64_ops_benchmarks,babybear_u64_extension_ops_benchmarks, babybear_p3_ops_benchmarks,babybear_extension_ops_benchmarks_p3,mersenne31_extension_ops_benchmarks,mersenne31_ops_benchmarks,
    starkfield_ops_benchmarks,u64_goldilocks_ops_benchmarks,u64_goldilocks_montgomery_ops_benchmarks,mersenne31_mont_ops_benchmarks
);
criterion_main!(field_benches);
