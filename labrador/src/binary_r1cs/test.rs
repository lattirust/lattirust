use tracing_subscriber::fmt::format;
use tracing_subscriber::fmt::format::FmtSpan;

use lattirust_arithmetic::ring::Zq1;
use lattirust_arithmetic::ring::ntt::ntt_prime;
use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;
use relations::{test_completeness_with_init, test_soundness_with_init};
use relations::r1cs::Size;
use relations::reduction::Reduction;

use crate::binary_r1cs::ReductionBinaryR1CSPrincipalRelation;
use crate::binary_r1cs::util::BinaryR1CSCRS;

const Q: u64 = ntt_prime::<64>(32);
const D: usize = 64;
type F = Zq1<Q>;

type R = Pow2CyclotomicPolyRingNTT<F, D>;

fn init() {
    let _ = tracing_subscriber::fmt::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
        .event_format(format().compact())
        .with_env_filter("none,labrador=trace")
        .try_init();
}

type TestReduction = ReductionBinaryR1CSPrincipalRelation<R>;

const TEST_SIZE: Size = Size {
    num_constraints: D * 4,
    num_instance_variables: D,
    num_witness_variables: D * 3,
};

test_completeness_with_init!(
    TestReduction,
    BinaryR1CSCRS::new(TEST_SIZE.num_constraints, TEST_SIZE.num_instance_variables+TEST_SIZE.num_witness_variables),
    TEST_SIZE,
    init
);

test_soundness_with_init!(
    TestReduction,
    BinaryR1CSCRS::new(TEST_SIZE.num_constraints, TEST_SIZE.num_instance_variables+TEST_SIZE.num_witness_variables),
    TEST_SIZE,
    init
);
