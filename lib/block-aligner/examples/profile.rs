use block_aligner::scan_block::*;
use block_aligner::scores::*;
use simulate_seqs::*;

use std::hint::black_box;

fn run(len: usize, k: usize) {
    let mut rng = StdRng::seed_from_u64(1234);
    let r = rand_str(len, &AMINO_ACIDS, &mut rng);
    let q = rand_mutate(&r, k, &AMINO_ACIDS, &mut rng);
    let r = PaddedBytes::from_bytes::<AAMatrix>(&r, 2048);
    let q = PaddedBytes::from_bytes::<AAMatrix>(&q, 2048);
    let run_gaps = Gaps { open: -11, extend: -1 };
    let mut a = Block::<true, true>::new(q.len(), r.len(), 32);

    for _i in 0..10000 {
        a.align(&q, &r, &BLOSUM62, run_gaps, 32..=32, 1000);
        black_box(a.res());
    }
}

fn main() {
    run(10000, 1000);
}
