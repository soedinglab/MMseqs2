use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;
use simulate_seqs::*;

use std::str;

fn consistent(i: usize, j: usize, cigar: &Cigar) -> bool {
    let mut curr_i = 0;
    let mut curr_j = 0;

    for i in 0..cigar.len() {
        let op_len = cigar.get(i);
        match op_len.op {
            Operation::M => {
                curr_i += op_len.len;
                curr_j += op_len.len;
            },
            Operation::I => {
                curr_i += op_len.len;
            },
            _ => {
                curr_j += op_len.len;
            }
        }
    }

    curr_i == i && curr_j == j
}

fn test(iter: usize, len: usize, k: usize, insert_len: Option<usize>) -> usize {
    let mut wrong = 0usize;
    let mut rng = StdRng::seed_from_u64(1234);

    for _i in 0..iter {
        let r = rand_str(len, &AMINO_ACIDS, &mut rng);
        let q = match insert_len {
            Some(len) => rand_mutate_insert(&r, k, &AMINO_ACIDS, len, &mut rng),
            None => rand_mutate(&r, k, &AMINO_ACIDS, &mut rng)
        };

        let r_padded = PaddedBytes::from_bytes::<AAMatrix>(&r, 2048);
        let q_padded = PaddedBytes::from_bytes::<AAMatrix>(&q, 2048);
        let run_gaps = Gaps { open: -11, extend: -1 };

        let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), 2048);
        block_aligner.align(&q_padded, &r_padded, &BLOSUM62, run_gaps, 32..=2048, 0);
        let scan_score = block_aligner.res().score;
        let mut scan_cigar = Cigar::new(q.len(), r.len());
        block_aligner.trace().cigar(q.len(), r.len(), &mut scan_cigar);

        if !consistent(q.len(), r.len(), &scan_cigar) {
            wrong += 1;

            println!(
                "score: {}\nq: {}\nr: {}\ncigar: {}",
                scan_score,
                str::from_utf8(&q).unwrap(),
                str::from_utf8(&r).unwrap(),
                scan_cigar
            );
        }
    }

    wrong
}

fn main() {
    let iters = [100, 100, 100];
    let lens = [10, 100, 1000];
    let rcp_ks = [10.0, 5.0, 2.0];
    let inserts = [false, true];

    let mut total_wrong = 0usize;
    let mut total = 0usize;

    for (&len, &iter) in lens.iter().zip(&iters) {
        for &rcp_k in &rcp_ks {
            for &insert in &inserts {
                let insert_len = if insert { Some(len / 10) } else { None };
                let wrong = test(iter, len, ((len as f64) / rcp_k) as usize, insert_len);
                println!(
                    "\nlen: {}, k: {}, insert: {}, iter: {}, wrong: {}\n",
                    len,
                    ((len as f64) / rcp_k) as usize,
                    insert,
                    iter,
                    wrong
                );
                total_wrong += wrong;
                total += iter;
            }
        }
    }

    println!("\ntotal: {}, wrong: {}", total, total_wrong);
    println!("Done!");
}
