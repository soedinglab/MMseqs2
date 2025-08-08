use bio::scores::blosum62;

use block_aligner::scan_block::*;
use block_aligner::scores::*;
use simulate_seqs::*;

use std::{env, str, cmp};

fn test(iter: usize, len: usize, k: usize, verbose: bool) -> (usize, f64, i32, i32, usize) {
    let mut wrong = 0usize;
    let mut wrong_avg = 0f64;
    let mut wrong_min = i32::MAX;
    let mut wrong_max = i32::MIN;
    let mut diff_idx = 0usize;
    let mut rng = StdRng::seed_from_u64(1234);

    for _i in 0..iter {
        let mut r = rand_str(len, &AMINO_ACIDS, &mut rng);
        let q = rand_mutate_suffix(&mut r, k, &AMINO_ACIDS, 500, &mut rng);

        let r_padded = PaddedBytes::from_bytes::<AAMatrix>(&r, 2048);
        let q_padded = PaddedBytes::from_bytes::<AAMatrix>(&q, 2048);
        let run_gaps = Gaps { open: -11, extend: -1 };

        let slow_res = slow_align(&q, &r, 50);

        let mut block_aligner = Block::<false, true>::new(q.len(), r.len(), 64);
        block_aligner.align(&q_padded, &r_padded, &BLOSUM62, run_gaps, 32..=64, 50);
        let scan_res = block_aligner.res();

        if slow_res.0 != scan_res.score {
            wrong += 1;
            let score_diff = slow_res.0 - scan_res.score;
            wrong_avg += (score_diff as f64) / (slow_res.0 as f64);
            wrong_min = cmp::min(wrong_min, score_diff);
            wrong_max = cmp::max(wrong_max, score_diff);

            if verbose {
                println!(
                    "slow: (score: {}, i: {}, j: {}),\nours: (score: {}, i: {}, j: {})\nq: {}\nr: {}",
                    slow_res.0,
                    slow_res.1,
                    slow_res.2,
                    scan_res.score,
                    scan_res.query_idx,
                    scan_res.reference_idx,
                    str::from_utf8(&q).unwrap(),
                    str::from_utf8(&r).unwrap()
                );
            }
        }

        if slow_res.1 != scan_res.query_idx || slow_res.2 != scan_res.reference_idx {
            diff_idx += 1;

            if verbose {
                println!(
                    "slow: (i: {}, j: {}),\nours: (i: {}, j: {})\nq: {}\nr: {}",
                    slow_res.1,
                    slow_res.2,
                    scan_res.query_idx,
                    scan_res.reference_idx,
                    str::from_utf8(&q).unwrap(),
                    str::from_utf8(&r).unwrap()
                );
            }
        }
    }

    (wrong, wrong_avg / (wrong as f64), wrong_min, wrong_max, diff_idx)
}

fn main() {
    let arg1 = env::args().skip(1).next();
    let verbose = arg1.is_some() && arg1.unwrap() == "-v";
    let iters = [100, 100, 100];
    let lens = [10, 100, 1000];
    let rcp_ks = [10.0, 5.0, 2.0];

    let mut total_wrong = 0usize;
    let mut total = 0usize;
    let mut total_diff_idx = 0usize;

    for (&len, &iter) in lens.iter().zip(&iters) {
        for &rcp_k in &rcp_ks {
            let (wrong, wrong_avg, wrong_min, wrong_max, diff_idx) = test(iter, len, ((len as f64) / rcp_k) as usize, verbose);
            println!(
                "\nlen: {}, k: {}, iter: {}, wrong: {}, wrong % error: {}, wrong min: {}, wrong max: {}, diff idx: {}\n",
                len,
                ((len as f64) / rcp_k) as usize,
                iter,
                wrong,
                wrong_avg,
                wrong_min,
                wrong_max,
                diff_idx
            );
            total_wrong += wrong;
            total += iter;
            total_diff_idx += diff_idx;
        }
    }

    println!("\ntotal: {}, wrong: {}, diff idx: {}", total, total_wrong, total_diff_idx);
    println!("Done!");
}

#[allow(non_snake_case)]
fn slow_align(q: &[u8], r: &[u8], x_drop: i32) -> (i32, usize, usize) {
    let gap_open = -11;
    let gap_extend = -1;
    let idx = |i: usize, j: usize| { i + j * (q.len() + 1) };

    let mut D = vec![i32::MIN; (q.len() + 1) * (r.len() + 1)];
    let mut R = vec![i32::MIN; (q.len() + 1) * (r.len() + 1)];
    let mut C = vec![i32::MIN; (q.len() + 1) * (r.len() + 1)];
    D[idx(0, 0)] = 0;

    let mut best_max = i32::MIN;
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 0..=q.len() {
        let mut max = i32::MIN;
        let mut max_j = 0;
        for j in 0..=r.len() {
            if D[idx(i, j)] != i32::MIN {
                continue;
            }
            R[idx(i, j)] = if i == 0 { i32::MIN } else { cmp::max(
                R[idx(i - 1, j)].saturating_add(gap_extend),
                D[idx(i - 1, j)].saturating_add(gap_open)
            ) };
            C[idx(i, j)] = if j == 0 { i32::MIN } else { cmp::max(
                C[idx(i, j - 1)].saturating_add(gap_extend),
                D[idx(i, j - 1)].saturating_add(gap_open)
            ) };
            D[idx(i, j)] = cmp::max(
                if i == 0 || j == 0 || i > q.len() || j > r.len() { i32::MIN } else {
                    D[idx(i - 1, j - 1)].saturating_add(blosum62(q[i - 1], r[j - 1]))
                },
                cmp::max(R[idx(i, j)], C[idx(i, j)])
            );
            if D[idx(i, j)] > max {
                max = D[idx(i, j)];
                max_j = j;
            }
        }
        if max > best_max {
            best_max = max;
            best_i = i;
            best_j = max_j;
        }
        if max < best_max - x_drop {
            break;
        }
    }

    (best_max, best_i, best_j)
}
