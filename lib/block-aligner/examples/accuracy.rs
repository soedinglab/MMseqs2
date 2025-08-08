//use parasailors::{Matrix, *};

use bio::alignment::pairwise::*;
use bio::scores::blosum62;

use block_aligner::scan_block::*;
use block_aligner::scores::*;
use simulate_seqs::*;

use std::{env, str, cmp};

fn test(iter: usize, len: usize, k: usize, slow: bool, insert_len: Option<usize>, nuc: bool, max_size: usize, verbose: bool) -> (usize, f64, i32, i32) {
    let mut wrong = 0usize;
    let mut wrong_avg = 0f64;
    let mut wrong_min = i32::MAX;
    let mut wrong_max = i32::MIN;
    let mut rng = StdRng::seed_from_u64(1234);
    let nw = |a, b| if a == b { 1 } else { -1 };

    for _i in 0..iter {
        let r = rand_str(len, if nuc { &NUC } else { &AMINO_ACIDS }, &mut rng);
        let q = match insert_len {
            Some(len) => rand_mutate_insert(&r, k, if nuc { &NUC } else { &AMINO_ACIDS }, len, &mut rng),
            None => rand_mutate(&r, k, if nuc { &NUC } else { &AMINO_ACIDS }, &mut rng)
        };

        // rust-bio
        let bio_score = if nuc {
            let mut bio_aligner = Aligner::with_capacity(q.len(), r.len(), -1, -1, &nw);
            bio_aligner.global(&q, &r).score
        } else {
            let mut bio_aligner = Aligner::with_capacity(q.len(), r.len(), -10, -1, &blosum62);
            bio_aligner.global(&q, &r).score
        };

        // parasailors
        /*
        let matrix = Matrix::new(MatrixType::Blosum62);
        let profile = Profile::new(&q, &matrix);
        let parasail_score = global_alignment_score(&profile, &r, 11, 1);
        */

        // ours
        let scan_score = if slow {
            slow_align(&q, &r)
        } else {
            if nuc {
                let run_gaps = Gaps { open: -2, extend: -1 };
                let r_padded = PaddedBytes::from_bytes::<NucMatrix>(&r, 2048);
                let q_padded = PaddedBytes::from_bytes::<NucMatrix>(&q, 2048);
                let mut block_aligner = Block::<false, false>::new(q.len(), r.len(), max_size);
                block_aligner.align(&q_padded, &r_padded, &NW1, run_gaps, 32..=max_size, 0);
                block_aligner.res().score
            } else {
                let run_gaps = Gaps { open: -11, extend: -1 };
                let r_padded = PaddedBytes::from_bytes::<AAMatrix>(&r, 2048);
                let q_padded = PaddedBytes::from_bytes::<AAMatrix>(&q, 2048);
                let mut block_aligner = Block::<false, false>::new(q.len(), r.len(), max_size);
                block_aligner.align(&q_padded, &r_padded, &BLOSUM62, run_gaps, 32..=max_size, 0);
                block_aligner.res().score
            }
        };

        if bio_score != scan_score {
            wrong += 1;
            let score_diff = bio_score - scan_score;
            wrong_avg += (score_diff as f64) / (bio_score as f64);
            wrong_min = cmp::min(wrong_min, score_diff);
            wrong_max = cmp::max(wrong_max, score_diff);

            if verbose {
                println!(
                    "bio: {}, ours: {}\nq: {}\nr: {}",
                    bio_score,
                    scan_score,
                    str::from_utf8(&q).unwrap(),
                    str::from_utf8(&r).unwrap()
                );
            }
        }
    }

    (wrong, wrong_avg / (wrong as f64), wrong_min, wrong_max)
}

fn main() {
    let arg1 = env::args().skip(1).next();
    let slow = false;
    let nuc = true;
    let verbose = arg1.is_some() && arg1.unwrap() == "-v";
    let iters = [100, 100, 10];
    let lens = [100, 1000, 10000];
    let rcp_ks = [10.0, 5.0, 2.0];
    let inserts = [false, true];
    let max_sizes = [32, 2048];

    let mut total_wrong = 0usize;
    let mut total = 0usize;

    println!("\nlen, k, insert, iter, max size, wrong, wrong % error, wrong min, wrong max\n");

    for (&len, &iter) in lens.iter().zip(&iters) {
        for &rcp_k in &rcp_ks {
            for &insert in &inserts {
                for &max_size in &max_sizes {
                    let insert_len = if insert { Some(len / 10) } else { None };
                    let (wrong, wrong_avg, wrong_min, wrong_max) = test(iter, len, ((len as f64) / rcp_k) as usize, slow, insert_len, nuc, max_size, verbose);
                    println!(
                        "\n{}, {}, {}, {}, {}, {}, {}, {}, {}\n",
                        len,
                        ((len as f64) / rcp_k) as usize,
                        insert,
                        iter,
                        max_size,
                        wrong,
                        wrong_avg,
                        wrong_min,
                        wrong_max
                    );
                    total_wrong += wrong;
                    total += iter;
                }
            }
        }
    }

    println!("\n# total: {}, wrong: {}", total, total_wrong);
    println!("# Done!");
}

// Scalar version of the block aligner algorithm for testing
// purposes. May not exactly match the implementation of the
// vectorized version.
#[allow(non_snake_case)]
fn slow_align(q: &[u8], r: &[u8]) -> i32 {
    let mut block_width = 16usize;
    let mut block_height = 16usize;
    let block_grow = 16usize;
    let max_size = 256usize;
    let i_step = 4usize;
    let j_step = 4usize;
    let mut y_drop = 3i32;
    let y_drop_grow = 2i32;

    let mut D = vec![i32::MIN; (q.len() + 1 + max_size) * (r.len() + 1 + max_size)];
    let mut R = vec![i32::MIN; (q.len() + 1 + max_size) * (r.len() + 1 + max_size)];
    let mut C = vec![i32::MIN; (q.len() + 1 + max_size) * (r.len() + 1 + max_size)];
    D[0 + 0 * (q.len() + 1 + max_size)] = 0;
    let mut i = 0usize;
    let mut j = 0usize;
    let mut dir = 0;
    let mut best_max = 0;

    //println!("start");

    loop {
        let max = match dir {
            0 => { // right
                calc_block(q, r, &mut D, &mut R, &mut C, i, j, block_width, block_height, max_size, -11, -1)
            },
            _ => { // down
                calc_block(q, r, &mut D, &mut R, &mut C, i, j, block_width, block_height, max_size, -11, -1)
            }
        };

        if i + block_height > q.len() && j + block_width > r.len() {
            break;
        }

        let right_max = block_max(&D, q.len() + 1 + max_size, i, j + block_width - 1, 1, block_height);
        let down_max = block_max(&D, q.len() + 1 + max_size, i + block_height - 1, j, block_width, 1);
        best_max = cmp::max(best_max, max);

        if block_width < max_size && cmp::max(right_max, down_max) < best_max - y_drop {
            block_width += block_grow;
            block_height += block_grow;
            y_drop += y_drop_grow;
            //println!("i: {}, j: {}, w: {}", i, j, block_width);
            continue;
        }

        if j + block_width > r.len() {
            i += i_step;
            dir = 1;
        } else if i + block_height > q.len() {
            j += j_step;
            dir = 0;
        } else {
            if down_max > right_max {
                i += i_step;
                dir = 1;
            } else if right_max > down_max {
                j += j_step;
                dir = 0;
            } else {
                j += j_step;
                dir = 0;
            }
        }
    }

    D[q.len() + r.len() * (q.len() + 1 + max_size)]
}

#[allow(non_snake_case)]
fn block_max(D: &[i32], col_len: usize, start_i: usize, start_j: usize, block_width: usize, block_height: usize) -> i32 {
    let mut max = i32::MIN;
    for i in start_i..start_i + block_height {
        for j in start_j..start_j + block_width {
            max = cmp::max(max, D[i + j * col_len]);
        }
    }
    max
}

#[allow(non_snake_case)]
fn calc_block(q: &[u8], r: &[u8], D: &mut [i32], R: &mut [i32], C: &mut [i32], start_i: usize, start_j: usize, block_width: usize, block_height: usize, max_size: usize, gap_open: i32, gap_extend: i32) -> i32 {
    let idx = |i: usize, j: usize| { i + j * (q.len() + 1 + max_size) };
    let mut max = i32::MIN;

    for i in start_i..start_i + block_height {
        for j in start_j..start_j + block_width {
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
            max = cmp::max(max, D[idx(i, j)]);
        }
    }

    max
}
