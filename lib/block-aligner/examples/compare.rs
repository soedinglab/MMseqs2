use block_aligner::scan_block::*;
use block_aligner::scores::*;

use std::{env, cmp};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn test(file_name: &str, max_size: usize, x_drop: i32) -> (usize, usize, f64, usize, f64) {
    let reader = BufReader::new(File::open(file_name).unwrap());
    let mut count = 0;
    let mut other_better = 0;
    let mut other_better_avg = 0f64;
    let mut us_better = 0;
    let mut us_better_avg = 0f64;
    //let mut slow_better = 0;
    //let mut slow_equal = 0;

    for line in reader.lines() {
        let line = line.unwrap();
        let mut row = line.split_ascii_whitespace().take(5);
        let q = row.next().unwrap().to_ascii_uppercase();
        let r = row.next().unwrap().to_ascii_uppercase();
        let other_score = row.next().unwrap().parse::<i32>().unwrap();
        let _other_i = row.next().unwrap().parse::<usize>().unwrap();
        let _other_j = row.next().unwrap().parse::<usize>().unwrap();

        //let x_drop = 100;
        //let x_drop = 50;
        //let matrix = NucMatrix::new_simple(2, -3);
        let matrix = NucMatrix::new_simple(1, -1);
        let r_padded = PaddedBytes::from_bytes::<NucMatrix>(r.as_bytes(), 2048);
        let q_padded = PaddedBytes::from_bytes::<NucMatrix>(q.as_bytes(), 2048);
        //let run_gaps = Gaps { open: -5, extend: -1 };
        let run_gaps = Gaps { open: -2, extend: -1 };

        // ours
        let mut block_aligner = Block::<true, true>::new(q.len(), r.len(), max_size);
        block_aligner.align(&q_padded, &r_padded, &matrix, run_gaps, 32..=max_size, x_drop);
        let scan_res = block_aligner.res();
        let scan_score = scan_res.score;

        if scan_score > other_score {
            us_better += 1;
            us_better_avg += ((scan_score - other_score) as f64) / (scan_score as f64);
        }

        if scan_score < other_score {
            other_better += 1;
            other_better_avg += ((other_score - scan_score) as f64) / (other_score as f64);

            /*let slow_score = slow_align(q.as_bytes(), r.as_bytes(), x_drop);
            if slow_score > other_score {
                slow_better += 1;
            }
            if slow_score == other_score {
                slow_equal += 1;
            }
            println!("ours: {}, other: {}, slow: {}", scan_score, other_score, slow_score);*/
        }

        count += 1;
    }
    //println!("slow better: {}, slow equal: {}", slow_better, slow_equal);

    (count, other_better, other_better_avg / (other_better as f64), us_better, us_better_avg / (us_better as f64))
}

fn main() {
    let mut args = env::args().skip(1);
    let other_file = args.next().expect("Pass in the path to a tab-separated file to compare to!");
    let x_drop = args.next().expect("Pass in an X-drop threshold!").parse::<i32>().unwrap();
    let max_sizes = [32, 64];

    println!("max size, total, other better, other % better, us better, us % better");

    for &max_size in &max_sizes {
        let (count, other_better, other_better_avg, us_better, us_better_avg) = test(&other_file, max_size, x_drop);

        println!(
            "\n{}, {}, {}, {}, {}, {}",
            max_size,
            count,
            other_better,
            other_better_avg,
            us_better,
            us_better_avg
        );
    }

    println!("# Done!");
}

// Scalar version of the block aligner algorithm for testing
// purposes. May not exactly match the implementation of the
// vectorized version.
//
// Also possible to simulate diagonal adaptive banding methods.
#[allow(dead_code)]
#[allow(non_snake_case)]
fn slow_align(q: &[u8], r: &[u8], x_drop: i32) -> i32 {
    let block_size = 32usize;
    let step = 8usize;
    //let step = 1usize;

    let mut D = vec![i32::MIN; (q.len() + 1 + block_size) * (r.len() + 1 + block_size)];
    let mut R = vec![i32::MIN; (q.len() + 1 + block_size) * (r.len() + 1 + block_size)];
    let mut C = vec![i32::MIN; (q.len() + 1 + block_size) * (r.len() + 1 + block_size)];
    D[0 + 0 * (q.len() + 1 + block_size)] = 0;
    //let max = calc_block(q, r, &mut D, &mut R, &mut C, 0, 0, block_size, block_size, block_size, -2, -1);
    let mut i = 0usize;
    let mut j = 0usize;
    let mut dir = 0;
    //let mut best_max = max;
    let mut best_max = 0;

    loop {
        let max = match dir {
            0 => { // right
                calc_block(q, r, &mut D, &mut R, &mut C, i, j, block_size, block_size, block_size, -2, -1)
                //calc_diag(q, r, &mut D, &mut R, &mut C, i, j, block_size, -2, -1)
            },
            _ => { // down
                calc_block(q, r, &mut D, &mut R, &mut C, i, j, block_size, block_size, block_size, -2, -1)
                //calc_diag(q, r, &mut D, &mut R, &mut C, i, j, block_size, -2, -1)
            }
        };

        //let max = block_max(&D, q.len() + 1 + block_size, i + block_size / 2 - 1, j + block_size / 2, 1, 1);
        let right_max = block_sum(&D, q.len() + 1 + block_size, i, j + block_size - 1, 1, step);
        let down_max = block_sum(&D, q.len() + 1 + block_size, i + block_size - 1, j, step, 1);
        best_max = cmp::max(best_max, max);

        if max < best_max - x_drop {
            return best_max;
        }

        if i + block_size > q.len() && j + block_size > r.len() {
            return best_max;
        }

        if j + block_size > r.len() {
            i += step;
            dir = 1;
            continue;
        }
        if i + block_size > q.len() {
            j += step;
            dir = 0;
            continue;
        }
        if down_max > right_max {
            i += step;
            dir = 1;
        } else {
            j += step;
            dir = 0;
        }
    }
}

#[allow(dead_code)]
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

#[allow(dead_code)]
#[allow(non_snake_case)]
fn block_sum(D: &[i32], col_len: usize, start_i: usize, start_j: usize, block_width: usize, block_height: usize) -> i32 {
    let mut sum = 0;
    for i in start_i..start_i + block_height {
        for j in start_j..start_j + block_width {
            sum += D[i + j * col_len];
        }
    }
    sum
}

#[allow(dead_code)]
#[allow(non_snake_case)]
fn calc_diag(q: &[u8], r: &[u8], D: &mut [i32], R: &mut [i32], C: &mut [i32], start_i: usize, start_j: usize, block_size: usize, gap_open: i32, gap_extend: i32) -> i32 {
    let idx = |i: usize, j: usize| { i + j * (q.len() + 1 + block_size) };
    let mut max = i32::MIN;

    for off in 0..block_size {
        let i = start_i + block_size - 1 - off;
        let j = start_j + off;

        if D[idx(i, j)] != i32::MIN {
            max = cmp::max(max, D[idx(i, j)]);
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
                D[idx(i - 1, j - 1)].saturating_add(if q[i - 1] == r[j - 1] { 1 } else { -1 })
            },
            cmp::max(R[idx(i, j)], C[idx(i, j)])
        );
        max = cmp::max(max, D[idx(i, j)]);
    }

    max
}

#[allow(dead_code)]
#[allow(non_snake_case)]
fn calc_block(q: &[u8], r: &[u8], D: &mut [i32], R: &mut [i32], C: &mut [i32], start_i: usize, start_j: usize, block_width: usize, block_height: usize, block_size: usize, gap_open: i32, gap_extend: i32) -> i32 {
    let idx = |i: usize, j: usize| { i + j * (q.len() + 1 + block_size) };
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
                    D[idx(i - 1, j - 1)].saturating_add(if q[i - 1] == r[j - 1] { 1 } else { -1 })
                },
                cmp::max(R[idx(i, j)], C[idx(i, j)])
            );
            max = cmp::max(max, D[idx(i, j)]);
        }
    }

    max
}
