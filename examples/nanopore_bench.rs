#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use parasailors::{Matrix, *};

use block_aligner::scan_block::*;
use block_aligner::scores::*;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::{Instant, Duration};
use std::hint::black_box;
use std::iter;

use simulate_seqs::*;

static FILE_NAME: &str = "data/sequences.txt";
const ITER: usize = 10000;
const LEN: usize = 10000;
const K: usize = 1000;

fn get_data(file_name: Option<&str>) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut rng = StdRng::seed_from_u64(1234);

    if let Some(file_name) = file_name {
        let mut res = vec![];

        let reader = BufReader::new(File::open(file_name).unwrap());
        let all_lines = reader.lines().collect::<Vec<_>>();

        for lines in all_lines.chunks(2) {
            let r = lines[0].as_ref().unwrap().to_ascii_uppercase();
            let q = lines[1].as_ref().unwrap().to_ascii_uppercase();
            let mut r = r.as_bytes().to_owned();
            let mut q = q.as_bytes().to_owned();
            let extend_r = rand_str(100, &NUC, &mut rng);
            let extend_q = rand_str(100, &NUC, &mut rng);
            r.extend_from_slice(&extend_r);
            q.extend_from_slice(&extend_q);
            res.push((q, r));
        }

        res
    } else {
        let mut r = rand_str(LEN, &NUC, &mut rng);
        let mut q = rand_mutate(&r, K, &NUC, &mut rng);
        let extend_r = rand_str(500, &NUC, &mut rng);
        let extend_q = rand_str(500, &NUC, &mut rng);
        r.extend_from_slice(&extend_r);
        q.extend_from_slice(&extend_q);
        black_box(iter::repeat_with(|| (q.clone(), r.clone())).take(ITER).collect())
    }
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
#[allow(dead_code)]
fn bench_parasailors_nuc_core(file: bool, _trace: bool, _max_size: usize) -> (i32, Duration) {
    let file_data = get_data(if file { Some(&FILE_NAME) } else { None });
    let matrix = Matrix::new(MatrixType::IdentityWithPenalty);
    let data = file_data
        .iter()
        .map(|(q, r)| (parasailors::Profile::new(q, &matrix), r.to_owned()))
        .collect::<Vec<(parasailors::Profile, Vec<u8>)>>();

    let start = Instant::now();
    let mut temp = 0i32;
    for (p, r) in &data {
        temp = temp.wrapping_add(global_alignment_score(p, r, 2, 1));
    }
    (temp, start.elapsed())
}

fn bench_scan_nuc_core(_file: bool, trace: bool, max_size: usize) -> (i32, Duration) {
    let file_data = get_data(None);
    let x_drop = 100;
    let matrix = NucMatrix::new_simple(2, -3);
    let data = file_data
        .iter()
        .map(|(q, r)| (PaddedBytes::from_bytes::<NucMatrix>(q, 2048), PaddedBytes::from_bytes::<NucMatrix>(r, 2048)))
        .collect::<Vec<(PaddedBytes, PaddedBytes)>>();
    let bench_gaps = Gaps { open: -5, extend: -1 };

    let start = Instant::now();
    let mut temp = 0i32;
    for (q, r) in &data {
        if trace {
            let mut a = Block::<true, true>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &matrix, bench_gaps, 32..=max_size, x_drop);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        } else {
            let mut a = Block::<false, true>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &matrix, bench_gaps, 32..=max_size, x_drop);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        }
    }
    (temp, start.elapsed())
}

fn bench_scan_nuc_file(_file: bool, trace: bool, max_size: usize) -> (i32, Duration) {
    let file_data = get_data(Some(&FILE_NAME));
    let x_drop = 50;
    let data = file_data
        .iter()
        .map(|(q, r)| (PaddedBytes::from_bytes::<NucMatrix>(q, 2048), PaddedBytes::from_bytes::<NucMatrix>(r, 2048)))
        .collect::<Vec<(PaddedBytes, PaddedBytes)>>();
    let bench_gaps = Gaps { open: -2, extend: -1 };

    let start = Instant::now();
    let mut temp = 0i32;
    for (q, r) in &data {
        if trace {
            let mut a = Block::<true, true>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &NW1, bench_gaps, 32..=max_size, x_drop);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        } else {
            let mut a = Block::<false, true>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &NW1, bench_gaps, 32..=max_size, x_drop);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        }
    }
    (temp, start.elapsed())
}

fn time(f: fn(bool, bool, usize) -> (i32, Duration), file: bool, trace: bool, max_size: usize) -> Duration {
    let (temp, duration) = f(file, trace, max_size);
    black_box(temp);
    duration
}

fn main() {
    for _i in 0..3 {
        let _d = time(bench_scan_nuc_file, true, false, 32);
    }

    println!("# time (s)");
    println!("algorithm, dataset, time");

    let d = time(bench_scan_nuc_file, true, false, 32);
    let nanopore_time = d.as_secs_f64();
    println!("ours (no trace 32-32), nanopore 25kbp, {}", nanopore_time);
    let d = time(bench_scan_nuc_core, false, false, 32);
    let random_time = d.as_secs_f64();
    println!("ours (no trace 32-32), random, {}", random_time);

    let d = time(bench_scan_nuc_file, true, true, 32);
    let nanopore_time = d.as_secs_f64();
    println!("ours (trace 32-32), nanopore 25kbp, {}", nanopore_time);
    let d = time(bench_scan_nuc_core, false, true, 32);
    let random_time = d.as_secs_f64();
    println!("ours (trace 32-32), random, {}", random_time);

    let d = time(bench_scan_nuc_file, true, true, 64);
    let nanopore_time = d.as_secs_f64();
    println!("ours (trace 32-64), nanopore 25kbp, {}", nanopore_time);
    let d = time(bench_scan_nuc_core, false, true, 64);
    let random_time = d.as_secs_f64();
    println!("ours (trace 32-64), random, {}", random_time);

    /*#[cfg(not(target_arch = "wasm32"))]
    {
        let d = time(bench_parasailors_nuc_core, true);
        let nanopore_time = d.as_secs_f64();
        println!("parasail, nanopore 25kbp, {}", nanopore_time);
        let d = time(bench_parasailors_nuc_core, false);
        let random_time = d.as_secs_f64();
        println!("parasail, random, {}", random_time);
    }*/
}
