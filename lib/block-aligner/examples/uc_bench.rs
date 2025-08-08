#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use parasailors::{Matrix, *};

use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::{Instant, Duration};
use std::hint::black_box;

static FILE_NAMES: [[&str; 7]; 2] = [
    [
        "data/uc30_30_40.m8",
        "data/uc30_40_50.m8",
        "data/uc30_50_60.m8",
        "data/uc30_60_70.m8",
        "data/uc30_70_80.m8",
        "data/uc30_80_90.m8",
        "data/uc30_90_100.m8"
    ],
    /*[
        "data/merged_clu_aln_30_40.m8",
        "data/merged_clu_aln_40_50.m8",
        "data/merged_clu_aln_50_60.m8",
        "data/merged_clu_aln_60_70.m8",
        "data/merged_clu_aln_70_80.m8",
        "data/merged_clu_aln_80_90.m8",
        "data/merged_clu_aln_90_100.m8"
    ],*/
    [
        "data/uc30_0.95_30_40.m8",
        "data/uc30_0.95_40_50.m8",
        "data/uc30_0.95_50_60.m8",
        "data/uc30_0.95_60_70.m8",
        "data/uc30_0.95_70_80.m8",
        "data/uc30_0.95_80_90.m8",
        "data/uc30_0.95_90_100.m8"
    ]
];

fn get_data(file_names: &[&str]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut res = vec![];

    for file_name in file_names {
        let reader = BufReader::new(File::open(file_name).unwrap());

        for line in reader.lines() {
            let line = line.unwrap();
            let mut last_two = line.split_ascii_whitespace().rev().take(2);
            let r = last_two.next().unwrap().to_ascii_uppercase();
            let q = last_two.next().unwrap().to_ascii_uppercase();

            res.push((q.as_bytes().to_owned(), r.as_bytes().to_owned()));
        }
    }

    res
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_parasailors_aa_core(idx: usize, _trace: bool, _min_size: usize, _max_size: usize) -> (i32, Duration) {
    let file_data = get_data(&FILE_NAMES[idx]);
    let matrix = Matrix::new(MatrixType::Blosum62);
    let data = file_data
        .iter()
        .map(|(q, r)| (parasailors::Profile::new(q, &matrix), r.to_owned()))
        .collect::<Vec<(parasailors::Profile, Vec<u8>)>>();

    let start = Instant::now();
    let mut temp = 0i32;
    for (p, r) in &data {
        temp = temp.wrapping_add(global_alignment_score(p, r, 11, 1));
    }
    (temp, start.elapsed())
}

fn bench_scan_aa_core(idx: usize, trace: bool, min_size: usize, max_size: usize) -> (i32, Duration) {
    let file_data = get_data(&FILE_NAMES[idx]);
    let data = file_data
        .iter()
        .map(|(q, r)| (PaddedBytes::from_bytes::<AAMatrix>(q, 2048), PaddedBytes::from_bytes::<AAMatrix>(r, 2048)))
        .collect::<Vec<(PaddedBytes, PaddedBytes)>>();
    let bench_gaps = Gaps { open: -11, extend: -1 };

    let start = Instant::now();
    let mut temp = 0i32;
    for (q, r) in &data {
        if trace {
            let mut a = Block::<true, false>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &BLOSUM62, bench_gaps, min_size..=max_size, 0);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
            let mut cigar = Cigar::new(q.len(), r.len());
            a.trace().cigar(q.len(), r.len(), &mut cigar);
            temp = temp.wrapping_add(cigar.len() as i32);
        } else {
            let mut a = Block::<false, false>::new(q.len(), r.len(), max_size);
            a.align(&q, &r, &BLOSUM62, bench_gaps, min_size..=max_size, 0);
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        }
    }
    (temp, start.elapsed())
}

fn time(f: fn(usize, bool, usize, usize) -> (i32, Duration), idx: usize, trace: bool, min_size: usize, max_size: usize) -> Duration {
    let (temp, duration) = f(idx, trace, min_size, max_size);
    black_box(temp);
    duration
}

fn main() {
    for _i in 0..2 {
        let _d = time(bench_scan_aa_core, 1, false, 32, 32);
    }

    println!("# time (s)");
    println!("algorithm, dataset, size, time");

    let d = time(bench_scan_aa_core, 0, false, 32, 32);
    let uc30_time = d.as_secs_f64();
    println!("ours (no trace), uc30, 32-32, {}", uc30_time);
    let d = time(bench_scan_aa_core, 1, false, 32, 32);
    let uc30_95_time = d.as_secs_f64();
    println!("ours (no trace), uc30 0.95, 32-32, {}", uc30_95_time);

    let d = time(bench_scan_aa_core, 0, false, 32, 256);
    let uc30_time = d.as_secs_f64();
    println!("ours (no trace), uc30, 32-256, {}", uc30_time);
    /*let d = time(bench_scan_aa_core, 1);
    println!("scan merged time (s): {}", d.as_secs_f64());*/
    let d = time(bench_scan_aa_core, 1, false, 32, 256);
    let uc30_95_time = d.as_secs_f64();
    println!("ours (no trace), uc30 0.95, 32-256, {}", uc30_95_time);

    let d = time(bench_scan_aa_core, 0, false, 256, 256);
    let uc30_time = d.as_secs_f64();
    println!("ours (no trace), uc30, 256-256, {}", uc30_time);
    let d = time(bench_scan_aa_core, 1, false, 256, 256);
    let uc30_95_time = d.as_secs_f64();
    println!("ours (no trace), uc30 0.95, 256-256, {}", uc30_95_time);

    let d = time(bench_scan_aa_core, 0, true, 32, 256);
    let uc30_time = d.as_secs_f64();
    println!("ours (trace), uc30, 32-256, {}", uc30_time);
    let d = time(bench_scan_aa_core, 1, true, 32, 256);
    let uc30_95_time = d.as_secs_f64();
    println!("ours (trace), uc30 0.95, 32-256, {}", uc30_95_time);

    #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
    {
        let d = time(bench_parasailors_aa_core, 0, false, 0, 0);
        let uc30_time = d.as_secs_f64();
        println!("parasail, uc30, full, {}", uc30_time);
        /*let d = time(bench_parasailors_aa_core, 1);
        println!("parasail merged time (s): {}", d.as_secs_f64());*/
        let d = time(bench_parasailors_aa_core, 1, false, 0, 0);
        let uc30_95_time = d.as_secs_f64();
        println!("parasail, uc30 0.95, full, {}", uc30_95_time);
    }
}
