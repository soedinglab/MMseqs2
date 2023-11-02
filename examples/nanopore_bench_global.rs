#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use parasailors::{Matrix, *};

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use rust_wfa2::aligner::*;

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use edlib_rs::edlibrs::*;

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use ksw2_sys::*;

use block_aligner::percent_len;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use std::hint::black_box;

fn get_data(file_name: &str) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut res = vec![];

    let reader = BufReader::new(File::open(file_name).unwrap());
    let all_lines = reader.lines().collect::<Vec<_>>();

    for lines in all_lines.chunks(2) {
        let r = lines[0].as_ref().unwrap().to_ascii_uppercase();
        let q = lines[1].as_ref().unwrap().to_ascii_uppercase();
        let r = r.as_bytes().to_owned();
        let q = q.as_bytes().to_owned();
        res.push((q, r));
    }

    res
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_parasailors(file: &str) -> f64 {
    let file_data = get_data(file);
    let matrix = Matrix::create("ACGNT", 2, -4);
    let data = file_data
        .iter()
        .map(|(q, r)| (parasailors::Profile::new(q, &matrix), r.to_owned()))
        .collect::<Vec<(parasailors::Profile, Vec<u8>)>>();

    let mut total_time = 0f64;
    let mut temp = 0i32;
    for (p, r) in &data {
        let start = Instant::now();
        let res = global_alignment_score(p, r, 6, 2);
        total_time += start.elapsed().as_secs_f64();
        temp = temp.wrapping_add(res);
    }
    black_box(temp);
    total_time
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_wfa2(file: &str, use_heuristic: bool) -> f64 {
    let data = get_data(file);

    let mut wfa = WFAlignerGapAffine::new(4, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh);

    let mut total_time = 0f64;
    let mut temp = 0i32;
    for (q, r) in &data {
        if use_heuristic {
            wfa.set_heuristic(Heuristic::WFadaptive(10, percent_len(q.len().max(r.len()), 0.01) as i32, 1));
        } else {
            wfa.set_heuristic(Heuristic::None);
        }
        let start = Instant::now();
        wfa.align_end_to_end(&q, &r);
        total_time += start.elapsed().as_secs_f64();
        temp = temp.wrapping_add(wfa.score());
        temp = temp.wrapping_add(wfa.cigar().len() as i32);
    }
    black_box(temp);
    total_time
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_edlib(file: &str) -> f64 {
    let data = get_data(file);

    let mut total_time = 0f64;
    let mut temp = 0i32;
    for (q, r) in &data {
        let mut config = EdlibAlignConfigRs::default();
        config.task = EdlibAlignTaskRs::EDLIB_TASK_PATH;
        let start = Instant::now();
        let res = edlibAlignRs(&q, &r, &config);
        total_time += start.elapsed().as_secs_f64();
        temp = temp.wrapping_add(res.editDistance);
        temp = temp.wrapping_add(res.alignment.unwrap().len() as i32);
    }
    black_box(temp);
    total_time
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_ksw2(file: &str, band_width_percent: f32) -> f64 {
    let lut = {
        let mut l = [0u8; 128];
        l[b'A' as usize] = 0;
        l[b'C' as usize] = 1;
        l[b'G' as usize] = 2;
        l[b'T' as usize] = 3;
        l[b'N' as usize] = 4;
        l
    };
    let file_data = get_data(file);
    let data = file_data
        .iter()
        .map(|(q, r)| (q.iter().map(|&c| lut[c as usize]).collect(), r.iter().map(|&c| lut[c as usize]).collect()))
        .collect::<Vec<(Vec<u8>, Vec<u8>)>>();
    let matrix = {
        let mut m = [0i8; 5 * 5];
        m[0] = 2;
        m[1] = -4;
        m
    };
    let mut res: ksw_extz_t = unsafe { std::mem::zeroed() };

    let mut total_time = 0f64;
    let mut temp = 0i32;
    for (q, r) in &data {
        let band_width = percent_len(q.len().max(r.len()), band_width_percent) as i32;
        let start = Instant::now();
        unsafe {
            ksw_extz2_sse(std::ptr::null_mut(), q.len() as i32, q.as_ptr(), r.len() as i32, r.as_ptr(), 5, matrix.as_ptr(), 4, 2, band_width, -1, 0, 0 /* score only = 1 */, &mut res);
        }
        total_time += start.elapsed().as_secs_f64();
        temp = temp.wrapping_add(res.score);
        temp = temp.wrapping_add(res.n_cigar);
    }
    black_box(temp);
    total_time
}

fn bench_ours(file: &str, trace: bool, max_size: usize, block_grow: bool) -> f64 {
    let file_data = get_data(file);
    let data = file_data
        .iter()
        .map(|(q, r)| (PaddedBytes::from_bytes::<NucMatrix>(q, max_size), PaddedBytes::from_bytes::<NucMatrix>(r, max_size)))
        .collect::<Vec<(PaddedBytes, PaddedBytes)>>();
    let bench_gaps = Gaps { open: -6, extend: -2 };
    let matrix = NucMatrix::new_simple(2, -4);

    let mut total_time = 0f64;
    let mut temp = 0i32;
    for (q, r) in &data {
        let max_len = q.len().max(r.len());
        let max_percent = if block_grow { 0.1 } else { 0.01 };

        if trace {
            let mut a = Block::<true, false>::new(q.len(), r.len(), max_size);
            let mut cigar = Cigar::new(q.len(), r.len());
            let start = Instant::now();
            a.align(&q, &r, &matrix, bench_gaps, percent_len(max_len, 0.01)..=percent_len(max_len, max_percent), 0);
            a.trace().cigar(a.res().query_idx, a.res().reference_idx, &mut cigar);
            total_time += start.elapsed().as_secs_f64();
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
            temp = temp.wrapping_add(cigar.len() as i32); // prevent optimizations
        } else {
            let mut a = Block::<false, false>::new(q.len(), r.len(), max_size);
            let start = Instant::now();
            a.align(&q, &r, &matrix, bench_gaps, percent_len(max_len, 0.01)..=percent_len(max_len, max_percent), 0);
            total_time += start.elapsed().as_secs_f64();
            temp = temp.wrapping_add(a.res().score); // prevent optimizations
        }
    }
    black_box(temp);
    total_time
}

fn main() {
    let files = ["data/real.illumina.b10M.txt", "data/real.ont.b10M.txt", "data/seq_pairs.10kbps.5000.txt", "data/seq_pairs.50kbps.10000.txt"];
    let names = ["illumina", "nanopore 1kbp", "nanopore <10kbp", "nanopore <50kbp"];
    let max_sizes = [[32, 32], [32, 128], [128, 1024], [512, 8192]];
    let band_widths = [0.01, 0.1];
    let run_parasail_arr = [true, true, true, false];

    println!("# time (s)");
    println!("dataset, algorithm, time");

    for (((file, name), max_size), &run_parasail) in files.iter().zip(&names).zip(&max_sizes).zip(&run_parasail_arr) {
        for (&s, &g) in max_size.iter().zip(&[false, true]) {
            let t = bench_ours(file, true, s, g);
            println!("{}, ours ({}), {}", name, if g { "1%-10%" } else { "1%-1%" }, t);
        }

        #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
        {
            let t = bench_edlib(file);
            println!("{}, edlib, {}", name, t);
        }

        #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
        {
            for &b in &band_widths {
                let t = bench_ksw2(file, b);
                println!("{}, ksw_extz2_sse ({}%), {}", name, (b * 100.0).round() as usize, t);
            }
        }

        #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
        {
            let t = bench_wfa2(file, false);
            println!("{}, wfa2, {}", name, t);

            let t = bench_wfa2(file, true);
            println!("{}, wfa2 adaptive, {}", name, t);
        }

        if run_parasail {
            #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
            {
                let t = bench_parasailors(file);
                println!("{}, parasail, {}", name, t);
            }
        }
    }
}
