#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
use parasailors;

use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::usize;
use std::time::{Instant, Duration};
use std::hint::black_box;

fn bench_ours(pairs: &[(AAProfile, PaddedBytes)], min_size: usize, max_size: usize) -> Duration {
    let start = Instant::now();

    for (r, q) in pairs {
        let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), max_size);
        block_aligner.align_profile(q, r, min_size..=max_size, 0);
        let mut cigar = Cigar::new(q.len(), r.len());
        block_aligner.trace().cigar(q.len(), r.len(), &mut cigar);
        black_box(block_aligner.res());
        black_box(cigar);
    }

    start.elapsed()
}

#[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
fn bench_parasail(pairs: &[(Vec<u8>, Vec<u8>)], gap_open: i8, gap_extend: i8) -> Duration {
    let matrix = parasailors::Matrix::new(parasailors::MatrixType::Blosum62);

    let start = Instant::now();

    for (r, q) in pairs {
        let p = parasailors::Profile::new(&q, &matrix);
        black_box(parasailors::global_alignment_score(&p, r, -(gap_open + gap_extend) as i32, -gap_extend as i32));
    }

    start.elapsed()
}

static MAP: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

fn get_pairs(file_name: &str, padding: usize, gap_open: i8, gap_extend: i8) -> (Vec<(AAProfile, PaddedBytes)>, Vec<(Vec<u8>, Vec<u8>)>) {
    let mut reader = BufReader::new(File::open(file_name).unwrap());
    let mut seq_string = String::new();
    let mut pssm_string = String::new();
    let mut pssm_pairs = Vec::new();
    let mut cns_pairs = Vec::new();

    loop {
        seq_string.clear();
        let len = reader.read_line(&mut seq_string).unwrap();
        if len == 0 {
            break;
        }
        let seq = seq_string.trim_end();
        pssm_string.clear();
        reader.read_line(&mut pssm_string).unwrap();
        let pssm = pssm_string.trim_end();
        let len = pssm.len() - 1;
        let cns = pssm[1..].as_bytes().to_owned();
        let mut r = AAProfile::new(len, padding, gap_extend);
        let q = PaddedBytes::from_str::<AAMatrix>(&seq[1..], padding);

        for i in 0..len + 1 {
            pssm_string.clear();
            reader.read_line(&mut pssm_string).unwrap();
            let pssm = pssm_string.trim_end();
            if i == 0 {
                continue;
            }

            for (j, s) in pssm.split_whitespace().skip(2).enumerate() {
                let c = MAP[j];
                let s = s.parse::<i8>().unwrap();
                r.set(i, c, s);
            }

            r.set_gap_open_C(i, gap_open);
            r.set_gap_close_C(i, 0);
            r.set_gap_open_R(i, gap_open);
        }

        pssm_pairs.push((r, q));
        cns_pairs.push((cns, seq[1..].as_bytes().to_owned()));
    }

    (pssm_pairs, cns_pairs)
}

fn main() {
    let file_name = "data/scop/pairs.pssm";
    let min_sizes = [32, 32, 32, 128];
    let max_sizes = [32, 64, 128, 128];
    let gap_open = -10;
    let gap_extend = -1;

    let (pssm_pairs, cns_pairs) = get_pairs(file_name, 2048, gap_open, gap_extend);

    println!("size, time");

    bench_ours(&pssm_pairs, min_sizes[0], max_sizes[0]);

    for (&min_size, &max_size) in min_sizes.iter().zip(&max_sizes) {
        let duration = bench_ours(&pssm_pairs, min_size, max_size);
        println!("{}-{}, {}", min_size, max_size, duration.as_secs_f64());
    }

    #[cfg(not(any(feature = "simd_wasm", feature = "simd_neon")))]
    {
        let duration = bench_parasail(&cns_pairs, gap_open, gap_extend);
        println!("parasail, {}", duration.as_secs_f64());
    }

    println!("# Done!");
}
