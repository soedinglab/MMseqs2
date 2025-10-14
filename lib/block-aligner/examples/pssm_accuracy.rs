use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;

use std::fs::File;
use std::io::{BufRead, Write, BufReader, BufWriter};
use std::usize;

fn run_ours(r: &AAProfile, q: &PaddedBytes, min_size: &[usize], max_size: &[usize]) -> Vec<i32> {
    let mut scores = Vec::new();

    for i in 0..min_size.len() {
        let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), max_size[i]);
        block_aligner.align_profile(q, r, min_size[i]..=max_size[i], 0);
        let score = block_aligner.res().score;
        let mut cigar = Cigar::new(q.len(), r.len());
        block_aligner.trace().cigar(q.len(), r.len(), &mut cigar);
        //println!("{} {}", score, cigar);
        scores.push(score);
    }

    scores
}

static MAP: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

fn test(file_name: &str, out_file_name: &str, min_size: &[usize], max_size: &[usize], padding: usize, gap_open: i8, gap_extend: i8) {
    let mut reader = BufReader::new(File::open(file_name).unwrap());
    let mut writer = BufWriter::new(File::create(out_file_name).unwrap());
    let mut seq_string = String::new();
    let mut pssm_string = String::new();
    let mut matches = vec![0usize; min_size.len()];

    println!("size, correct");
    write!(writer, "size, seq len, profile len, pred score, true score\n").unwrap();

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

        let scores = run_ours(&r, &q, min_size, max_size);
        for i in 0..matches.len() {
            write!(
                writer,
                "{}-{}, {}, {}, {}, {}\n",
                min_size[i],
                max_size[i],
                q.len(),
                r.len(),
                scores[i],
                scores[scores.len() - 1]
            ).unwrap();
            matches[i] += if scores[i] == scores[scores.len() - 1] { 1 } else { 0 };
        }
    }

    for (i, m) in matches.iter().enumerate() {
        println!("{}-{}, {}", min_size[i], max_size[i], m);
    }

    println!("# compared to {}-{}", min_size[min_size.len() - 1], max_size[max_size.len() - 1]);
}

fn main() {
    let file_name = "data/scop/pairs.pssm";
    let out_file_name = "data/pssm_accuracy.csv";
    let min_sizes = [32, 32, 32, 128, 2048];
    let max_sizes = [32, 64, 128, 128, 2048];
    let gap_open = -10;
    let gap_extend = -1;

    test(file_name, out_file_name, &min_sizes, &max_sizes, 2048, gap_open, gap_extend);

    println!("# Done!");
}
