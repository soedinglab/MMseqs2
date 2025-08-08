use bio::alignment::pairwise::*;
use bio::alignment::{Alignment, AlignmentOperation};
use bio::scores::blosum62;

use block_aligner::scan_block::*;
use block_aligner::scores::*;
use block_aligner::cigar::*;

use std::{env, cmp};
use std::fs::File;
use std::io::{BufRead, Write, BufReader, BufWriter};
use std::usize;

fn test(file_name: &str, min_size: usize, max_size: usize, string: &str, verbose: bool, wrong: &mut [usize], wrong_avg: &mut [f64], count: &mut [usize], writer: &mut impl Write) -> (f64, usize, usize, f64) {
    let reader = BufReader::new(File::open(file_name).unwrap());
    let mut length_sum = 0f64;
    let mut length_min = usize::MAX;
    let mut length_max = usize::MIN;
    let mut dp_fraction = 0f64;

    for line in reader.lines() {
        let line = line.unwrap();
        let mut last_two = line.split_ascii_whitespace().rev().take(2);
        let r = last_two.next().unwrap().to_ascii_uppercase();
        let q = last_two.next().unwrap().to_ascii_uppercase();

        // rust-bio
        let mut bio_aligner = Aligner::with_capacity(q.len(), r.len(), -10, -1, &blosum62);
        let bio_alignment = bio_aligner.global(q.as_bytes(), r.as_bytes());
        let bio_score = bio_alignment.score;
        let seq_identity = seq_id(&bio_alignment);
        let id_idx = cmp::min((seq_identity * 10.0) as usize, 9);
        let indels = indels(&bio_alignment, cmp::max(q.len(), r.len()));

        let r_padded = PaddedBytes::from_bytes::<AAMatrix>(r.as_bytes(), 2048);
        let q_padded = PaddedBytes::from_bytes::<AAMatrix>(q.as_bytes(), 2048);
        let run_gaps = Gaps { open: -11, extend: -1 };

        // ours
        let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), max_size);
        block_aligner.align(&q_padded, &r_padded, &BLOSUM62, run_gaps, min_size..=max_size, 0);
        let scan_res = block_aligner.res();
        let scan_score = scan_res.score;

        write!(
            writer,
            "{}, {}-{}, {}, {}, {}, {}, {}\n",
            string,
            min_size,
            max_size,
            q.len(),
            r.len(),
            seq_identity,
            scan_score,
            bio_score
        ).unwrap();

        if bio_score != scan_score {
            wrong[id_idx] += 1;
            wrong_avg[id_idx] += ((bio_score - scan_score) as f64) / (bio_score as f64);

            if verbose {
                let mut cigar = Cigar::new(scan_res.query_idx, scan_res.reference_idx);
                block_aligner.trace().cigar(scan_res.query_idx, scan_res.reference_idx, &mut cigar);
                let (a_pretty, b_pretty) = cigar.format(q.as_bytes(), r.as_bytes());
                println!(
                    "seq id: {}, max indel len: {}, bio: {}, ours: {}\nq (len = {}): {}\nr (len = {}): {}\nbio pretty:\n{}\nours pretty:\n{}\n{}",
                    seq_identity,
                    indels,
                    bio_score,
                    scan_score,
                    q.len(),
                    q,
                    r.len(),
                    r,
                    bio_alignment.pretty(q.as_bytes(), r.as_bytes()),
                    a_pretty,
                    b_pretty
                );
            }
        }

        count[id_idx] += 1;
        length_sum += (q.len() + r.len()) as f64;
        length_min = cmp::min(length_min, cmp::min(q.len(), r.len()));
        length_max = cmp::max(length_max, cmp::max(q.len(), r.len()));

        let computed = block_aligner.trace().blocks().iter().map(|b| (b.width as f64) * (b.height as f64)).sum::<f64>();
        dp_fraction += computed / (((q.len() + 1) as f64) * ((r.len() + 1) as f64));
    }

    (length_sum, length_min, length_max, dp_fraction)
}

fn indels(a: &Alignment, len: usize) -> f64 {
    let mut indels = 0;

    for &op in &a.operations {
        if op == AlignmentOperation::Ins
            || op == AlignmentOperation::Del {
            indels += 1;
        }
    }
    (indels as f64) / (len as f64)
}

// BLAST sequence identity
fn seq_id(a: &Alignment) -> f64 {
    let mut matches = 0;

    for &op in &a.operations {
        if op == AlignmentOperation::Match {
            matches += 1;
        }
    }

    (matches as f64) / (a.operations.len() as f64)
}

fn main() {
    let arg1 = env::args().skip(1).next();
    let verbose = arg1.is_some() && arg1.unwrap() == "-v";
    let file_names_arr = [
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
        ],
        [
            "data/uc30_30_40.m8",
            "data/uc30_40_50.m8",
            "data/uc30_50_60.m8",
            "data/uc30_60_70.m8",
            "data/uc30_70_80.m8",
            "data/uc30_80_90.m8",
            "data/uc30_90_100.m8"
        ]
    ];
    let strings = [/*"merged_clu_aln", */"uc30_0.95", "uc30"];
    let min_sizes = [32, 32, 256];
    let max_sizes = [32, 256, 256];

    let out_file_name = "data/uc_accuracy.csv";
    let mut writer = BufWriter::new(File::create(out_file_name).unwrap());
    write!(writer, "dataset, size, query len, reference len, seq id, pred score, true score\n").unwrap();

    println!("# seq identity is lower bound (inclusive)");
    println!("dataset, size, seq identity, count, wrong, wrong % error");

    for (file_names, string) in file_names_arr.iter().zip(&strings) {
        for (&min_size, &max_size) in min_sizes.iter().zip(&max_sizes) {
            let mut wrong = [0usize; 10];
            let mut wrong_avg = [0f64; 10];
            let mut count = [0usize; 10];
            let mut length_avg = 0f64;
            let mut length_min = usize::MAX;
            let mut length_max = usize::MIN;
            let mut dp_fraction = 0f64;

            for file_name in file_names {
                let (len_sum, len_min, len_max, dp_fract) = test(file_name, min_size, max_size, string, verbose, &mut wrong, &mut wrong_avg, &mut count, &mut writer);
                length_avg += len_sum;
                length_min = cmp::min(length_min, len_min);
                length_max = cmp::max(length_max, len_max);
                dp_fraction += dp_fract;
            }

            length_avg /= (count.iter().sum::<usize>() * 2) as f64;
            dp_fraction /= count.iter().sum::<usize>() as f64;

            for i in 0..10 {
                println!(
                    "{}, {}-{}, {}, {}, {}, {}",
                    string,
                    min_size,
                    max_size,
                    (i as f64) / 10.0,
                    count[i],
                    wrong[i],
                    (wrong_avg[i] as f64) / (wrong[i] as f64)
                );
            }

            println!(
                "\n# total: {}, wrong: {}, wrong % error: {}, length avg: {}, length min: {}, length max: {}, dp fraction: {}\n",
                count.iter().sum::<usize>(),
                wrong.iter().sum::<usize>(),
                wrong_avg.iter().sum::<f64>() / (wrong.iter().sum::<usize>() as f64),
                length_avg,
                length_min,
                length_max,
                dp_fraction
            );
        }
    }

    println!("# Done!");
}
