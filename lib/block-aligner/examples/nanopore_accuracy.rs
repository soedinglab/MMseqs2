#[cfg(not(feature = "simd_avx2"))]
fn main() {}

#[cfg(feature = "simd_avx2")]
fn test(file_name: &str, max_size: usize, name: &str, verbose: bool, writer: &mut impl std::io::Write) -> (usize, usize, usize, f64, usize, f64) {
    use parasailors::{Matrix, *};

    use rust_wfa2::aligner::*;

    use bio::alignment::distance::simd::levenshtein;

    use block_aligner::percent_len;
    use block_aligner::scan_block::*;
    use block_aligner::scores::*;

    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let mut wrong = 0usize;
    let mut min_size_wrong = 0usize;
    let mut wfa_wrong = 0usize;
    let mut wrong_avg = 0f64;
    let mut count = 0usize;
    let mut seq_id_avg = 0f64;
    let reader = BufReader::new(File::open(file_name).unwrap());
    let all_lines = reader.lines().collect::<Vec<_>>();

    for lines in all_lines.chunks(2) {
        let r = lines[0].as_ref().unwrap().to_ascii_uppercase();
        let q = lines[1].as_ref().unwrap().to_ascii_uppercase();

        let correct_score;

        if r.len().max(q.len()) < 15000 {
            // parasail
            let matrix = Matrix::create("ACGNT", 2, -4);
            let profile = parasailors::Profile::new(q.as_bytes(), &matrix);
            let parasail_score = global_alignment_score(&profile, r.as_bytes(), 6, 2);
            correct_score = parasail_score;
        } else {
            // parasail is not accurate enough, so use block aligner with large fixed block size
            let len = 8192;
            let r_padded = PaddedBytes::from_bytes::<NucMatrix>(r.as_bytes(), len);
            let q_padded = PaddedBytes::from_bytes::<NucMatrix>(q.as_bytes(), len);
            let run_gaps = Gaps { open: -6, extend: -2 };
            let matrix = NucMatrix::new_simple(2, -4);
            let mut block_aligner = Block::<false, false>::new(q.len(), r.len(), len);
            block_aligner.align(&q_padded, &r_padded, &matrix, run_gaps, len..=len, 0);
            let scan_score = block_aligner.res().score;
            correct_score = scan_score;
        }

        let r_padded = PaddedBytes::from_bytes::<NucMatrix>(r.as_bytes(), max_size);
        let q_padded = PaddedBytes::from_bytes::<NucMatrix>(q.as_bytes(), max_size);
        let run_gaps = Gaps { open: -6, extend: -2 };
        let matrix = NucMatrix::new_simple(2, -4);

        // ours
        let mut block_aligner = Block::<false, false>::new(q.len(), r.len(), max_size);
        let max_len = q.len().max(r.len());
        block_aligner.align(&q_padded, &r_padded, &matrix, run_gaps, percent_len(max_len, 0.01)..=percent_len(max_len, 0.1), 0);
        let scan_score = block_aligner.res().score;

        write!(
            writer,
            "{}, {}, {}, {}, {}\n",
            name,
            q.len(),
            r.len(),
            scan_score,
            correct_score
        ).unwrap();

        if correct_score != scan_score {
            wrong += 1;
            wrong_avg += ((correct_score - scan_score) as f64) / (correct_score as f64);

            if verbose {
                let edit_dist = levenshtein(q.as_bytes(), r.as_bytes());
                println!(
                    "parasail: {}, ours: {}, edit dist: {}\nq (len = {}): {}\nr (len = {}): {}",
                    correct_score,
                    scan_score,
                    edit_dist,
                    q.len(),
                    q,
                    r.len(),
                    r
                );
            }
        }

        block_aligner.align(&q_padded, &r_padded, &matrix, run_gaps, percent_len(max_len, 0.01)..=percent_len(max_len, 0.01), 0);
        let min_size_score = block_aligner.res().score;
        if min_size_score != correct_score {
            min_size_wrong += 1;
        }

        let wfa_adaptive_score = {
            let mut wfa = WFAlignerGapAffine::new(4, 4, 2, AlignmentScope::Score, MemoryModel::MemoryHigh);
            wfa.set_heuristic(Heuristic::WFadaptive(10, 50, 1));
            wfa.align_end_to_end(q.as_bytes(), r.as_bytes());
            wfa.score()
        };
        let wfa_score = {
            let mut wfa = WFAlignerGapAffine::new(4, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh);
            wfa.set_heuristic(Heuristic::None);
            wfa.align_end_to_end(q.as_bytes(), r.as_bytes());
            let cigar = wfa.cigar();
            let matches = cigar.bytes().filter(|&c| c == b'M').count();
            let seq_id = (matches as f64) / (cigar.len() as f64);
            seq_id_avg += seq_id;
            wfa.score()
        };
        if wfa_adaptive_score != wfa_score {
            wfa_wrong += 1;
        }

        count += 1;
    }

    (wrong, min_size_wrong, wfa_wrong, wrong_avg / (wrong as f64), count, seq_id_avg / (count as f64))
}

#[cfg(feature = "simd_avx2")]
fn main() {
    use std::env;
    use std::fs::File;
    use std::io::{Write, BufWriter};

    let arg1 = env::args().skip(1).next();
    let verbose = arg1.is_some() && arg1.unwrap() == "-v";
    let paths = ["data/real.illumina.b10M.txt", "data/real.ont.b10M.txt", "data/seq_pairs.10kbps.5000.txt", "data/seq_pairs.50kbps.10000.txt"];
    let names = ["illumina", "nanopore 1kbp", "nanopore <10kbp", "nanopore <50kbp"];
    let max_size = [32, 128, 1024, 8192];

    let out_file_name = "data/nanopore_accuracy.csv";
    let mut writer = BufWriter::new(File::create(out_file_name).unwrap());
    write!(writer, "dataset, query len, reference len, pred score, true score\n").unwrap();

    println!("\ndataset, total, wrong, wrong % error, min size wrong, wfa wrong");

    for ((path, name), &max_size) in paths.iter().zip(&names).zip(&max_size) {
        let (wrong, min_size_wrong, wfa_wrong, wrong_avg, count, seq_id_avg) = test(path, max_size, name, verbose, &mut writer);
        println!("\n{}, {}, {}, {}, {}, {}", name, count, wrong, wrong_avg, min_size_wrong, wfa_wrong);
        println!("# {} seq id avg: {}", name, seq_id_avg);
    }

    println!("# Done!");
}
