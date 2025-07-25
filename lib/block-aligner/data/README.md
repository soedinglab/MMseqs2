# Data
Many testing and benchmark programs require large files of sequence data
that should be placed in this directory.

Below are instructions for how to download the necessary data. Make sure
you are in this directory (`cd data`).

## 25kbp Nanopore data
This data is from the difference recurrence [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2014-8)
by Suzuki and Kasahara.

1. `curl -OL https://github.com/Daniel-Liu-c0deb0t/diff-bench-paper/releases/download/v1.0/sequences.txt.gz`
2. `gunzip sequences.txt.gz`

Since these reads are filtered to only have gaps smaller than 20bp, it is not representative of typical reads. Therefore,
this dataset will be rarely used.

## \<10kbp and \<50kbp Nanopore data
This data is from the BiWFA [repository](https://github.com/smarco/BiWFA-paper/tree/main/evaluation/data)
and reformatted.

1. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/seq_pairs.10kbps.5000.txt.gz`
2. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/seq_pairs.50kbps.10000.txt.gz`
3. `gunzip seq_pairs.10kbps.5000.txt.gz`
4. `gunzip seq_pairs.50kbps.10000.txt.gz`

These files contain pairs of reads that are alignable.

## Illumina and 1kbp Nanopore data
This data is from the Wavefront Aligner [paper](https://academic.oup.com/bioinformatics/article/37/4/456/5904262).

1. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/real.illumina.b10M.txt.gz`
2. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/real.ont.b10M.txt.gz`
3. `gunzip real.illumina.b10M.txt.gz`
4. `gunzip real.ont.b10M.txt.gz`

The Illumina, 1kbp Nanopore, and 25kbp Nanopore datasets are just a list of reads, where every two reads
form a pair that is alignable.

## Uniclust30 data
This data is generated with [mmseqs2](https://github.com/soedinglab/MMseqs2)
and the [Uniclust30](https://uniclust.mmseqs.com/) dataset.
Two datasets with two different coverages percentages are used: `0.8`
(default in `mmseqs2`) and `0.95`. Using a higher coverage helps gather
sequences that are "globally alignable", as `mmseqs2` uses local alignment.
The dataset with the lower coverage percent is expected to be more challenging.

Scripts for generating the data: [`0.8` coverage](uc30_pairwise_aln.sh)
and [`0.95` coverage](uc30_0.95_pairwise_aln.sh).

1. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/uc30.tar.gz`
2. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/uc30_0.95.tar.gz`
3. `tar -xvf uc30.tar.gz`
4. `tar -xvf uc30_0.95.tar.gz`

## SCOP PSSM data
This data is generated with `mmseqs2` and the [SCOPe](https://scop.berkeley.edu/astral/ver=2.01) dataset.
This data is used for aligning sequences to profiles (position-specific scoring matrices) of protein domains.

1. `mkdir scop && cd scop`
2. `curl -OL https://github.com/Daniel-Liu-c0deb0t/block-aligner/releases/download/v0.0.0/scop.tar.gz`
3. `tar -xvf scop.tar.gz`
