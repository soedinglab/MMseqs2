# EM Abundance Notes 04/13

This document summarizes the outputs implemented in `src/util/EM_abundnace.cpp`.

## Overview

`mmseqs abundance` reads the reclassified alignment DB produced by `mmseqs reclassify` and produces two possible outputs:

- Per-target abundance table (default).
- Kraken-style report when `--taxonomy 1` is set.

The input alignment DB must include the posterior probability for each hit. The current implementation reads the posterior from the `seqId` field, or from an extra trailing column if present.

## Example usage

- mmseqs abundance queryDB targetDB newDB abundance.tsv
- mmseqs abundance queryDB targetDB newDB abundance.report --taxonomy 1

Notes:
- `newDB` should be created by `mmseqs reclassify`.
- When `--taxonomy 1` is used, `targetDB_mapping` and `targetDB_taxonomy` are required.

## Per-target abundance table

For each target, the command reports:
- target key and identifier
- abundance percentage
- coverage confidence
- drop flag based on low-abundance filtering
- mapped length and target length

This table is written to the output path (typically `.tsv`).

Drop handling:
- The per-target table includes all targets and marks low-abundance filtered ones in `Drop(y/n)`.
- When `--taxonomy 1` is used, dropped targets are removed before building Kraken/Bracken reports.

## Kraken-style report

When `--taxonomy 1` is set, the output is a Kraken-style report with fields:

- percent of reads in the clade
- clade read count
- direct (taxon) read count
- rank code
- taxid
- name

The report is written to the output path (typically `.report`).

Notes on Kraken compatibility:
- The report format mirrors Kraken's `--report` layout, but values are not guaranteed to be identical to Kraken outputs.
- Counts are derived from EM abundance by converting percent to expected reads (with rounding), then aggregating clade/direct counts via the taxonomy tree.
- Ordering and handling of missing or unclassified taxa can differ from Kraken's implementation.

## Abundance from posterior

For each target $t$, abundance is computed from posteriors as:

$$
\hat{A}_t = \frac{1}{|Q|} \sum_q R_{qt}
$$

`R_{qt}` is the posterior for hit $(q,t)$ in `newDB`. Abundances are then converted to percentages and filtered by a low-abundance tail cutoff.
