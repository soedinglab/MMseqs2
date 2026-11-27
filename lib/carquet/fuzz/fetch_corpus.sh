#!/usr/bin/env bash
#
# Fetch the Apache parquet-testing data files and use them as fuzzing seeds.
#
# These are real-world Parquet files (varied types, encodings, codecs,
# nested/nullable, multi-row-group, page-indexed, plus intentionally malformed
# inputs under bad_data/). They are a far more diverse seed source than any
# hand-written generator, so we don't commit them — fetch on demand here.
#
# The files land in fuzz/external/parquet-testing (gitignored). run_fuzzer.py's
# ensure_corpus() picks them up automatically for the reader and page_filter
# targets the next time those corpora are (re)seeded.
#
# Usage:
#   fuzz/fetch_corpus.sh            # clone (or update) into fuzz/external/
#
# Note: some files intentionally exercise edge cases or unsupported features
# and may crash a buggy build — that is the point. A crashing seed is a finding.

set -euo pipefail
cd "$(dirname "$0")"

DEST="external/parquet-testing"
REPO="https://github.com/apache/parquet-testing.git"

if [ -d "$DEST/.git" ]; then
    echo "Updating $DEST ..."
    git -C "$DEST" pull --ff-only
else
    echo "Cloning $REPO into $DEST ..."
    mkdir -p external
    rm -rf "$DEST"
    git clone --depth 1 "$REPO" "$DEST"
fi

n=$(find "$DEST" -name '*.parquet' | wc -l | tr -d ' ')
echo "Done: $n .parquet files under $DEST"
echo "Seed dirs used by run_fuzzer.py:"
echo "  data/, data/geospatial/, bad_data/  -> reader, page_filter corpora"
