#!/bin/sh -e
# Shared-filesystem parallel linclust.
#
# Runs the two linclust passes across many independent worker processes that
# coordinate only through files: no MPI, no rank argument, no node-to-node
# communication. Every worker of a stage runs a byte-identical command line, so a
# stage maps onto a Slurm array job and workers may join late, die, or be
# restarted.
#
#   RUNNER="srun -n 64" ./linclustparallel.sh input.fasta clusters.tsv tmp
#
# Two things about this script are load-bearing and should not be "simplified":
#
#   1. **The parameters are derived once and passed consistently.** The stages do
#      not re-derive them. `--cov-mode` must reach the *grouping* stage, because
#      assignGroup uses it to decide edge orientation -- stock passes it to
#      `kmermatcher`, which does extraction and grouping together, whereas here
#      those are separate commands. `--kmer-per-seq 21` must be given explicitly
#      because the standalone default derives 20, and `--min-seq-id` selects the
#      k-mer length. Mismatching any of these silently produces a different,
#      wrong clustering rather than an error.
#   2. **Pass 2 runs on a densely re-keyed representative database.** Every stage
#      assumes dense keys, so `createrepdb` re-keys rather than sub-setting, and
#      `translatecluster` maps the result back before the merge.
#
# Restart behaviour, which differs by stage and matters on a walltime-limited
# queue:
#   - The multi-worker stages are **never skipped**. Their work queue is the
#     source of truth, so re-running is how they resume; if everything is already
#     done the workers see that and exit immediately. Guarding them on an output
#     file would skip a half-finished stage.
#   - The single-node stages are not resumable, so they are guarded on their
#     output -- but they write to a temporary name and rename on success, so a
#     half-written file is never mistaken for a finished one. Redoing one costs
#     about an hour at 1e11, comfortably inside a 24 h walltime.
[ -z "$MMSEQS" ] && MMSEQS=mmseqs
[ -z "$RUNNER" ] && RUNNER=""
[ -z "$THREADS" ] && THREADS=$(nproc 2>/dev/null || echo 8)
[ -z "$MIN_SEQ_ID" ] && MIN_SEQ_ID=0.9
[ -z "$COV" ] && COV=0.8
[ -z "$COV_MODE" ] && COV_MODE=1
[ -z "$EVAL" ] && EVAL=0.001
[ -z "$SPLIT_MEMORY_LIMIT" ] && SPLIT_MEMORY_LIMIT=0
[ -z "$SCRATCH_BUDGET" ] && SCRATCH_BUDGET=0

notExists() { [ ! -f "$1" ]; }
fail() { echo "Error: $1"; exit 1; }

# Bytes currently sitting on the scratch filesystem for this run.
#
# The wave count has to be derived against what is actually free, not against the
# whole budget: pass 2 runs with all of pass 1's surviving output still on disk,
# and sizing it as though the filesystem were empty is what made a 1e11 run derive
# a single wave and then peak at ~1.9x its ceiling. Measuring beats modelling here
# -- it needs no per-stage accounting and stays right when the pipeline changes.
scratchUsed() {
    du -sb "$TMP" "$DB" 2>/dev/null | awk '{s += $1} END {print s + 0 "B"}'
}

# Deletes an intermediate whose consumers have all finished. Nothing downstream
# reopens these, and they are the bulk of peak scratch: at 100M the candidate
# edges alone are 21 GB and the pass-1 alignments 3 GB.
dropIntermediate() {
    [ -n "$KEEP_INTERMEDIATE" ] && return 0
    rm -rf "$@"
}

# Extraction waves. A wave re-extracts every k-mer but keeps only its own slice of
# partition space, so peak scratch is the whole shuffle divided by the wave count,
# paid for with that many passes over the sequences. How many are needed follows
# from --scratch-budget and is decided by the map, not here; the map writes it into
# the shuffle manifest, and wave 0 is always valid, so wave 0 runs first and the
# rest of the loop reads its count off the manifest. Each wave is reduced before
# the next is mapped, and the reduce drops the buckets it consumed -- running two
# waves' buckets at once is exactly what the budget said would not fit.
#
# The alignment does *not* run per wave. Waves partition k-mer space while edges
# are bucketed by representative, so every wave's surviving edges land in the same
# buckets and the align runs once, afterwards, over their union.
waveCount() { awk '$1 == "waveCount" { print $2 }' "$1/coord/shuffle.info"; }

mapReduceWaves() {
    # $1 sequence DB, $2 k-mer dir, $3 edge dir, $4... extra map arguments
    _db="$1"; _kmer="$2"; _edges="$3"; shift 3
    # shellcheck disable=SC2086
    _used=$(scratchUsed)
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcherparallel "$_db" "$_kmer" $KMER_COMMON "$@" --kmer-wave 0 \
        --scratch-used "$_used" \
        || fail "kmermatcherparallel died"
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmerreduceparallel "$_db" "$_kmer" "$_edges" $REDUCE_PAR --kmer-wave 0 \
        || fail "kmerreduceparallel died"
    _waves=$(waveCount "$_kmer")
    [ -z "$_waves" ] && fail "no wave count in $_kmer/coord/shuffle.info"
    _w=1
    while [ "$_w" -lt "$_waves" ]; do
        echo "--- extraction wave $_w of $_waves ---"
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcherparallel "$_db" "$_kmer" $KMER_COMMON "$@" --kmer-wave $_w \
            --scratch-used "$_used" \
            || fail "kmermatcherparallel (wave $_w) died"
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmerreduceparallel "$_db" "$_kmer" "$_edges" $REDUCE_PAR --kmer-wave $_w \
            || fail "kmerreduceparallel (wave $_w) died"
        _w=$((_w + 1))
    done
}

[ "$#" -ne 3 ] && { echo "usage: [RUNNER=\"srun -n N\"] $0 <i:fasta|sequenceDB> <o:clusterTsv> <tmpDir>"; exit 1; }
INPUT="$1"
OUT="$2"
TMP="$3"
mkdir -p "$TMP"

# Shared by both passes. --min-seq-id belongs here because it selects the k-mer
# length and alphabet; coverage does not, because this stage only extracts k-mers.
KMER_COMMON="--alph-size aa:21,nucl:5 --min-seq-id $MIN_SEQ_ID --kmer-per-seq 21 \
 --mask 0 --mask-prob 0.9 --mask-lower-case 0 --mask-n-repeat 0 -k 0 --max-seq-len 65535 \
 --hash-shift 67 --ignore-multi-kmer 0 \
 --split-memory-limit $SPLIT_MEMORY_LIMIT --scratch-budget $SCRATCH_BUDGET --threads $THREADS"
REDUCE_PAR="-c $COV --cov-mode $COV_MODE --include-adjacency 1 --num-adjacency 3 \
 --split-memory-limit $SPLIT_MEMORY_LIMIT --threads $THREADS"
ALIGN_PAR="--min-seq-id $MIN_SEQ_ID --min-aln-len 0 --seq-id-mode 0 -e $EVAL -c $COV \
 --cov-mode $COV_MODE --threads $THREADS"

# 0. Sequence database with dense, length-ranked keys. Every later stage depends
#    on that ordering: it is what makes the greedy a single forward sweep.
DB="$INPUT"
if notExists "$INPUT.dbtype"; then
    DB="$TMP/db"
    # Guarded on the finalize sentinel, not on .dbtype: createdbparallel writes
    # .dbtype before the text indices, so a death in that window leaves a database
    # that looks complete and is not. The sentinel is written last.
    if notExists "$DB.coord/finalize.done"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" createdbparallel "$INPUT" "$DB" --threads $THREADS \
            || fail "createdbparallel died"
    fi
fi

# ---- pass 1, over the whole database --------------------------------------
mapReduceWaves "$DB" "$TMP/kmer1" "$TMP/edges1" \
    --spaced-kmer-mode 0 --kmer-per-seq-scale aa:0.000,nucl:0.200

# shellcheck disable=SC2086
    $RUNNER "$MMSEQS" alignparallel "$DB" "$TMP/edges1" "$TMP/aln1" $ALIGN_PAR \
        || fail "alignparallel (pass 1) died"

# Single node from here to the end of the pass: the greedy sweep is sequential by
# necessity, which is what makes it exact.
if notExists "$TMP/clu1.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" greedycluster "$DB" "$TMP/aln1" "$TMP/clu1.tsv.tmp" --threads $THREADS \
        || fail "greedycluster (pass 1) died"
    mv -f "$TMP/clu1.tsv.tmp" "$TMP/clu1.tsv"
fi
# Both are dead once clu1.tsv exists, and together they are the largest thing pass
# 1 leaves behind for pass 2 to be sized around.
dropIntermediate "$TMP/edges1" "$TMP/aln1"

# ---- pass 2, over the representatives --------------------------------------
# Re-keyed densely rather than sub-set, because every stage addresses keys as
# array offsets. The key map translates the result back below.
if notExists "$TMP/rep.keymap"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createrepdb "$DB" "$TMP/clu1.tsv" "$TMP/rep" --threads $THREADS \
        || fail "createrepdb died"
fi

mapReduceWaves "$TMP/rep" "$TMP/kmer2" "$TMP/edges2" \
    --spaced-kmer-mode 1 --kmer-per-seq-scale aa:0.100,nucl:0.100

# The filter gate: a representative may only take a member if every sequence of
# that member's pass-1 cluster also aligns to it. Needs the original database and
# the key map, since the clustering it consults is in original keys.
# shellcheck disable=SC2086
    $RUNNER "$MMSEQS" alignparallel "$TMP/rep" "$TMP/edges2" "$TMP/aln2" $ALIGN_PAR \
        --filter-cludb-file "$TMP/clu1.tsv" --filter-seqdb-file "$DB" \
        --key-map "$TMP/rep.keymap" \
        || fail "alignparallel (pass 2) died"

if notExists "$TMP/clu2_sub.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" greedycluster "$TMP/rep" "$TMP/aln2" "$TMP/clu2_sub.tsv.tmp" --threads $THREADS \
        || fail "greedycluster (pass 2) died"
    mv -f "$TMP/clu2_sub.tsv.tmp" "$TMP/clu2_sub.tsv"
fi
dropIntermediate "$TMP/edges2" "$TMP/aln2" "$TMP/kmer2"

if notExists "$TMP/clu2.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatecluster "$TMP/clu2_sub.tsv" "$TMP/rep.keymap" "$TMP/clu2.tsv.tmp" \
        --threads $THREADS --split-memory-limit $SPLIT_MEMORY_LIMIT \
        || fail "translatecluster died"
    mv -f "$TMP/clu2.tsv.tmp" "$TMP/clu2.tsv"
fi

# ---- fold pass 2 into pass 1 ------------------------------------------------
# Keys, not accessions: every stage addresses sequences as dense keys, so this is
# the form the merge produces. It is translated below.
if notExists "$TMP/clu.keys.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" mergeclusterparallel "$DB" "$TMP/clu1.tsv" "$TMP/clu2.tsv" "$TMP/clu.keys.tsv.tmp" \
        --threads $THREADS --split-memory-limit $SPLIT_MEMORY_LIMIT \
        || fail "mergeclusterparallel died"
    mv -f "$TMP/clu.keys.tsv.tmp" "$TMP/clu.keys.tsv"
fi

# ---- back to accession space ------------------------------------------------
# Stock reaches this with createtsv, which needs a result database and a resident
# id->name table. Here it is a streaming join against the .lookup. Skipped when
# the database has none (--write-lookup 0), leaving the key-space result as the
# output rather than failing at the last step.
if notExists "$OUT"; then
    if [ -f "$DB.lookup" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" translatekeys "$TMP/clu.keys.tsv" "$DB.lookup" "$OUT.tmp" \
            --threads $THREADS --split-memory-limit $SPLIT_MEMORY_LIMIT \
            || fail "translatekeys died"
        mv -f "$OUT.tmp" "$OUT"
    else
        echo "No $DB.lookup; leaving the result in database keys"
        cp "$TMP/clu.keys.tsv" "$OUT"
    fi
fi

echo "Wrote $OUT"
