#!/bin/sh -e
# Shared-filesystem parallel linclust.
#
# Runs the two linclust passes across many independent worker processes that
# coordinate only through files: no MPI, no rank argument, no node-to-node
# communication. Every worker of a stage runs a byte-identical command line, so a
# stage maps onto a Slurm array job and workers may join late, die, or be
# restarted.
#
#   WORKER_RUNNER="srun -n 64" ./linclustparallel.sh input.fasta clusters.tsv tmp
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
# WORKER_RUNNER, not RUNNER: RUNNER already means the MPI runner to ten other
# workflows and to Parameters.cpp, which seeds par.runner from it.
[ -z "$WORKER_RUNNER" ] && WORKER_RUNNER=""
# getconf, not nproc: nproc is GNU-only and absent on macOS/BSD.
[ -z "$THREADS" ] && THREADS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)
[ -z "$MIN_SEQ_ID" ] && MIN_SEQ_ID=0.9
[ -z "$COV" ] && COV=0.8
[ -z "$COV_MODE" ] && COV_MODE=1
[ -z "$EVAL" ] && EVAL=0.001
[ -z "$SPLIT_MEMORY_LIMIT" ] && SPLIT_MEMORY_LIMIT=0
[ -z "$SCRATCH_BUDGET" ] && SCRATCH_BUDGET=0
# How many workers $WORKER_RUNNER starts. Passed to every stage that decides how
# many work units to make, so the allocation has something to claim: those counts
# otherwise come from --split-memory-limit alone, and a generous per-node limit
# makes *fewer* units and leaves most workers idle.
#
# Derived from the allocation when the caller did not say, because "not told" is
# the normal case and it is silent: the documented invocation on line 10 sets
# WORKER_RUNNER and nothing else, so every stage sized its work units for one
# worker while 64 ran, and every one of them reported success. Slurm exports both
# of these; outside Slurm the count stays 0 and the stages say so themselves.
if [ -z "$WORKERS" ]; then
    WORKERS="${SLURM_NTASKS:-${SLURM_JOB_NUM_NODES:-0}}"
fi
# Writes the k-mer and candidate-edge buckets as fixed-width structs instead of
# the packed block encoding. Roughly doubles the scratch those two intermediates
# need and exists only as a control: a run with RAW_RECORDS=1 must produce
# byte-identical output to one without, which is what separates an encoding
# defect from a semantic one.
[ -z "$RAW_RECORDS" ] && RAW_RECORDS=0
# Forces the reduce to group each partition in this many k-mer slices. 0 derives
# the count from the memory budget. Slicing is exact -- a slice is a pure function
# of the k-mer, so a group is never split -- so any value must give byte-identical
# output, which is how that exactness is tested.
[ -z "$REDUCE_SLICES" ] && REDUCE_SLICES=0
# Which stage groups this invocation runs, space separated, in the order they
# appear below:
#
#   p1  createdbparallel, the pass-1 k-mer waves, alignparallel      multi-worker
#   s1  greedycluster (pass 1), createrepdb                          one node
#   p2  the pass-2 k-mer waves, alignparallel                        multi-worker
#   s2  greedycluster (pass 2)                                       one node
#   p3  translatecluster, mergeclusterparallel                       multi-worker
#   s3  translatekeys                                                one node
#
# The default is all of them, which is the single-allocation behaviour and the
# only one a run outside a batch system wants.
#
# Splitting them exists for the opposite case. The p-phases divide across the
# allocation; the s-phases are one process on one node however many were asked
# for, so in a single N-node job the other N-1 sit idle through them and are
# billed for it. At 1e10 the s-phases are ~1.9 h of a 4.7 h run on 32 nodes:
# 39% of the allocation buys nothing. Submitting the phases as a dependency
# chain, the p-phases on N nodes and the s-phases on one, makes the bill flat in
# N instead of rising with it.
#
# Every phase is idempotent, so the split is safe rather than merely cheaper.
# The multi-worker stages resume from their work queues and exit at once when
# the queue is already drained; the single-node stages are guarded on their
# output. Running a phase twice, or running the whole pipeline in one job after
# a partial chain, therefore costs a queue scan and nothing else.
[ -z "$LINCLUST_PHASES" ] && LINCLUST_PHASES="p1 s1 p2 s2 p3 s3"
for _phase in $LINCLUST_PHASES; do
    case "$_phase" in
        p1|s1|p2|s2|p3|s3) ;;
        *) echo "Error: unknown phase '$_phase' in LINCLUST_PHASES. Expected a subset of \
'p1 s1 p2 s2 p3 s3'." >&2; exit 1 ;;
    esac
done
# Substring match on a padded copy, so "p1" cannot match inside another token.
phase() { case " $LINCLUST_PHASES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

notExists() { [ ! -f "$1" ]; }
# To stderr, not stdout: scratchUsed() runs inside a command substitution, so a
# failure message written to stdout would be captured as the measurement instead
# of being shown, and the run would die with no explanation at all.
fail() { echo "Error: $1" >&2; exit 1; }

# Bytes currently sitting on the scratch filesystem for this run.
#
# The wave count has to be derived against what is actually free, not against the
# whole budget: pass 2 runs with all of pass 1's surviving output still on disk,
# and sizing it as though the filesystem were empty is what made a 1e11 run derive
# a single wave and then peak at ~1.9x its ceiling. Measuring beats modelling here
# -- it needs no per-stage accounting and stays right when the pipeline changes.
#
# du -sk rather than -sb: -b is GNU-only. Kibibytes are ample precision against a
# budget in terabytes. printf "%.0f" keeps awk from rendering large sums in
# scientific notation, which --scratch-used would reject.
#
# Fails closed. This used to be one `du ... | awk ...` pipeline without pipefail:
# when du failed -- an unreadable path, a vanished directory, a stale mount -- awk
# still succeeded and printed "0B", so the caller derived its wave count as though
# the filesystem were empty and planned a run that could not fit. Reproduced with a
# missing path: exit status 0, used=0B. So du runs on its own, its status is
# checked, and a partial measurement is refused rather than rounded down to zero.
#
# $DB is a prefix, not a path, so the index, .index.bin, .lookup and headers have
# to be counted too -- at 1e11 the index alone is over a terabyte of what this
# measures. When $DB lives under $TMP (the usual case: the workflow builds it at
# $TMP/db) only $TMP is passed: GNU du de-duplicates repeated paths, but BSD/macOS
# du counts the database twice, and the surrounding code is otherwise careful to
# stay portable.
scratchUsed() {
    # The path list is built in the positional parameters rather than in a string,
    # so a directory with a space in it is one argument and not two.
    set -- "$TMP"
    case "$DB" in
        "$TMP"/*) ;;  # already inside $TMP; passing it again double-counts on BSD du
        *)
            for _p in "$DB"*; do
                if [ -e "$_p" ]; then set -- "$@" "$_p"; fi
            done
            ;;
    esac
    # `|| _du=""` rather than testing $? afterwards: under `sh -e` a command
    # substitution that exits non-zero aborts the script at the assignment, so the
    # check would never run. Discarding a partial measurement is the point --
    # anything less than a complete figure has to be refused, not rounded down.
    _du=$(du -sk "$@" 2>/dev/null) || _du=""
    if [ -z "$_du" ]; then
        fail "cannot measure the scratch already in use under $TMP.
--scratch-budget is a hard ceiling on this run, and a measurement that silently
returned zero would let the k-mer shuffle be planned as though the filesystem were
empty. Check that $TMP and $DB* are readable."
    fi
    printf '%s\n' "$_du" | awk '{s += $1} END {printf "%.0fB\n", s * 1024}'
}

# Deletes an intermediate whose consumers have all finished. Nothing downstream
# reopens these, and they are the bulk of peak scratch: at 100M the candidate
# edges alone are 21 GB and the pass-1 alignments 3 GB.
#
# KEEP_INTERMEDIATE=1 suppresses it. Env-only and deliberately not a parameter:
# it exists to inspect a run's intermediates while debugging, and keeping them
# breaks the scratch budget the run was sized against.
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
# `|| true` so a missing or unreadable manifest reaches the check below rather than
# aborting the script under `sh -e` with awk's own error, which is what used to
# happen and made the friendly message unreachable.
waveCount() { awk '$1 == "waveCount" { print $2 }' "$1/coord/shuffle.info" 2>/dev/null || true; }

mapReduceWaves() {
    # $1 sequence DB, $2 k-mer dir, $3 edge dir, $4... extra map arguments
    _db="$1"; _kmer="$2"; _edges="$3"; shift 3
    # shellcheck disable=SC2086
    _used=$(scratchUsed)
    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" kmermatcherparallel "$_db" "$_kmer" $KMER_COMMON "$@" --kmer-wave 0 \
        --scratch-used "$_used" \
        || fail "kmermatcherparallel died"
    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" kmerreduceparallel "$_db" "$_kmer" "$_edges" $REDUCE_PAR --kmer-wave 0 \
        || fail "kmerreduceparallel died"
    _waves=$(waveCount "$_kmer")
    [ -z "$_waves" ] && fail "no wave count in $_kmer/coord/shuffle.info"
    _w=1
    while [ "$_w" -lt "$_waves" ]; do
        echo "--- extraction wave $_w of $_waves ---"
        # shellcheck disable=SC2086
        $WORKER_RUNNER "$MMSEQS" kmermatcherparallel "$_db" "$_kmer" $KMER_COMMON "$@" --kmer-wave $_w \
            --scratch-used "$_used" \
            || fail "kmermatcherparallel (wave $_w) died"
        # shellcheck disable=SC2086
        $WORKER_RUNNER "$MMSEQS" kmerreduceparallel "$_db" "$_kmer" "$_edges" $REDUCE_PAR --kmer-wave $_w \
            || fail "kmerreduceparallel (wave $_w) died"
        _w=$((_w + 1))
    done
}

[ "$#" -ne 3 ] && { echo "usage: [WORKER_RUNNER=\"srun -n N\"] $0 <i:fasta|sequenceDB> <o:clusterTsv> <tmpDir>"; exit 1; }
INPUT="$1"
OUT="$2"
TMP="$3"
# Refused rather than skipped, as linclust.sh:13 does. Guarding the final write on
# `notExists "$OUT"` instead meant a stale file survived a full run and the script
# still reported success, so the user believed they had clustered and had not.
# Resume is keyed on state under $TMP, never on $OUT.
# -e, not -f: a directory at $OUT passes a -f test, and the `mv -f` at the end
# then files the result *inside* it and reports success.
[ -e "$OUT" ] && fail "$OUT exists already"
mkdir -p "$TMP"

# Shared by both passes. --min-seq-id belongs here because it selects the k-mer
# length and alphabet; coverage does not, because this stage only extracts k-mers.
KMER_COMMON="--alph-size aa:21,nucl:5 --min-seq-id $MIN_SEQ_ID --kmer-per-seq 21 \
 --mask 0 --mask-prob 0.9 --mask-lower-case 0 --mask-n-repeat 0 -k 0 --max-seq-len 65535 \
 --hash-shift 67 --ignore-multi-kmer 0 \
 --split-memory-limit $SPLIT_MEMORY_LIMIT --scratch-budget $SCRATCH_BUDGET \
 --raw-records $RAW_RECORDS --workers $WORKERS --threads $THREADS"
REDUCE_PAR="-c $COV --cov-mode $COV_MODE --include-adjacency 1 --num-adjacency 3 \
 --split-memory-limit $SPLIT_MEMORY_LIMIT --raw-records $RAW_RECORDS \
 --reduce-slices $REDUCE_SLICES --workers $WORKERS --threads $THREADS"
# --alignment-mode 3 is SCORE_COV_SEQID, which is what stock linclust forces
# (Parameters.cpp, ALIGN2CLUST_PAR). It has to be passed rather than assumed: the
# workflow process sets it as a default, but alignparallel is a separate process
# and used its own FAST_AUTO -- which at --min-seq-id 0 degrades to SCORE_COV and
# at -c 0 to SCORE_ONLY, reporting the score-per-column estimate as the identity.
# Both are silent, and both make this pipeline cluster differently from stock.
ALIGN_PAR="--min-seq-id $MIN_SEQ_ID --min-aln-len 0 --seq-id-mode 0 -e $EVAL -c $COV \
 --cov-mode $COV_MODE --alignment-mode 3 --threads $THREADS"

# 0. Sequence database with dense, length-ranked keys. Every later stage depends
#    on that ordering: it is what makes the greedy a single forward sweep.
DB="$INPUT"
# .lookup rather than .index.bin as the second test: createdbparallel writes the
# dense index in its plan phase but .dbtype only in finalize, so a user database
# whose build was killed by a walltime has .index.bin and no .dbtype, and one
# killed slightly earlier has neither. Either way the guard below has to fire
# rather than hand a half-built database to the next stage as though it were
# FASTA. The workflow's own build is covered by its finalize.done sentinel.
if [ -f "$INPUT.dbtype" ] && notExists "$INPUT.index.bin"; then
    # Checked here rather than left to the first stage. Reaching for `createdb`
    # first is the natural thing to do, and its databases are not dense or
    # length-ranked, which every stage below assumes. The stage does say so, but
    # under a runner that message lands in one worker's log after the run has
    # already started; this fails before anything is launched.
    fail "$INPUT was not built by createdbparallel: no $INPUT.index.bin.
Every stage addresses sequences as dense, length-ranked keys, which is what makes
the greedy a single forward sweep, and a createdb database has neither. Either
pass the FASTA and let this build the database, or build it with
'mmseqs createdbparallel'."
fi
if notExists "$INPUT.dbtype"; then
    DB="$TMP/db"
    # Guarded on the finalize sentinel, not on .dbtype: createdbparallel writes
    # .dbtype before the text indices, so a death in that window leaves a database
    # that looks complete and is not. The sentinel is written last.
    #
    # Every phase needs $DB resolved, but only p1 may build it: a later phase
    # reaching this with no database is a chain submitted out of order, and
    # rebuilding it there would run a multi-worker stage on the one node the
    # s-phases are given.
    if notExists "$DB.coord/finalize.done" && ! phase p1; then
        fail "$DB has no finished database and this invocation does not run phase \
p1, which builds it (LINCLUST_PHASES=\"$LINCLUST_PHASES\").
Run the phases in order: p1 s1 p2 s2 p3 s3."
    fi
    if notExists "$DB.coord/finalize.done"; then
        # shellcheck disable=SC2086
        # --write-text-index 0: no stage of this pipeline reads a text .index --
        # they all address entries through the dense .index.bin -- and generating
        # it is one snprintf and one buffered write per entry, single-threaded, for
        # both the sequence and the header database. Measured at 177 ns per line,
        # that is ~98 h and ~36 TB at 1e12 for output nothing here consumes. Build
        # the database with `mmseqs createdbparallel --write-text-index 1` if a
        # stock MMseqs2 tool has to open it afterwards.
        $WORKER_RUNNER "$MMSEQS" createdbparallel "$INPUT" "$DB" $VERBOSITY --threads $THREADS \
            --write-text-index 0 --write-header-db 0 \
            || fail "createdbparallel died"
    fi
fi

# Checked here, not by alignparallel, which is where it used to surface. That is
# the far end of pass 1, so a nucleotide input paid createdbparallel plus the whole
# map and reduce -- hours and the entire k-mer shuffle on scratch -- before being
# told it was never supported. Stock routes nucleotides to linclust v1, which has
# no v2 counterpart to port.
case "$("$MMSEQS" dbtype "$DB" 2>/dev/null || echo unknown)" in
    *Nucleotide*)
        fail "$DB holds nucleotide sequences, which this pipeline does not support.
linclust clusters nucleotides with its v1 path (rescorediagonal), and this is the
parallel form of the v2 path, which is protein-only upstream."
        ;;
esac

# ---- pass 1, over the whole database --------------------------------------
# Guarded as a whole on clu1.tsv. The stages inside are individually resumable, but
# the pass is not re-enterable once its k-mer buckets have been consumed and its
# intermediates dropped: the reduce would then run over an empty shuffle, emit zero
# edges and exit 0, and a greedycluster over that produces a clustering of
# singletons rather than an error.
if phase p1 && notExists "$TMP/clu1.tsv"; then
    mapReduceWaves "$DB" "$TMP/kmer1" "$TMP/edges1" \
        --spaced-kmer-mode 0 --kmer-per-seq-scale aa:0.000,nucl:0.200

    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" alignparallel "$DB" "$TMP/edges1" "$TMP/aln1" $ALIGN_PAR $VERBOSITY \
        || fail "alignparallel (pass 1) died"
fi

# Single node: the greedy sweep is sequential by necessity, which is what makes
# it exact. Its own phase so an allocation sized for the k-mer waves does not
# have to be held through it.
if phase s1 && notExists "$TMP/clu1.tsv"; then
    if notExists "$TMP/aln1/coord/align.info"; then
        fail "no pass-1 alignments in $TMP/aln1. Phase p1 has not finished."
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" greedycluster "$DB" "$TMP/aln1" "$TMP/clu1.tsv.tmp" $VERBOSITY --threads $THREADS \
        || fail "greedycluster (pass 1) died"
    mv -f "$TMP/clu1.tsv.tmp" "$TMP/clu1.tsv"
fi
# Dead once clu1.tsv exists, and together the largest thing pass 1 leaves behind
# for pass 2 to be sized around. kmer1 is already gone: the reduce unlinks each
# wave's buckets as it consumes them.
#
# Conditioned on clu1.tsv rather than run unconditionally, which is what it used
# to be. With the phases split, p1 ends with the alignments still being the only
# copy of pass 1's work and greedycluster not yet run: dropping them there would
# delete the input of the very next phase, and the pass is guarded as
# unre-enterable precisely so that it cannot be rebuilt.
if [ -f "$TMP/clu1.tsv" ]; then
    dropIntermediate "$TMP/edges1" "$TMP/aln1"
fi

# ---- pass 2, over the representatives --------------------------------------
# Re-keyed densely rather than sub-set, because every stage addresses keys as
# array offsets. The key map translates the result back below.
# Single node, and in the same phase as the greedy sweep it consumes: both are
# one process reading the whole key space, so a chain that separated them would
# pay a queue wait for nothing.
if phase s1 && notExists "$TMP/rep.keymap"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createrepdb "$DB" "$TMP/clu1.tsv" "$TMP/rep" $VERBOSITY --threads $THREADS \
        --split-memory-limit $SPLIT_MEMORY_LIMIT \
        || fail "createrepdb died"
fi

# Guarded as a whole, for the same reason as pass 1.
if phase p2 && notExists "$TMP/clu2_sub.tsv"; then
    if notExists "$TMP/rep.keymap"; then
        fail "no representative database at $TMP/rep. Phase s1 has not finished."
    fi
    mapReduceWaves "$TMP/rep" "$TMP/kmer2" "$TMP/edges2" \
        --spaced-kmer-mode 1 --kmer-per-seq-scale aa:0.100,nucl:0.100

    # The filter gate: a representative may only take a member if every sequence of
    # that member's pass-1 cluster also aligns to it. Needs the original database
    # and the key map, since the clustering it consults is in original keys.
    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" alignparallel "$TMP/rep" "$TMP/edges2" "$TMP/aln2" $ALIGN_PAR \
        $VERBOSITY \
        --filter-cludb-file "$TMP/clu1.tsv" --filter-seqdb-file "$DB" \
        --key-map "$TMP/rep.keymap" \
        || fail "alignparallel (pass 2) died"
fi

# Single node, like its pass-1 counterpart.
if phase s2 && notExists "$TMP/clu2_sub.tsv"; then
    if notExists "$TMP/aln2/coord/align.info"; then
        fail "no pass-2 alignments in $TMP/aln2. Phase p2 has not finished."
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" greedycluster "$TMP/rep" "$TMP/aln2" "$TMP/clu2_sub.tsv.tmp" $VERBOSITY \
        --threads $THREADS \
        || fail "greedycluster (pass 2) died"
    mv -f "$TMP/clu2_sub.tsv.tmp" "$TMP/clu2_sub.tsv"
fi
# Conditioned for the same reason as pass 1's: p2 leaves aln2 as the only copy
# of its work, and s2 is what consumes it.
if [ -f "$TMP/clu2_sub.tsv" ]; then
    dropIntermediate "$TMP/edges2" "$TMP/aln2" "$TMP/kmer2"
fi

if phase p3 && notExists "$TMP/clu2.tsv"; then
    if notExists "$TMP/clu2_sub.tsv"; then
        fail "no pass-2 clustering at $TMP/clu2_sub.tsv. Phase s2 has not finished."
    fi
    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" translatecluster "$TMP/clu2_sub.tsv" "$TMP/rep.keymap" \
        "$TMP/clu2.tsv.tmp" \
        $VERBOSITY --split-memory-limit $SPLIT_MEMORY_LIMIT --workers $WORKERS \
        || fail "translatecluster died"
    # After the runner has reaped every worker: the last phase is a parallel
    # pwrite into this path by all of them.
    mv -f "$TMP/clu2.tsv.tmp" "$TMP/clu2.tsv"
    rm -rf "$TMP/clu2.tsv.tmp.coord"
fi

# ---- fold pass 2 into pass 1 ------------------------------------------------
# Keys, not accessions: every stage addresses sequences as dense keys, so this is
# the form the merge produces. It is translated below.
if phase p3 && notExists "$TMP/clu.keys.tsv"; then
    # shellcheck disable=SC2086
    $WORKER_RUNNER "$MMSEQS" mergeclusterparallel "$DB" "$TMP/clu1.tsv" "$TMP/clu2.tsv" \
        "$TMP/clu.keys.tsv.tmp" \
        $VERBOSITY --split-memory-limit $SPLIT_MEMORY_LIMIT --workers $WORKERS \
        || fail "mergeclusterparallel died"
    # After the runner has reaped every worker, so the file is complete: the last
    # phase of the merge is a parallel pwrite into this path by all of them.
    mv -f "$TMP/clu.keys.tsv.tmp" "$TMP/clu.keys.tsv"
    rm -rf "$TMP/clu.keys.tsv.tmp.coord"
fi

# ---- back to accession space ------------------------------------------------
# Stock reaches this with createtsv, which needs a result database and a resident
# id->name table. Here it is a streaming join against the .lookup. Skipped when
# the database has none (--write-lookup 0), leaving the key-space result as the
# output rather than failing at the last step.
if ! phase s3; then
    echo "Phases '$LINCLUST_PHASES' done; $OUT is written by phase s3"
    exit 0
fi
if notExists "$TMP/clu.keys.tsv"; then
    fail "no merged clustering at $TMP/clu.keys.tsv. Phase p3 has not finished."
fi
if [ -f "$DB.lookup" ]; then
    # shellcheck disable=SC2086
    # --spill-prefix under $TMP. translatekeys derives its spill files from the
    # output path when it is not told otherwise, which puts every intermediate it
    # writes -- at 1e11 the first fixed-width spill alone is ~1.6 TB -- on the
    # *output* filesystem, outside everything --scratch-budget accounts for. Only
    # the finished file is published next to $OUT, by the rename below.
    "$MMSEQS" translatekeys "$TMP/clu.keys.tsv" "$DB.lookup" "$OUT.tmp" $VERBOSITY \
        --threads $THREADS --split-memory-limit $SPLIT_MEMORY_LIMIT \
        --spill-prefix "$TMP/translatekeys" \
        || fail "translatekeys died"
    mv -f "$OUT.tmp" "$OUT"
else
    echo "No $DB.lookup; leaving the result in database keys"
    # Copied to a temporary name and renamed, like every other output here: an
    # interrupted cp straight onto $OUT leaves a truncated file behind.
    cp "$TMP/clu.keys.tsv" "$OUT.tmp" || fail "cannot copy the key-space result"
    mv -f "$OUT.tmp" "$OUT"
fi

# Only on success, and only the run's own intermediates. The hashed directory
# itself and its `latest` symlink are left in place, as every other workflow does,
# so --force-reuse still resolves; `rm -rf "$TMP"` used to leave `latest` dangling.
# The per-stage intermediates are already dropped as they die (dropIntermediate
# above), which is what keeps the run inside --scratch-budget.
if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -rf "$TMP/kmer1" "$TMP/kmer2" "$TMP/edges1" "$TMP/edges2" "$TMP/aln1" "$TMP/aln2"
    rm -f "$TMP/clu1.tsv" "$TMP/clu2.tsv" "$TMP/clu2_sub.tsv" "$TMP/clu.keys.tsv"
    rm -f "$TMP/rep" "$TMP/rep".* "$TMP/rep_h" "$TMP/rep_h".*
    rm -rf "$TMP/db.coord"
    rm -f "$TMP/db" "$TMP/db".* "$TMP/db_h" "$TMP/db_h".*
fi

echo "Wrote $OUT"
