#!/bin/sh -e
# One shot linear clustering. Run it once on every machine with the same NODES and its own NODE.
# Guards are on the per machine markers, not on the output every machine publishes into.
fail() { echo "$1" >&2; exit 1; }
notExists() { [ ! -f "$1" ]; }

# how many machines published their part, so a wait can speak when that changes and not every second
finishedCount() {
    _n=0; _i=0
    while [ "$_i" -lt "$2" ]; do
        [ -f "$1.$_i.done" ] && _n=$((_n + 1))
        _i=$((_i + 1))
    done
    echo "$_n"
}

# waits for every machine's marker, because a machine that has not started looks like one with nothing
waitForAll() {
    _path="$1"; _n="$2"; _waited=0; _seen=""
    while true; do
        _have=$(finishedCount "$_path" "$_n")
        [ "$_have" -eq "$_n" ] && return 0
        [ "$_have" != "$_seen" ] && echo "$_have of $_n machines finished $_path"
        _seen="$_have"
        [ "$_waited" -ge "$WAIT_LIMIT" ] && fail "waited ${WAIT_LIMIT}s for $_path, $_have of $_n machines finished"
        sleep 1
        _waited=$((_waited + 1))
    done
}

# a counter one machine rewrites in place, because on NFS a name that never existed stays cached as
# missing while a name that does gets revalidated on open
readProgress() {
    _pv=0
    read -r _pv 2>/dev/null < "$1" || _pv=0
    _pv=$(printf '%s' "$_pv" | tr -cd '0-9')
    _pv=${_pv#"${_pv%%[!0]*}"}
    printf '%s' "${_pv:-0}"
}

# every machine has to reach the count, because one that has not started looks like one with nothing

# the first round speaks on the console; the rest keep their say in a per machine file
quietOnRetry() {
    if [ -z "$_said" ]; then
        "$@"; _rc=$?
        echo "the later rounds report to $TMP/lin8clust.$NODE.log"
        return $_rc
    fi
    "$@" >> "$TMP/lin8clust.$NODE.log"
}

# arguments when called as a command, IN/OUT/TMP in the environment when called on its own
if [ "$#" -ge 3 ]; then
    OUT="$(eval echo "\${$(($# - 1))}")"
    TMP="$(eval echo "\${$#}")"
    IN=""
    _i=1
    while [ "$_i" -le "$(($# - 2))" ]; do
        IN="$IN $(eval echo "\${$_i}")"
        _i=$((_i + 1))
    done
fi
[ -z "$IN" ] && fail "no input sequence file given"
[ -z "$OUT" ] && fail "no output cluster database given"
[ -z "$TMP" ] && fail "no tmp directory given"
[ ! -d "$TMP" ] && mkdir -p "$TMP"
NODES="${NODES:-1}"
NODE="${NODE:-0}"
REP_RANK_BLOCKS="${REP_RANK_BLOCKS:-1024}"
THREADS="${THREADS:-1}"
SEQID="${SEQID:-0.9}"
# the redundancy pass may fold at any threshold at or above the target, and 0.9 is the loosest
# the original allows; the workflow sets this, so this is for a script run on its own
HASHSEQID="${HASHSEQID:-$(awk -v s="$SEQID" 'BEGIN{print (s>0.9)?s:0.9}')}"
COV="${COV:-0.8}"
COVMODE="${COVMODE:-0}"
REPSEQ="${REPSEQ:-0}"
WAIT_LIMIT="${WAIT_LIMIT:-86400}"
# how many blocks one wave step takes, so the assignment bitmap is carried in and out once for all
# the switched representative comes out of the alignments, so they have to have been kept
[ -n "$SWITCH_CONSENSUS_REP" ] && [ -z "$INCLUDE_ALIGN_FILES" ] && \
    fail "--switch-consensus-rep needs --include-align-files 1, which keeps the alignments it scores"
ALIGN_TEXT=""
[ -n "$INCLUDE_ALIGN_FILES" ] && ALIGN_TEXT="--include-align-files 1"

mkdir -p "$TMP/kmer" "$TMP/pairs" "$TMP/pref" "$TMP/aln" "$TMP/clu_accepted"

# every machine samples its own process tree from /proc while the passes run, one line an interval,
# appended beside the output so a rerun continues the same file; a summary a phase closes it
MONITOR_INTERVAL="${MONITOR_INTERVAL:-5}"
MONITOR_PID=""
stopMonitor() {
    if [ -n "$MONITOR_PID" ]; then
        kill "$MONITOR_PID" 2>/dev/null || true
        wait "$MONITOR_PID" 2>/dev/null || true
        MONITOR_PID=""
    fi
    return 0
}
if [ "$MONITOR_INTERVAL" -gt 0 ] 2>/dev/null; then
    "$MMSEQS" lin8-monitor "$OUT.monitor.$NODE.tsv" --monitor-pid "$$" \
        --monitor-interval "$MONITOR_INTERVAL" -v 1 &
    MONITOR_PID=$!
    trap stopMonitor EXIT
fi

# two rounds: every machine writes a histogram, then places sequences, so this is a retry loop
_waited=0
_said=""
while notExists "$TMP/db.$NODE.done" || notExists "$TMP/db.dbtype"; do
    # shellcheck disable=SC2086
    "$MMSEQS" lin8-createdb $IN "$TMP/db" --node-count "$NODES" --node-id "$NODE" \
        --threads "$THREADS" ${CREATEDB_PAR}
    if notExists "$TMP/db.$NODE.done" || notExists "$TMP/db.dbtype"; then
        _have=$(finishedCount "$TMP/db" "$NODES")
        [ "$_have" != "$_said" ] && echo "$_have of $NODES machines placed their sequences"
        _said="$_have"
        [ "$_waited" -ge "$WAIT_LIMIT" ] && fail "waited ${WAIT_LIMIT}s for the database"
        sleep 1
        _waited=$((_waited + 1))
    fi
done

# the same two rounds, ending in the bitmap
if [ -n "$CLUSTHASH" ]; then
    # each machine folds once; machine 0 then waits inside the binary and merges. A rerun that
    # finds its shard already there re-enters cheaply and node 0 still produces the bitmap.
    if notExists "$TMP/db.clusthash_kept"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-clusthash "$TMP/db" "$TMP/hash" --node-count "$NODES" --node-id "$NODE" \
            --threads "$THREADS" --min-seq-id "$HASHSEQID" ${HASH_PAR} || fail "lin8-clusthash died"
    fi
    _waited=0
    while notExists "$TMP/db.clusthash_kept"; do
        [ "$_waited" -ge "$WAIT_LIMIT" ] && fail "waited ${WAIT_LIMIT}s for the redundancy pass"
        sleep 1
        _waited=$((_waited + 1))
    done

    [ "$NODES" -gt 1 ] && waitForAll "$TMP/hash" "$NODES"
fi

if notExists "$TMP/kmer/out.$NODE.done"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lin8-extractkmers "$TMP/db" "$TMP/kmer/out" --node-count "$NODES" --node-id "$NODE" \
        --threads "$THREADS" --min-seq-id "$SEQID" ${EXTRACT_PAR} || fail "lin8-extractkmers died"
fi

[ "$NODES" -gt 1 ] && waitForAll "$TMP/kmer/out" "$NODES"

if notExists "$TMP/pairs/pairs.$NODE.done"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lin8-assignedpairs "$TMP/db" "$TMP/kmer/out" "$TMP/pairs/pairs" \
        --node-count "$NODES" --node-id "$NODE" --threads "$THREADS" --pair-splits "$REP_RANK_BLOCKS" \
        -c "$COV" --cov-mode "$COVMODE" ${GROUP_PAR} || fail "lin8-assignedpairs died"
fi

[ "$NODES" -gt 1 ] && waitForAll "$TMP/pairs/pairs" "$NODES"

if notExists "$TMP/pref/pref.$NODE.done"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lin8-pref "$TMP/pairs/pairs" "$TMP/db" "$TMP/pref/pref" \
        --node-count "$NODES" --node-id "$NODE" --threads "$THREADS" ${PREF_PAR} || fail "lin8-pref died"
fi

[ "$NODES" -gt 1 ] && waitForAll "$TMP/pref/pref" "$NODES"

# The wave: block k only needs what was decided up to k-1, and deciding is one thread because a
# greedy assignment is one order. One machine fuses the two halves; many machines barrier here.
if [ "$NODES" -eq 1 ]; then
    if notExists "$TMP/clu_accepted/clu_accepted.0.$((REP_RANK_BLOCKS - 1))"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-align2clust "$TMP/db" "$TMP/pref/pref" "$TMP/aln/aln" \
            "$TMP/clu_accepted/clu_accepted" --node-count 1 --node-id 0 \
            --threads "$THREADS" --min-seq-id "$SEQID" -c "$COV" --cov-mode "$COVMODE" \
            ${ALIGN_TEXT} ${ALIGN_PAR} || fail "lin8-align2clust died"
    fi
    R="$REP_RANK_BLOCKS"
else
    R=0
fi
# aligning may start this far before deciding caught up, so the two overlap; a block the bitmap has
# not seen yet only costs the pairs deciding then throws away
# 0 or unset takes the whole wave in one invocation a machine; the decider trails block by block
PAIR_SPLIT_COUNT="${PAIR_SPLIT_COUNT:-0}"
{ [ "$PAIR_SPLIT_COUNT" -le 0 ] || [ "$PAIR_SPLIT_COUNT" -gt "$REP_RANK_BLOCKS" ]; } && PAIR_SPLIT_COUNT=$REP_RANK_BLOCKS
LOOKAHEAD="${PAIR_SPLIT_LOOKAHEAD:-8}"
[ "$NODES" -eq 1 ] && LOOKAHEAD=0
if [ "$NODES" -gt 1 ]; then
    [ "$NODE" -eq 0 ] && notExists "$TMP/clu_accepted/clu_accepted.progress" && \
        printf '%020d\n' 0 > "$TMP/clu_accepted/clu_accepted.progress"
fi
MULTI_PID=""
# one block speaks for the wave; the rest repeat it and stay quiet. Errors are stderr and always show
_said=""
# a coarse heartbeat: a chunked wave reports at ~16 milestones, not once per chunk
_beat=$((REP_RANK_BLOCKS / 16)); [ "$_beat" -lt 1 ] && _beat=1
_beaten=-1
while [ "$R" -lt "$REP_RANK_BLOCKS" ]; do
    _last=$((R + PAIR_SPLIT_COUNT - 1))
    [ "$_last" -ge "$REP_RANK_BLOCKS" ] && _last=$((REP_RANK_BLOCKS - 1))
    if notExists "$TMP/clu_accepted/clu_accepted.0.$_last"; then
        if [ "$NODE" -eq 0 ]; then
            # the decider starts first and trails the aligners on the markers it waits for
            # shellcheck disable=SC2086
            quietOnRetry "$MMSEQS" lin8-align2clustmulti "$TMP/aln/aln" "$TMP/pref/pref" \
                "$TMP/clu_accepted/clu_accepted" --pair-split "$R" \
                --pair-split-count "$PAIR_SPLIT_COUNT" ${ALIGN_TEXT} ${ASSIGN_PAR} &
            MULTI_PID=$!
        fi
        if notExists "$TMP/aln/aln.$_last.$NODE.done"; then
            # shellcheck disable=SC2086
            quietOnRetry "$MMSEQS" lin8-align2clust "$TMP/db" "$TMP/pref/pref" "$TMP/aln/aln" \
                "$TMP/clu_accepted/clu_accepted" \
                --node-count "$NODES" --node-id "$NODE" --pair-split "$R" \
                --pair-split-count "$PAIR_SPLIT_COUNT" --pair-split-lookahead "$LOOKAHEAD" \
                --threads "$THREADS" --min-seq-id "$SEQID" -c "$COV" --cov-mode "$COVMODE" \
                ${ALIGN_TEXT} ${ALIGN_PAR} || fail "lin8-align2clust died"
        fi
        if [ "$NODE" -eq 0 ] && [ -n "$MULTI_PID" ]; then
            wait "$MULTI_PID" || fail "lin8-align2clustmulti died"
            MULTI_PID=""
        fi
        if [ "$_last" -eq "$((REP_RANK_BLOCKS - 1))" ] || [ "$((_last / _beat))" -ne "$_beaten" ]; then
            echo "blocks up to $((_last + 1)) of $REP_RANK_BLOCKS aligned"
            _beaten=$((_last / _beat))
        fi
        if [ "$NODE" -ne 0 ]; then
            # the next block needs the bitmap up to the lookahead, not up to this block
            _need=$((_last + 1 - LOOKAHEAD))
            [ "$_need" -lt 0 ] && _need=0
            _waited=0
            _seenp=""
            while true; do
                _at=$(readProgress "$TMP/clu_accepted/clu_accepted.progress")
                [ "$_at" -ge "$_need" ] && break
                [ "$_at" != "$_seenp" ] && echo "$_at of $REP_RANK_BLOCKS blocks decided, this machine needs $_need"
                _seenp="$_at"
                [ "$_waited" -ge "$WAIT_LIMIT" ] && fail "waited ${WAIT_LIMIT}s for repRankBlock $((_need - 1)) to be decided"
                sleep 1
                _waited=$((_waited + 1))
            done
        fi
    fi
    _said=1
    R=$((_last + 1))
done

# the lookahead let every machine run past the last decided block, and they all hold what it needs
if [ "$NODE" -ne 0 ] && [ "$NODES" -gt 1 ]; then
    _waited=0
    while notExists "$TMP/clu_accepted/clu_accepted.0.$((REP_RANK_BLOCKS - 1))"; do
        [ "$_waited" -ge "$WAIT_LIMIT" ] && fail "waited ${WAIT_LIMIT}s for every repRankBlock to be decided"
        sleep 1
        _waited=$((_waited + 1))
    done
fi

# the redundant sequences go back in on the pair stream, and the clustering database is made from that
if [ "$NODE" -eq 0 ]; then
    CLUDB="$TMP/clu_accepted/clu_accepted"
    if [ -n "$CLUSTHASH" ]; then
        mkdir -p "$TMP/clu_accepted_plus_redundant"
        if notExists "$TMP/clu_accepted_plus_redundant/clu_accepted_plus_redundant"; then
            # shellcheck disable=SC2086
            "$MMSEQS" lin8-mergehashredundancy "$CLUDB" "$TMP/hash" \
                "$TMP/clu_accepted_plus_redundant/clu_accepted_plus_redundant" ${EXPAND_PAR} \
                || fail "lin8-mergehashredundancy died"
        fi
        CLUDB="$TMP/clu_accepted_plus_redundant/clu_accepted_plus_redundant"
    fi
    if notExists "$OUT.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-createclusterdb "$CLUDB" "$OUT" ${CLUSTERDB_PAR} || fail "lin8-createclusterdb died"
    fi

    # the alignments name the member that speaks for its cluster best, and the marker outlives
    # REMOVE_TMP because a rerun must not switch what is already switched
    if [ -n "$SWITCH_CONSENSUS_REP" ] && notExists "$TMP/switched.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-pickrepprofile "$TMP/db" "$TMP/clu_accepted/clu_accepted" "$OUT" \
            "$OUT.switched" --threads "$THREADS" -c "$COV" ${PICKREP_PAR} || fail "lin8-pickrepprofile died"
        rm -f "$OUT" "$OUT".[0-9]*
        for _part in "$OUT.switched" "$OUT.switched".[0-9]*; do
            [ -f "$_part" ] || continue
            if [ "$_part" = "$OUT.switched" ]; then mv -f "$_part" "$OUT"; else
                mv -f "$_part" "$OUT.${_part##*.}"
            fi
        done
        mv -f "$OUT.switched.index" "$OUT.index"
        mv -f "$OUT.switched.dbtype" "$OUT.dbtype"
        touch "$TMP/switched.done"
    fi

    # naming is its own pass, and it runs the way linclust names its answer: always
    if notExists "$OUT.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-createtsv "$TMP/db" "$OUT" "$OUT.tsv" ${TSV_PAR} || fail "lin8-createtsv died"
    fi

    if [ "$REPSEQ" -eq 1 ] && notExists "$OUT.rep.fasta"; then
        # shellcheck disable=SC2086
        "$MMSEQS" lin8-createrepseqfasta "$TMP/db" "$OUT" "$OUT.rep.fasta" ${REPSEQ_PAR} \
            || fail "lin8-createrepseqfasta died"
    fi
fi

# what is left here is what a rerun would have needed; every machine clears its own
if [ -n "$REMOVE_TMP" ]; then
    rm -f "$TMP/clu_accepted/clu_accepted.align_assigned_$NODE" "$TMP/kmer/out.$NODE."* "$TMP/pairs/pairs.$NODE."* \
          "$TMP/pref/pref.$NODE."* "$TMP/aln/aln."*".$NODE"* "$TMP/aln/aln.$NODE."* \
          "$TMP/aln/aln_text.$NODE."* "$TMP/lin8clust.$NODE.log"
    if [ "$NODE" -eq 0 ]; then
        rm -rf "$TMP/kmer" "$TMP/pairs" "$TMP/pref" "$TMP/aln" "$TMP/clu_accepted" \
               "$TMP/clu_accepted_plus_redundant"
        # the database stays until the clustering is named, because a cluster is named by rank
        rm -f "$TMP/hash" "$TMP/hash."* "$TMP/lin8clust.sh"
        [ "$REPSEQ" -eq 0 ] && rm -f "$TMP/db" "$TMP/db".[0-9]* "$TMP/db.hist."* "$TMP/db.runs"* \
            "$TMP/db.files" "$TMP/db.dbtype" "$TMP/db_h."*
    fi
fi
echo "machine $NODE finished"
exit 0
