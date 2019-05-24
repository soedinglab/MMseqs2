#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

hasCommand awk
hasCommand wc

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <resultDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "TMP directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$(abspath "$2")"
RESULT="$3"
TMP_PATH="$4"

PROFILEDB="${TMP_PATH}/profileDB"
if notExists "${PROFILEDB}.index"; then
    # symlink the profile DB that can be reduced at every iteration the search
    ln -s "${TARGET}" "${PROFILEDB}"
    ln -s "${TARGET}.dbtype" "${PROFILEDB}.dbtype"
    cp -f "${TARGET}.index" "${PROFILEDB}.index"

    echo "${AVAIL_DISK}" > "${PROFILEDB}.meta"
else
    read -r AVAIL_DISK < "${PROFILEDB}.meta"
fi

NUM_PROFILES=$(wc -l < "${PROFILEDB}.index")

PREV_MAX_SEQS=""
STEP=0
# MAX_STEPS is set by the workflow
# shellcheck disable=SC2153
while [ "${STEP}" -lt "${MAX_STEPS}" ] && [ "${NUM_PROFILES}" -gt 0 ]; do
    if [ -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
        # restore values from previous run, in case it was aborted
        read -r NUM_PROFILES PREV_MAX_SEQS < "${TMP_PATH}/aln_${STEP}.checkpoint"
        continue
    fi

    # Disk usage allowance not set by the user (i.e. AVAIL_DISK = 0),
    # Compute it for optimal usage
    if [ "${AVAIL_DISK}" -eq 0 ]; then
        CURRENT_AVAIL_DISK_SPACE=$(($("$MMSEQS" diskspaceavail "${TMP_PATH}")/2))
        # Compute the max number of sequence according to the number of profiles
        # 90 bytes/query-result line max.
        MAX_SEQS="$((1024*CURRENT_AVAIL_DISK_SPACE/NUM_PROFILES/90))"
    else
        MAX_SEQS="$((1024*AVAIL_DISK/NUM_PROFILES/90))"
    fi

    if notExists "${TMP_PATH}/pref.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" prefilter "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" \
            --max-seqs "${MAX_SEQS}" --prev-max-seqs "${PREV_MAX_SEQS}" ${PREFILTER_PAR} \
            || fail "prefilter died"
        if [ "${PREV_MAX_SEQS}" = "" ]; then
            PREV_MAX_SEQS="${MAX_SEQS}"
        else
            PREV_MAX_SEQS="${PREV_MAX_SEQS},${MAX_SEQS}"
        fi
        touch "${TMP_PATH}/pref.done"
    fi

    if notExists "${TMP_PATH}/pref_count.index" || notExists "${TMP_PATH}/pref_keep.list"; then
        # shellcheck disable=SC2086
        "$MMSEQS" result2stats "$PROFILEDB" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/pref_count" \
            --stat linecount ${THREADS_COMP_PAR} \
            || fail "result2stats died"
    fi

    if notExists "${TMP_PATH}/pref_keep.list"; then
        # remove profiles that do not provide more hits
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/pref_count" "${TMP_PATH}/pref_keep" \
            --filter-column 1 --comparison-operator ge --comparison-value "${MAX_SEQS}" ${THREADS_COMP_PAR} \
            || fail "filterdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_count" ${VERBOSITY_PAR}
        awk '$3 > 1 { print $1 }' "${TMP_PATH}/pref_keep.index" > "${TMP_PATH}/pref_keep.list"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_keep"
    fi

    if notExists "${TMP_PATH}/aln.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" align "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "align died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY_PAR}
        touch "${TMP_PATH}/aln.done"
    fi

    if notExists "${TMP_PATH}/aln_swap.done"; then
        # note: the evalue has been corrected for inverted search by the workflow caller
        # shellcheck disable=SC2086
        "$MMSEQS" swapresults "${TARGET}" "${INPUT}" "${TMP_PATH}/aln" "${TMP_PATH}/aln_swap" ${SWAP_PAR} \
            || fail "swapresults died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY_PAR}
        touch "${TMP_PATH}/aln_swap.done"
    fi

    MERGED="${TMP_PATH}/aln_swap"
    if [ -f "${TMP_PATH}/aln_merged.index" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "${INPUT}" "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" "${TMP_PATH}/aln_swap" ${VERBOSITY_PAR} \
            || fail "mergedbs died"
        # shellcheck disable=SC2086
        # rmdb of aln_merged to avoid conflict with unmerged dbs: aln_merged.0, .1...
        "$MMSEQS" rmdb "${TMP_PATH}/aln_merged" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_swap" ${VERBOSITY_PAR}
        MERGED="${TMP_PATH}/aln_merged"
    fi

    # keep only the top max-seqs hits according to the default alignment sorting criteria
    # shellcheck disable=SC2086
    "$MMSEQS" sortresult "${MERGED}" "${TMP_PATH}/aln_merged_trunc" ${SORTRESULT_PAR} \
        || fail "sortresult died"

    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${TMP_PATH}/aln_merged_trunc" "${TMP_PATH}/aln_merged" ${VERBOSITY_PAR}
    
    awk 'NR == FNR { f[$1] = 0; next; } $1 in f { print; }' "${TMP_PATH}/pref_keep.list" "${PROFILEDB}.index" > "${PROFILEDB}.index.tmp"
    mv -f "${PROFILEDB}.index.tmp" "${PROFILEDB}.index"

    # shellcheck disable=SC2046
    NUM_PROFILES=$(wc -l < "${PROFILEDB}.index")
    rm -f "${TMP_PATH}/pref.done" "${TMP_PATH}/aln.done" "${TMP_PATH}/pref_keep.list"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_swap" ${VERBOSITY_PAR}
    rm -f "${TMP_PATH}/aln_swap.done"
    printf "%d\\t%s\\n" "${NUM_PROFILES}" "${PREV_MAX_SEQS}" > "${TMP_PATH}/aln_${STEP}.checkpoint"

    STEP="$((STEP+1))"
done

# shellcheck disable=SC2086
"$MMSEQS" mvdb "${TMP_PATH}/aln_merged" "${RESULT}" ${VERBOSITY_PAR}

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    STEP=0
    while [ "${STEP}" -lt "${MAX_STEPS}" ]; do
        if [ -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
            rm -f "${TMP_PATH}/aln_${STEP}.checkpoint"
        fi
        STEP="$((STEP+1))"
    done
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${PROFILEDB}" ${VERBOSITY_PAR}
    rm -f "${PROFILEDB}.meta"
    rm -f "$TMP_PATH/searchslicedtargetprofile.sh"
fi


