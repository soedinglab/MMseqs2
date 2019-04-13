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
hasCommand join
hasCommand sort

# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <resultDB> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3.dbtype" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "TMP directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$(abspath "$2")"
RESULT="$3"
TMP_PATH="$4"

PROFILEDB="${TMP_PATH}/profileDB"
if notExists "${TMP_PATH}/profileDB"; then
    # symlink the profile DB that can be reduced at every iteration the search
    ln -s "${TARGET}" "${PROFILEDB}"
    ln -s "${TARGET}.dbtype" "${PROFILEDB}.dbtype"
    sort -k1,1 "${TARGET}.index" > "${PROFILEDB}.index"

    echo "${AVAIL_DISK}" > "${PROFILEDB}.meta"
else
    read -r AVAIL_DISK < "${PROFILEDB}.meta"
fi

# echo call to trim whitespace wc produces
# shellcheck disable=SC2046,SC2005
NUM_PROFILES="$(echo $(wc -l < "${PROFILEDB}.index"))"

OFFSET=0
STEP=0
# MAX_STEPS is set by the workflow
# shellcheck disable=SC2153
while [ "${STEP}" -lt "${MAX_STEPS}" ] && [ "${NUM_PROFILES}" -gt 0 ]; do
    if [ -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
        # restore values from previous run, in case it was aborted
        read -r NUM_PROFILES OFFSET < "${TMP_PATH}/aln_${STEP}.checkpoint"
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

    # Max result to go into the prefilter list
    SEARCH_LIM="$((OFFSET+MAX_SEQS))"

    if notExists "${TMP_PATH}/pref.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" prefilter "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" \
            --max-seqs "${SEARCH_LIM}" --offset-result "${OFFSET}" ${PREFILTER_PAR} \
            || fail "prefilter died"
        touch "${TMP_PATH}/pref.done"
    fi

    if notExists "${TMP_PATH}/pref_count" || notExists "${TMP_PATH}/pref_keep.list"; then
        # shellcheck disable=SC2086
        "$MMSEQS" result2stats "$PROFILEDB" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/pref_count" \
            --stat linecount ${THREADS_PAR} \
            || fail "result2stats died"
    fi

    if notExists "${TMP_PATH}/pref_keep.list"; then
        # remove profiles that do not provide more hits
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/pref_count" "${TMP_PATH}/pref_keep" \
            --filter-column 1 --comparison-operator ge --comparison-value "${MAX_SEQS}" ${THREADS_PAR} \
            || fail "filterdb died"
        "$MMSEQS" rmdb "${TMP_PATH}/pref_count"
        awk '$3 > 1 { print $1 }' "${TMP_PATH}/pref_keep.index" | sort -k1,1 > "${TMP_PATH}/pref_keep.list"
         "$MMSEQS" rmdb "${TMP_PATH}/pref_keep"
    fi

    if notExists "${TMP_PATH}/aln.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" align "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "align died"
        "$MMSEQS" rmdb "${TMP_PATH}/pref"
        touch "${TMP_PATH}/aln.done"
    fi

    if notExists "${TMP_PATH}/aln_swap.done"; then
        # note: the evalue has been corrected for inverted search by MMseqs wrapper
        # shellcheck disable=SC2086
        "$MMSEQS" swapresults "${TARGET}" "${INPUT}" "${TMP_PATH}/aln" "${TMP_PATH}/aln_swap" ${SWAP_PAR} \
            || fail "swapresults died"
        "$MMSEQS" rmdb "${TMP_PATH}/aln"
        touch "${TMP_PATH}/aln_swap.done"
    fi

    MERGED="${TMP_PATH}/aln_swap"
    if [ -f "${TMP_PATH}/aln_merged" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "${INPUT}" "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" "${TMP_PATH}/aln_swap" ${VERBOSITY_PAR} \
            || fail "mergedbs died"
        "$MMSEQS" mvdb "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged"
        "$MMSEQS" rmdb "${TMP_PATH}/aln_swap"
        MERGED="${TMP_PATH}/aln_merged"
    fi

    # keep only the top max-seqs hits according to the default alignment sorting criteria
    # shellcheck disable=SC2086
    "$MMSEQS" sortresult "${MERGED}" "${TMP_PATH}/aln_merged_trunc" ${SORTRESULT_PAR} \
        || fail "sortresult died"

    "$MMSEQS" mvdb "${TMP_PATH}/aln_merged_trunc" "${TMP_PATH}/aln_merged"
    
    join "${TMP_PATH}/pref_keep.list" "${PROFILEDB}.index" > "${PROFILEDB}.index.tmp"
    mv -f "${PROFILEDB}.index.tmp" "${PROFILEDB}.index"

    # keep for the prefilter only the next hits
    OFFSET="$SEARCH_LIM"
    # shellcheck disable=SC2046,SC2005
    NUM_PROFILES="$(echo $(wc -l < "${PROFILEDB}.index"))"
    rm -f "${TMP_PATH}/pref.done" "${TMP_PATH}/aln.done" "${TMP_PATH}/pref_keep.list"
    "$MMSEQS" rmdb "${TMP_PATH}/aln_swap"
    rm -f "${TMP_PATH}/aln_swap.done"
    printf "%d\\t%d\\n" "${NUM_PROFILES}" "${OFFSET}" > "${TMP_PATH}/aln_${STEP}.checkpoint"

    STEP="$((STEP+1))"
done

"$MMSEQS" mvdb "${TMP_PATH}/aln_merged" "${RESULT}"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    STEP=0
    while [ "${STEP}" -lt "${MAX_STEPS}" ]; do
        if [ -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
            "$MMSEQS" rmdb "${TMP_PATH}/aln_${STEP}.checkpoint"
        fi
        STEP="$((STEP+1))"
    done
    "$MMSEQS" rmdb "${TMP_PATH}/profileDB"
    rm -f "${TMP_PATH}/profileDB.meta"
    rm -f "$TMP_PATH/searchslicedtargetprofile.sh"
fi


