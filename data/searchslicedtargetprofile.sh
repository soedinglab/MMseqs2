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

hasCommand awk
hasCommand wc
hasCommand join
hasCommand sort

# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <resultDB> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "TMP directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$2"
RESULTS="$3"
TMP_PATH="$4"

PROFILEDB="${TMP_PATH}/profileDB"
if notExists "${TMP_PATH}/profileDB"; then
    # symlink the profile DB that can be reduced at every iteration the search
    ln -s "${TARGET}" "${PROFILEDB}"
    ln -s "${TARGET}.dbtype" "${PROFILEDB}.dbtype"
    sort -k1,1 "${TARGET}.index" > "${PROFILEDB}.index"

    echo "${AVAIL_MEM}" > "${PROFILEDB}.params"
else
    read -r AVAIL_MEM < "${PROFILEDB}.params"
fi

# shellcheck disable=SC2046,SC2005
NUM_PROFILES="$(echo $(wc -l < "${PROFILEDB}.index"))"

OFFSET=0
STEP=0
# MAX_STEPS is set by the workflow
# shellcheck disable=SC2153
while [ "$STEP" -lt "$MAX_STEPS" ] && [ "$NUM_PROFILES" -gt 0 ]; do
    if [ ! -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
        # Compute the max number of sequence that are reasonable to swap according to the number of profiles
        # 90 bytes/query-result line max.
        MAX_SEQS="$((1024*(AVAIL_MEM/NUM_PROFILES)/90))"
    else
        # restore values from previous run, in case it was aborted
        read -r MAX_SEQS OFFSET < "${TMP_PATH}/aln_${STEP}.checkpoint"
    fi

    # Max result to go into the prefilter list
    SEARCH_LIM="$((OFFSET+MAX_SEQS))"

    if notExists "${TMP_PATH}/pref_${OFFSET}"; then
        # shellcheck disable=SC2086
        "$MMSEQS" prefilter "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref_${OFFSET}" \
            --max-seqs "${SEARCH_LIM}" --offset-result "${OFFSET}" ${PREFILTER_PAR} \
             || fail "prefilter died"
    fi

    if notExists "${TMP_PATH}/pref_swapped_${OFFSET}"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapresults "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref_${OFFSET}" "${TMP_PATH}/pref_swapped_${OFFSET}" ${SWAP_PAR} \
             || fail "swapresults died"
    fi

    if notExists "${TMP_PATH}/aln_${OFFSET}"; then
        # note: we recover the right evalue, since it is computed according to the original target db
        # shellcheck disable=SC2086
        "$MMSEQS" align "${INPUT}" "${TARGET}" "${TMP_PATH}/pref_swapped_${OFFSET}" "${TMP_PATH}/aln_${OFFSET}" ${ALIGNMENT_PAR} \
             || fail "align died"
    fi

    MERGED="${TMP_PATH}/aln_${OFFSET}"
    if [ -f "${TMP_PATH}/aln_merged" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "${INPUT}" "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_${OFFSET}" "${TMP_PATH}/aln_merged" ${THREADS_PAR} \
            || fail "mergedbs died"

        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" \
            --sort-entries 2 --filter-column 2 ${THREADS_PAR} \
            || fail "filterdb (sorting) died"

        rm -f "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged_new"
        MERGED="${TMP_PATH}/aln_merged"
    fi

    # keep only the top MAX_RESULTS_PER_QUERY hits according to evalue
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${MERGED}" "${TMP_PATH}/aln_merged_trunc" \
        --extract-lines "$MAX_RESULTS_PER_QUERY" ${THREADS_PAR} \
        || fail "filterdb (extract lines) died"
    mv -f "${TMP_PATH}/aln_merged_trunc" "${TMP_PATH}/aln_merged"
    mv -f "${TMP_PATH}/aln_merged_trunc" "${TMP_PATH}/aln_merged.index"

    # remove profiles that do not provide more hits
    # shellcheck disable=SC2086
    "$MMSEQS" result2stats "$PROFILEDB" "$INPUT" "${TMP_PATH}/pref_${OFFSET}"  "${TMP_PATH}/aln_merged_count" \
        --stat linecount ${THREADS_PAR} \
        || fail "result2stats died"

    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/aln_merged_count" "${TMP_PATH}/aln_merged_keep" \
        --filter-column 1 --comparison-operator ge --comparison-value "${MAX_SEQS}" ${THREADS_PAR} \
        || fail "filterdb died"

    awk '$3 > 1 { print $1 }' "${TMP_PATH}/aln_merged_keep.index" | sort -k1,1 > "${TMP_PATH}/aln_merged_keep.list"
    join "${TMP_PATH}/aln_merged_keep.list" "${PROFILEDB}.index" > "${PROFILEDB}.index.tmp"
    mv -f "${PROFILEDB}.index.tmp" "${PROFILEDB}.index"

    # keep for the prefilter only the next hits
    OFFSET="$SEARCH_LIM"

    # shellcheck disable=SC2046,SC2005
    NUM_PROFILES="$(echo $(wc -l < "${PROFILEDB}.index"))"

    printf "%d\\t%d\\n" "${MAX_SEQS}" "${OFFSET}" > "${TMP_PATH}/aln_${STEP}.checkpoint"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove intermediate temporary files"
        rm -f "${TMP_PATH}/aln_merged_keep" "${TMP_PATH}/aln_merged_keep.index" "${TMP_PATH}/aln_merged_keep.list"
        rm -f "${TMP_PATH}/aln_merged_count" "${TMP_PATH}/aln_merged_count.index"
    fi
    STEP="$((STEP+1))"
done

# Save the results
mv -f "${TMP_PATH}/aln_merged" "$RESULTS"
mv -f "${TMP_PATH}/aln_merged.index" "$RESULTS.index"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -f "$TMP_PATH/searchslicedtargetprofile.sh"
fi


