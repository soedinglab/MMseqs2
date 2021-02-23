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
if notExists "${PROFILEDB}.dbtype"; then
    # symlink the profile DB that can be reduced at every iteration the search
    ln -s "${TARGET}" "${PROFILEDB}"
    ln -s "${TARGET}.dbtype" "${PROFILEDB}.dbtype"
    cp -f "${TARGET}.index" "${PROFILEDB}.index"

    echo "${AVAIL_DISK}" > "${PROFILEDB}.meta"
else
    read -r AVAIL_DISK < "${PROFILEDB}.meta"
fi

TOTAL_NUM_PROFILES=$(wc -l < "${PROFILEDB}.index")
NUM_SEQS_THAT_SATURATE="$(wc -l < "${INPUT}.index")"
FIRST_INDEX_LINE=1
NUM_PROFS_IN_STEP=1
NUM_PREF_RESULTS_IN_ALL_PREV_STEPS=0

STEP=0

while [ "${FIRST_INDEX_LINE}" -le "${TOTAL_NUM_PROFILES}" ]; do
    if [ -f "${TMP_PATH}/aln_${STEP}.checkpoint" ]; then
        # restore values from previous run, in case it was aborted
        read -r FIRST_INDEX_LINE NUM_PREF_RESULTS_IN_ALL_PREV_STEPS < "${TMP_PATH}/aln_${STEP}.checkpoint"
    fi

    # predict NUM_SEQS_THAT_SATURATE as the average number of prefilter results per profile in previous steps
    # this allows one to increase NUM_PROFS_IN_STEP
    if [ "${NUM_PREF_RESULTS_IN_ALL_PREV_STEPS}" -gt 0 ]; then
        # BE MORE CAUTIOUS?
        NUM_PROFS_PROCESSED="$((FIRST_INDEX_LINE-1))"
        NUM_SEQS_THAT_SATURATE="$((NUM_PREF_RESULTS_IN_ALL_PREV_STEPS/NUM_PROFS_PROCESSED))"
        # this also assures no division by 0 later on
        if [ "${NUM_SEQS_THAT_SATURATE}" -lt 1 ]; then
            NUM_SEQS_THAT_SATURATE=1
        fi
    fi
    # prefilter res size (10 + 1 + 6 + 1 + 3 + 1) + 3 byte buffer
    RESSIZE=25
    # disk usage allowance not set by the user (i.e. AVAIL_DISK = 0), compute it for optimal usage
    if [ "${AVAIL_DISK}" -eq 0 ]; then
        CURRENT_AVAIL_DISK_SPACE=$(($("$MMSEQS" diskspaceavail "${TMP_PATH}")/2))
        # Compute the max number of profiles that can be processed
        # based on the number of hits that saturate
        NUM_PROFS_IN_STEP="$((CURRENT_AVAIL_DISK_SPACE/NUM_SEQS_THAT_SATURATE/RESSIZE))"
    else
        NUM_PROFS_IN_STEP="$((AVAIL_DISK/NUM_SEQS_THAT_SATURATE/RESSIZE))"
    fi

    # no matter what, process at least one profile...
    if [ "${NUM_PROFS_IN_STEP}" -lt 1 ]; then
        NUM_PROFS_IN_STEP=1
    fi

    # take a chunk of profiles from FIRST_INDEX_LINE to (FIRST_INDEX_LINE + NUM_PROFS_IN_STEP -1)
    LAST_INDEX_LINE_TO_PROCESS="$((FIRST_INDEX_LINE+NUM_PROFS_IN_STEP-1))"
    awk -v first="${FIRST_INDEX_LINE}" -v last="${LAST_INDEX_LINE_TO_PROCESS}" \
        "NR >= first && NR <= last { print; }" "${TARGET}.index" > "${PROFILEDB}.index"

    # prefilter current chunk
    if notExists "${TMP_PATH}/pref.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" prefilter "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
            || fail "prefilter died"
        touch "${TMP_PATH}/pref.done"
    fi

    # count the number of hits for all profiles in current step chunk
    if notExists "${TMP_PATH}/pref_count.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" result2stats "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/pref_count.tsv" \
            --stat linecount --tsv ${THREADS_COMP_PAR} \
            || fail "result2stats died"
    fi
    # update the starting point for the next step and the total number of pref results
    NUM_PREF_RESULTS_IN_STEP=$(awk '{sum+=$1;} END{print sum;}' "${TMP_PATH}/pref_count.tsv")
    rm -f "${TMP_PATH}/pref_count.tsv"

    NUM_PREF_RESULTS_IN_ALL_PREV_STEPS="$((NUM_PREF_RESULTS_IN_ALL_PREV_STEPS+NUM_PREF_RESULTS_IN_STEP))"
    FIRST_INDEX_LINE="$((FIRST_INDEX_LINE+NUM_PROFS_IN_STEP))"

    # align current step chunk
    if notExists "${TMP_PATH}/aln.done"; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" "${ALIGN_MODULE}" "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/aln"  ${ALIGNMENT_IT_PAR} \
            || fail "align died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
        touch "${TMP_PATH}/aln.done"
    fi


        # no matter what, process at least one profile...
    if [ "${FILTER_RESULT}" -eq 1 ]; then
        # shellcheck disable=SC2086
        ${RUNNER} "$MMSEQS" filterresult "${PROFILEDB}" "${INPUT}" "${TMP_PATH}/aln" "${TMP_PATH}/aln_filt"  ${FILTER_PAR} \
            || fail "align died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY} || fail "rmdb aln died"
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/aln_filt" "${TMP_PATH}/aln" ${VERBOSITY} || fail "mv aln_filt aln died"
        touch "${TMP_PATH}/aln.done"
    fi


    # merge swapped alignment of current chunk to previous steps
    if [ -f "${TMP_PATH}/aln_merged.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "${TARGET}" "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" "${TMP_PATH}/aln" ${VERBOSITY} \
            || fail "mergedbs died"
        # rmdb of aln_merged to avoid conflict with unmerged dbs: aln_merged.0, .1...
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_merged" ${VERBOSITY} || fail "rmdb aln_merged died"
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/aln_merged_new" "${TMP_PATH}/aln_merged" ${VERBOSITY} || fail "mv aln_merged_new aln_merged died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY} || fail "rmdb aln died"
    else
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/aln" "${TMP_PATH}/aln_merged" ${VERBOSITY} \
            || fail "mvdb died"
    fi

    STEP="$((STEP+1))"
    # update for the next step
    rm -f "${TMP_PATH}/pref.done" "${TMP_PATH}/aln.done"
    printf "%d\\t%s\\n" "${FIRST_INDEX_LINE}" "${NUM_PREF_RESULTS_IN_ALL_PREV_STEPS}" > "${TMP_PATH}/aln_${STEP}.checkpoint"

done


# swap alignment of current step chunk
if notExists "${TMP_PATH}/aln.done"; then
    # keep only the top max-seqs hits according to the default alignment sorting criteria
    # shellcheck disable=SC2086
    "$MMSEQS" align "${TARGET}" "${INPUT}" "${TMP_PATH}/aln_merged" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
       || fail "sortresult died"
    # rmdb of aln_merged to avoid conflict with unmerged dbs: aln_merged.0, .1...
       # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_merged" ${VERBOSITY} || fail "rmdb aln_merged died"
fi

# note: the evalue has been corrected for inverted search by the workflow caller
# shellcheck disable=SC2086
"$MMSEQS" swapresults  "${TARGET}" "${INPUT}" "${TMP_PATH}/aln" "${RESULT}" ${SWAPRES_PAR} \
    || fail "swapresults died"


if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_merged" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${PROFILEDB}" ${VERBOSITY}
    CURR_STEP=0
    while [ "${CURR_STEP}" -le "${STEP}" ]; do
        if [ -f "${TMP_PATH}/aln_${CURR_STEP}.checkpoint" ]; then
            rm -f "${TMP_PATH}/aln_${CURR_STEP}.checkpoint"
        fi
        CURR_STEP="$((CURR_STEP+1))"
    done
    rm -f "${PROFILEDB}.meta"
    rm -f "$TMP_PATH/searchslicedtargetprofile.sh"
fi
