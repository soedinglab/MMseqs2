#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
RESULTS="$3"
TMP_PATH="$4"

# call prefilter module
if notExists "${TMP_PATH}/pref.dbtype"; then
     # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "${INPUT}" "${2}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Prefilter died"
fi

if notExists "${TMP_PATH}/pref_swapped.dbtype"; then
     # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${INPUT}" "${2}" "${TMP_PATH}/pref" "${TMP_PATH}/pref_swapped" ${SWAP_PAR} \
        || fail "Swapresults pref died"
fi

# call alignment module
if notExists "$TMP_PATH/aln_swapped.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${2}" "${INPUT}" "${TMP_PATH}/pref_swapped" "${TMP_PATH}/aln_swapped" ${ALIGNMENT_PAR} \
        || fail "Alignment died"
fi

if notExists "${RESULTS}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${2}" "${INPUT}" "${TMP_PATH}/aln_swapped"  "${RESULTS}" ${SWAP_PAR} \
        || fail "Swapresults aln died"
fi


if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref_swapped" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_swapped" ${VERBOSITY}
    rm -f "${TMP_PATH}/searchtargetprofile.sh"
fi
