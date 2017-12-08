#!/bin/bash -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}

# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$(abspath "$1")"
TARGET="$(abspath "$2")"
RESULTS="$(abspath "$3")"
TMP_PATH="$(abspath "$4")"

# call prefilter module
if notExists "${TMP_PATH}/pref"; then
    $RUNNER $MMSEQS prefilter   "${INPUT}" "${2}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Prefilter died"
fi

if notExists "${TMP_PATH}/pref_swapped"; then
    $MMSEQS swapresults "${INPUT}" "${2}" "${TMP_PATH}/pref" "${TMP_PATH}/pref_swapped" ${SWAP_PAR} \
        || fail "Swapresults pref died"
fi

# call alignment module
if notExists "$TMP_PATH/aln_swapped"; then
    $RUNNER $MMSEQS align       "${2}" "${INPUT}" "${TMP_PATH}/pref_swapped" "${TMP_PATH}/aln_swapped" ${ALIGNMENT_PAR} \
        || fail "Alignment died"
fi

if notExists "$TMP_PATH/aln"; then
    $MMSEQS swapresults "${2}" "${INPUT}" "${TMP_PATH}/aln_swapped"  "${TMP_PATH}/aln" ${SWAP_PAR} \
        || fail "Swapresults aln died"
fi

# post processing
(mv -f "${TMP_PATH}/aln" "${RESULTS}"; mv -f "${TMP_PATH}/aln.index" "${RESULTS}.index") || fail "Could not move result to ${RESULTS}"

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -f "${TMP_PATH}/pref" "${TMP_PATH}/pref.index"
    rm -f "${TMP_PATH}/pref_swapped" "${TMP_PATH}/pref_swapped.index"
    rm -f "${TMP_PATH}/aln_swapped" "${TMP_PATH}/aln_swapped.index"
    rm -f "${TMP_PATH}/searchtargetprofile.sh"
fi
