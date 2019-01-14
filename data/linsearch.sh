#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}


#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERY="$1"
TARGET="$2"
TMP_PATH="$4"

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmersearch "$QUERY" "$TARGET" "${TMP_PATH}/pref" ${KMERSEARCH_PAR} \
        || fail "kmermatcher died"
fi

# 2. Local gapped sequence alignment.
if notExists "${TMP_PATH}/reverse_aln"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$TARGET" "$QUERY" "${TMP_PATH}/pref" "${TMP_PATH}/reverse_aln" ${ALIGNMENT_PAR} \
        || fail "Alignment step died"
fi

# 3. Local gapped sequence alignment.
if notExists "$3"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" swapresults "$TARGET" "$QUERY" "${TMP_PATH}/reverse_aln" "$3" ${SWAPRESULT_PAR} \
        || fail "Alignment step died"
fi


if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -f "${TMP_PATH}/pref" "${TMP_PATH}/pref.index"
    rm -f "${TMP_PATH}/reverse_aln" "${TMP_PATH}/reverse_aln.index"

    rm -f "${TMP_PATH}/linsearch.sh"
fi
