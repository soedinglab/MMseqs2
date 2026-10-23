#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

if [ -n "$CLUSTHASH" ]; then
    CLUDB="${TMP_PATH}/clu"
else
    CLUDB="$2"
fi

# 0. clusthash
if [ -n "$CLUSTHASH" ]; then
    if notExists "${TMP_PATH}/input_clusthash_clust.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" clusthashfast "$INPUT" "${TMP_PATH}/input_clusthash_clust" ${CLUSTHASHFAST_PAR} \
            || fail "clusthashfast died"
    fi

    INPUT="${TMP_PATH}/input_clusthash_clust_redundancy"
fi

# 1. k-mer matching
if notExists "${TMP_PATH}/pref.dbtype" && notExists "${CLUDB}.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
        || fail "kmermatcher died"
fi

if notExists "${CLUDB}.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align2clust "$INPUT" "${TMP_PATH}/pref" "$CLUDB" ${ALIGN2CLUST_PAR} \
        || fail "align2clust died"
    if [ -n "$REMOVE_TMP" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
        if [ -n "$CLUSTHASH" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_clust_redundancy" ${VERBOSITY}
        fi
    fi
fi

# 2. Merge clusthash pre-clusters with alignment clusters
if [ -n "$CLUSTHASH" ]; then
    if notExists "$2.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" mergeclusters "$SOURCE" "$2" \
            "${TMP_PATH}/input_clusthash_clust" "$CLUDB" ${MERGECLU_PAR} \
            || fail "mergeclusters died"
        if [ -n "$REMOVE_TMP" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_clust" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "$CLUDB" ${VERBOSITY}
        fi
    fi
fi

if [ -n "$REMOVE_TMP" ]; then
    rm -f "${TMP_PATH}/lin8clust.sh"
fi
