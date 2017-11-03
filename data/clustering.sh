#!/bin/bash -e
# Clustering workflow script
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
    elif [ -d $(dirname "$1") ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}

# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

# processing

INPUT="$(abspath $1)"
TMP_PATH="$(abspath $3)"

if notExists "${TMP_PATH}/aln_redundancy"; then
    $MMSEQS clusthash "$INPUT" "${TMP_PATH}/aln_redundancy" ${DETECTREDUNDANCY_PAR} \
        || fail "Fast filter step $STEP died"
fi

if notExists "${TMP_PATH}/clu_redundancy"; then
    $MMSEQS clust "$INPUT" "${TMP_PATH}/aln_redundancy" "${TMP_PATH}/clu_redundancy" ${CLUSTER_PAR} \
        || fail "Fast Cluster filter step $STEP died"
fi

if notExists "${TMP_PATH}/input_step_redundancy"; then
    $MMSEQS createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" \
        || fail "MMseqs order step $STEP died"
fi

ORIGINAL="$INPUT"
INPUT="${TMP_PATH}/input_step_redundancy"
# call prefilter module
if notExists "${TMP_PATH}/pref"; then
    $RUNNER $MMSEQS prefilter "$INPUT" "$INPUT" "${TMP_PATH}/pref" $PREFILTER_PAR \
        || fail "Prefilter died"
fi

# call alignment module
if notExists "${TMP_PATH}/aln"; then
    $RUNNER $MMSEQS align "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" $ALIGNMENT_PAR \
        || fail "Alignment died"
fi

# call cluster module
if notExists "${TMP_PATH}/clu_step0"; then
    $MMSEQS clust "$INPUT" "${TMP_PATH}/aln" "${TMP_PATH}/clu_step0" $CLUSTER_PAR \
        || fail "Clustering died"
fi

# merge clu_redundancy and clu
if notExists "$2"; then
    $MMSEQS mergeclusters "$ORIGINAL" "$2" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu_step0" \
        || fail "Merging of clusters has died"
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -f "${TMP_PATH}/pref" "${TMP_PATH}/pref.index"
    rm -f "${TMP_PATH}/aln" "${TMP_PATH}/aln.index"
    rm -f "${TMP_PATH}/clu_step0" "${TMP_PATH}/clu_step0.index"
    rm -f "${TMP_PATH}/order_redundancy"
    rm -f "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu_redundancy.index"
    rm -f "${TMP_PATH}/aln_redundancy" "${TMP_PATH}/aln_redundancy.index"
    rm -f "${TMP_PATH}/input_step_redundancy" "${TMP_PATH}/input_step_redundancy.index"
    rm -f "${TMP_PATH}/clustering.sh"
fi
