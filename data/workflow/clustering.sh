#!/bin/sh -e
# Clustering workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"

if notExists "${TMP_PATH}/aln_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clusthash "$INPUT" "${TMP_PATH}/aln_redundancy" ${DETECTREDUNDANCY_PAR} \
        || fail "Fast filter step $STEP died"
fi

if notExists "${TMP_PATH}/clu_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_redundancy" "${TMP_PATH}/clu_redundancy" ${CLUSTER_PAR} \
        || fail "Fast Cluster filter step $STEP died"
fi

if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
        || fail "MMseqs order step $STEP died"
fi

ORIGINAL="$INPUT"
INPUT="${TMP_PATH}/input_step_redundancy"
# call prefilter module
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "$INPUT" "$INPUT" "${TMP_PATH}/pref" $PREFILTER_PAR \
        || fail "Prefilter died"
fi

# call alignment module
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" $ALIGNMENT_PAR \
        || fail "Alignment died"
fi

# call cluster module
if notExists "${TMP_PATH}/clu_step0.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln" "${TMP_PATH}/clu_step0" $CLUSTER_PAR \
        || fail "Clustering died"
fi

# merge clu_redundancy and clu
# shellcheck disable=SC2086
"$MMSEQS" mergeclusters "$ORIGINAL" "$2" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu_step0" $MERGECLU_PAR \
        || fail "Merging of clusters has died"

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_step0" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
    rm -f "${TMP_PATH}/order_redundancy"
    rm -f "${TMP_PATH}/clustering.sh"
fi
