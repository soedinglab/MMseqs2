#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref"; then
    # shellcheck disable=SC2086
    "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
        || fail "kmermatcher died"
fi
# 2. Hamming distance pre-clustering
if notExists "${TMP_PATH}/pref_rescore1"; then
    # shellcheck disable=SC2086
    "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/pref_rescore1" ${HAMMING_PAR} \
        || fail "Rescore with hamming distance step died"
fi
if notExists "${TMP_PATH}/pre_clust"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "$INPUT" "${TMP_PATH}/pref_rescore1" "${TMP_PATH}/pre_clust" ${CLUSTER_PAR} \
        || fail "Pre-clustering step died"
fi

awk '{ print $1 }' "${TMP_PATH}/pre_clust.index" > "${TMP_PATH}/order_redundancy"
if notExists "${TMP_PATH}/input_step_redundancy"; then
    "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" \
        || fail "Createsubdb step died"
fi

if notExists "${TMP_PATH}/pref_filter1"; then
    "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${TMP_PATH}/pref" "${TMP_PATH}/pref_filter1" \
        || fail "Createsubdb step died"
fi

if notExists "${TMP_PATH}/pref_filter2"; then
    "$MMSEQS" filterdb "${TMP_PATH}/pref_filter1" "${TMP_PATH}/pref_filter2" --filter-file "${TMP_PATH}/order_redundancy" \
        || fail "Filterdb step died"
fi

INPUT="${TMP_PATH}/input_step_redundancy"
# 3. Ungapped alignment filtering
RESULTDB="${TMP_PATH}/pref_filter2"
if [ -n "$FILTER" ]; then
    if notExists "${TMP_PATH}/pref_rescore2"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_filter2" "${TMP_PATH}/pref_rescore2" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
    fi
    RESULTDB="${TMP_PATH}/pref_rescore2"
fi

# 4. Local gapped sequence alignment.
if notExists "${TMP_PATH}/aln"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$INPUT" "$INPUT" "$RESULTDB" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
        || fail "Alignment step died"
fi

# 5. Clustering using greedy set cover.
if notExists "${TMP_PATH}/clust"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln" "${TMP_PATH}/clust" ${CLUSTER_PAR} \
        || fail "Clustering step died"
fi
if notExists "${TMP_PATH}/clu"; then
    "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/pre_clust" "${TMP_PATH}/clust" \
        || fail "mergeclusters died"
fi

# post processing
mv -f "${TMP_PATH}/clu" "$2" || fail "Could not move result to $2"
mv -f "${TMP_PATH}/clu.index" "$2.index" || fail "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "${TMP_PATH}/pref" "${TMP_PATH}/pref.index"
 rm -f "${TMP_PATH}/pref_rescore1" "${TMP_PATH}/pref_rescore1.index"
 rm -f "${TMP_PATH}/pre_clust" "${TMP_PATH}/pre_clust.index"
 rm -f "${TMP_PATH}/input_step_redundancy" "${TMP_PATH}/input_step_redundancy.index"
 rm -f "${TMP_PATH}/pref_filter1" "${TMP_PATH}/pref_filter1.index"
 rm -f "${TMP_PATH}/pref_filter2" "${TMP_PATH}/pref_filter2.index"
 rm -f "${TMP_PATH}/aln" "${TMP_PATH}/aln.index"
 rm -f "${TMP_PATH}/clust" "${TMP_PATH}/clust.index"

 rm -f "${TMP_PATH}/linclust.sh"
fi
