#!/bin/sh -e
# reciprocal best hit workflow
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

A_DB="$1"
B_DB="$2"
RBH_RES="$3"
TMP_PATH="$4"

# search in both directions:
if [ ! -e "${TMP_PATH}/resAB.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${A_DB}" "${B_DB}" "${TMP_PATH}/resAB" "${TMP_PATH}/tempAB" ${SEARCH_A_B_PAR} \
        || fail "search A vs. B died"
fi

if [ ! -e "${TMP_PATH}/resBA.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${B_DB}" "${A_DB}" "${TMP_PATH}/resBA" "${TMP_PATH}/tempBA" ${SEARCH_B_A_PAR} \
        || fail "search B vs. A died"
fi


# sort A->B by decreasing bitscores:
if [ ! -e "${TMP_PATH}/resAB_sorted.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/resAB" "${TMP_PATH}/resAB_sorted" --sort-entries 2 --filter-column 2 ${THREADS_COMP_PAR} \
        || fail "sort resAB by bitscore died"
fi

# extract a single best hit in A->B direction (used to take best bitscore for A):
if [ ! -e "${TMP_PATH}/resA_best_B.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/resAB_sorted" "${TMP_PATH}/resA_best_B" --extract-lines 1 ${THREADS_COMP_PAR} \
        || fail "extract A best B died"
fi

# extract best hit(s) in B->A direction:
if [ ! -e "${TMP_PATH}/resB_best_A.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/resBA" "${TMP_PATH}/resB_best_A" --beats-first --filter-column 2 --comparison-operator e ${THREADS_COMP_PAR} \
        || fail "extract B best A died"
fi

# swap the direction of resB_best_A:
if [ ! -e "${TMP_PATH}/resB_best_A_swap.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${B_DB}" "${A_DB}" "${TMP_PATH}/resB_best_A" "${TMP_PATH}/resB_best_A_swap" ${THREADS_COMP_PAR} -e 100000000 \
        || fail "swap B best A died"
fi

# merge the best results:
if [ ! -e "${TMP_PATH}/res_best_merged.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" mergedbs "${TMP_PATH}/resA_best_B" "${TMP_PATH}/res_best_merged" "${TMP_PATH}/resA_best_B" "${TMP_PATH}/resB_best_A_swap" ${VERB_COMP_PAR} \
        || fail "merge best hits died"
fi

# sort by bitscore (decreasing order):
if [ ! -e "${TMP_PATH}/res_best_merged_sorted.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/res_best_merged" "${TMP_PATH}/res_best_merged_sorted" --sort-entries 2 --filter-column 2 ${THREADS_COMP_PAR} \
        || fail "sort by bitscore died"
fi

# identify the RBH pairs and write them to a result db:
if [ ! -e "${RBH_RES}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2rbh "${TMP_PATH}/res_best_merged_sorted" "${RBH_RES}" ${THREADS_COMP_PAR} \
        || fail "result2rbh died"
fi

if [ -n "$REMOVE_TMP" ]; then
    rm -rf "${TMP_PATH}/tempAB"
    rm -rf "${TMP_PATH}/tempBA"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resAB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resBA" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resA_best_B" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resB_best_A" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/res_best" ${VERBOSITY}
    rm -f "${TMP_PATH}/rbh.sh"
fi
