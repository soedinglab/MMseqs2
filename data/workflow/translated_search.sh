#!/bin/sh -e
# Translated search workflow
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
[ ! -f "$1.dbtype" ] && echo "$1 not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2 not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

TMP_PATH="$4"
QUERY="$1"
QUERY_ORF="$1"
if [ -n "$QUERY_NUCL" ]; then
    if notExists "${TMP_PATH}/q_orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$1" "${TMP_PATH}/q_orfs_aa" ${ORF_PAR} \
            || fail  "extract orfs step died"
    fi
    QUERY="${TMP_PATH}/q_orfs_aa"
    QUERY_ORF="${TMP_PATH}/q_orfs_aa"
fi

TARGET="$2"
TARGET_ORF="$2"
if [ -n "$TARGET_NUCL" ]; then
if [ -n "$NO_TARGET_INDEX" ]; then
    if notExists "${TMP_PATH}/t_orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$2" "${TMP_PATH}/t_orfs_aa" ${ORF_PAR} \
            || fail  "extract target orfs step died"
    fi
    TARGET="${TMP_PATH}/t_orfs_aa"
    TARGET_ORF="${TMP_PATH}/t_orfs_aa"
fi
fi

if [ -n "$QUERY_NUCL" ] && [ -n "${ORF_FILTER}" ]; then
    if notExists "${TMP_PATH}/q_orfs_aa_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" prefilter "${QUERY}" "${TARGET}" "${TMP_PATH}/q_orfs_aa_pref" --min-ungapped-score 3 -s 3 -k 6 --diag-score 0 --spaced-kmer-mode 0 --max-seqs 1 ${THREAD_COMP_PAR} \
            || fail "Reference search died"
    fi

    if notExists "${TMP_PATH}/q_orfs_aa_filter.dbtype"; then
        awk '$3 > 1 { print $1 }' "${TMP_PATH}/q_orfs_aa_pref.index" > "${TMP_PATH}/q_orfs_aa_filter.list"
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/q_orfs_aa_filter.list" "${QUERY}" "${TMP_PATH}/q_orfs_aa_filter" ${CREATESUBDB_PAR} \
            || fail "createsubdb died"
    fi
    QUERY="${TMP_PATH}/q_orfs_aa_filter"
    QUERY_ORF="${TMP_PATH}/q_orfs_aa_filter"
fi

mkdir -p "${TMP_PATH}/search"
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    "$SEARCH" "${QUERY}" "${TARGET}" "${TMP_PATH}/aln" "${TMP_PATH}/search" \
        || fail "Search step died"
fi

if notExists "$3.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$1" "$QUERY_ORF" "$2" "$TARGET_ORF" "${TMP_PATH}/aln"  "$3" ${OFFSETALIGNMENT_PAR} \
        || fail "Offset step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/q_orfs_aa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/t_orfs_aa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
    if [ -n "${REFERENCE_FILTER}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/q_orfs_aa_pref" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/q_orfs_aa_filter" ${VERBOSITY}
        rm -f "${TMP_PATH}/q_orfs_aa_filter.list"
    fi
    rm -f "${TMP_PATH}/translated_search.sh"
fi
