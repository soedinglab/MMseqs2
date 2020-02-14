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

QUERY="$1"
QUERY_ORF="$1"
if [ -n "$QUERY_NUCL" ]; then
    if notExists "$4/q_orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$1" "$4/q_orfs_aa" ${ORF_PAR} \
            || fail  "extract orfs step died"
    fi
    QUERY="$4/q_orfs_aa"
    QUERY_ORF="$4/q_orfs_aa"
fi

TARGET="$2"
TARGET_ORF="$2"
if [ -n "$TARGET_NUCL" ]; then
if [ -n "$NO_TARGET_INDEX" ]; then
    if notExists "$4/t_orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$2" "$4/t_orfs_aa" ${ORF_PAR} \
            || fail  "extract target orfs step died"
    fi
    TARGET="$4/t_orfs_aa"
    TARGET_ORF="$4/t_orfs_aa"
fi
fi

mkdir -p "$4/search"
if notExists "$4/aln.dbtype"; then
    # shellcheck disable=SC2086
    "$SEARCH" "${QUERY}" "${TARGET}" "$4/aln" "$4/search" \
        || fail "Search step died"
fi

if notExists "$3.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$1" "$QUERY_ORF" "$2" "$TARGET_ORF" "$4/aln"  "$3" ${OFFSETALIGNMENT_PAR} \
        || fail "Offset step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/q_orfs_aa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/t_orfs_aa" ${VERBOSITY}
fi


