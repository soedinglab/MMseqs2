#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] &&  echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] &&  echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] &&  echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";


QUERY="$1"
TARGET="$2"
TMP_PATH="$4"

if [ -n "$NEEDTARGETSPLIT" ]; then
    if notExists "$TMP_PATH/target_seqs_split.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$TARGET" "$TMP_PATH/target_seqs_split" ${SPLITSEQUENCE_PAR}  \
            || fail "Split sequence died"
    fi
    TARGET="$TMP_PATH/target_seqs_split"
fi

if [ -n "$EXTRACTFRAMES" ]; then
    if notExists "$TMP_PATH/query_seqs.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractframes "$QUERY" "$TMP_PATH/query_seqs" ${EXTRACT_FRAMES_PAR}  \
            || fail "Extractframes died"
    fi
    QUERY="$TMP_PATH/query_seqs"
fi

if [ -n "$NEEDQUERYSPLIT" ]; then
    if notExists "$TMP_PATH/query_seqs_split.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$QUERY" "$TMP_PATH/query_seqs_split" ${SPLITSEQUENCE_PAR}  \
        || fail "Split sequence died"
    fi
    QUERY="$TMP_PATH/query_seqs_split"
fi

mkdir -p "$4/search"
if notExists "$4/aln.dbtype"; then
    # search does not need a parameter because the environment variables will be set by the workflow
    # shellcheck disable=SC2086
    "$SEARCH" "${QUERY}" "${TARGET}" "$4/aln" "$4/search"  \
        || fail "Search step died"
fi

if notExists "$3.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$1" "${QUERY}" "$2" "${TARGET}" "$4/aln"  "$3" ${OFFSETALIGNMENT_PAR} \
        || fail "Offset step died"
fi


if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/q_orfs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/q_orfs_aa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/t_orfs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "$4/t_orfs_aa" ${VERBOSITY}
fi

