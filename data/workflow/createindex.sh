#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 2 ] && echo "Please provide <sequenceDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -d "$2" ] && echo "tmp directory $2 not found!" && mkdir -p "$2";

INPUT="$1"
if [ -n "$TRANSLATED" ]; then
    # 1. extract orf
    if notExists "$2/orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$INPUT" "$2/orfs_aa" ${ORF_PAR} \
            || fail "extractorfs died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$2/orfs_aa" "$INPUT" ${INDEX_PAR} \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$2/orfs_aa" ${VERBOSITY}
        rm -f "$2/createindex.sh"
    fi
elif [ -n "$LIN_NUCL" ] || [ -n "$NUCL" ]; then
      # 1. extract orf
    if notExists "$2/nucl_split_seq.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$INPUT" "$2/nucl_split_seq" ${SPLIT_SEQ_PAR} \
            || fail "splitsequence died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$2/nucl_split_seq" "$INPUT" ${INDEX_PAR} \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$2/nucl_split_seq" ${VERBOSITY}
        rm -f "$2/createindex.sh"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$INPUT" "$INPUT" ${INDEX_PAR} \
        || fail "indexdb died"
fi

