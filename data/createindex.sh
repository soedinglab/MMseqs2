#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check amount of input variables
[ "$#" -ne 2 ] && echo "Please provide <sequenceDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -d "$2" ] &&  echo "tmp directory $2 not found!" && mkdir -p "$2";

INPUT="$1"
if [ -n "$TRANSLATED" ]; then
    # 1. extract orf
    if notExists "$2/orfs.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$INPUT" "$2/orfs" $ORF_PAR \
            || fail "extractorfs died"
    fi

    if notExists "$2/orfs_aa.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" translatenucs "$2/orfs" "$2/orfs_aa" $TRANSLATE_PAR \
            || fail "translatenucs died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$2/orfs_aa" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        "$MMSEQS" rmdb "$2/orfs"
        "$MMSEQS" rmdb "$2/orfs_aa"
        rm -f "$2/createindex.sh"
    fi
elif [ -n "$NUCL" ]; then
      # 1. extract orf
    if notExists "$2/nucl_split_seq.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$INPUT" "$2/nucl_split_seq" $SPLIT_SEQ_PAR \
            || fail "splitsequence died"
    fi

    if notExists "$2/nucl_split_seq_rev.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractframes "$2/nucl_split_seq" "$2/nucl_split_seq_rev" $EXTRACT_FRAMES_PAR  \
            || fail "Extractframes died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$2/nucl_split_seq_rev.dbtype" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        rm -f "$2/nucl_split_seq" "$2/nucl_split_seq.index" "$2/nucl_split_seq.dbtype"
        rm -f "$2/nucl_split_seq_rev" "$2/nucl_split_seq_rev.index" "$2/nucl_split_seq_rev.dbtype"
        rm -f "$2/createindex.sh"
    fi
elif [ -n "$LIN_NUCL" ]; then
      # 1. extract orf
    if notExists "$2/nucl_split_seq.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" splitsequence "$INPUT" "$2/nucl_split_seq" $SPLIT_SEQ_PAR \
            || fail "splitsequence died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$2/nucl_split_seq" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        "$MMSEQS" rmdb "$2/nucl_split_seq"
        rm -f "$2/createindex.sh"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" $INDEXER "$INPUT" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"
fi

