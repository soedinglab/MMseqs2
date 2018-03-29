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
if [ -n "$NUCL" ]; then
    # 1. extract orf
    if notExists "$2/orfs"; then
        # shellcheck disable=SC2086
        "$MMSEQS" extractorfs "$INPUT" "$2/orfs" $ORF_PAR \
            || fail "extractorfs died"
    fi

    if notExists "$2/orfs_aa"; then
        # shellcheck disable=SC2086
        "$MMSEQS" translatenucs "$2/orfs" "$2/orfs_aa" $TRANSLATE_PAR \
            || fail "translatenucs died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" indexdb "$2/orfs_aa" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        rm -f "$2/orfs" "$2/orfs.index" "$2/orfs.dbtype"
        rm -f "$2/orfs_aa" "$2/orfs_aa.index" "$2/orfs_aa.dbtype"
        rm -f "$2/createindex.sh"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" indexdb "$INPUT" "$INPUT" $INDEX_PAR \
        || fail "indexdb died"
fi

