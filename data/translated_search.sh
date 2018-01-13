#!/bin/bash
# Translated search workflow
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERY="$1"
QUERY_ORF="$1"

if [ -n "$QUERY_NUCL" ]; then
    if notExists "$4/q_orfs"; then
        $MMSEQS extractorfs "$1" "$4/q_orfs" ${ORF_PAR} \
            || fail  "extract orfs step died"
    fi
    if notExists "$4/q_orfs_aa"; then
        $MMSEQS translatenucs "$4/q_orfs" "$4/q_orfs_aa" ${TRANSLATE_PAR} \
            || fail  "translate step died"
    fi
    QUERY="$4/q_orfs_aa"
    QUERY_ORF="$4/q_orfs"
fi

TARGET="$2"
TARGET_ORF="$2"
if [ -n "$TARGET_NUCL" ]; then
    if notExists "$4/t_orfs"; then
        $MMSEQS extractorfs "$2" "$4/t_orfs" ${ORF_PAR} \
            || fail  "extract target orfs step died"
    fi
    if notExists "$4/t_orfs_aa"; then
        $MMSEQS translatenucs "$4/t_orfs" "$4/t_orfs_aa" ${TRANSLATE_PAR} \
            || fail  "translate target step died"
    fi
    TARGET="$4/t_orfs_aa"
    TARGET_ORF="$4/t_orfs"
fi


mkdir -p "$4/search"
if notExists "$4/aln"; then
    $SEARCH "${QUERY}" "${TARGET}" "$4/aln" "$4/search" \
        || fail "Search step died"
fi
if notExists "$4/aln_offset"; then
    $MMSEQS offsetalignment "$QUERY_ORF" "$TARGET_ORF" "$4/aln"  "$4/aln_offset" \
        || fail "Offset step died"
fi
(mv -f "$4/aln_offset" "$3" && mv -f "$4/aln_offset.index" "$3.index") \
    || fail "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$4/q_orfs"    "$4/q_orfs.index"    "$4/q_orfs.dbtype"
  rm -f "$4/q_orfs_aa" "$4/q_orfs_aa.index" "$4/q_orfs_aa.dbtype"
  rm -f "$4/t_orfs"    "$4/t_orfs.index"    "$4/t_orfs.dbtype"
  rm -f "$4/t_orfs_aa" "$4/t_orfs_aa.index" "$4/t_orfs_aa.dbtype"
fi


