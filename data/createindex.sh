#!/bin/bash
# Clustering workflow script
checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
# check amount of input variables
[ "$#" -ne 2 ] && echo "Please provide <sequenceDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -d "$2" ] &&  echo "tmp directory $2 not found!" && mkdir -p "$2";

INPUT="$1"
# 1. extract orf
if [ -n "$NUCL" ]; then
    $MMSEQS extractorfs   "$1" "$2/orfs" $ORF_PAR
    $MMSEQS translatenucs "$2/orfs" "$2/orfs_aa" $TRANSLATE_PAR
    $MMSEQS indexdb "$2/orfs_aa" "$1" $INDEX_PAR

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        rm -f "$2/orfs" "$2/orfs.index" "$2/orfs.dbtype"
        rm -f "$2/orfs_aa" "$2/orfs_aa.index" "$2/orfs_aa.dbtype"
    fi
else
    $MMSEQS indexdb "$1" "$1" $INDEX_PAR
fi


