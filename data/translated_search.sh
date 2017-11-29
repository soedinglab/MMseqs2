#!/bin/bash
# Translated search workflow
# Helper functions
checkReturnCode () { 
	if [ $? -ne 0 ]; then
	    echo "$1"
	    exit 1
    fi
}
notExists () { 
	[ ! -f "$1" ] 
}

#pre processing
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";


if [ -n "$QUERY_NUCL" ]; then
    INPUT=$1
elif [ -n "$TARGET_NUCL" ]; then
    INPUT=$2
fi

notExists "$3/orfs"    && $MMSEQS extractorfs   "$INPUT" "$4/orfs"     ${ORF_PAR}       && checkReturnCode "extract orfs step died"
notExists "$3/orfs_aa" && $MMSEQS translatenucs "$4/orfs" "$4/orfs_aa" ${TRANSLATE_PAR} && checkReturnCode "translate step died"
# overwrite inpute database
QUERY=$1
QUERY_ORF=$1
TARGET=$2
TARGET_ORF=$2
if [ -n "$QUERY_NUCL" ]; then
    QUERY="$4/orfs_aa"
    QUERY_ORF="$4/orfs"
elif [ -n "$TARGET_NUCL" ]; then
    TARGET="$4/orfs_aa"
    TARGET_ORF="$4/orfs"
fi

mkdir -p "$4/search"
notExists "$4/aln" && ${SEARCH} ${QUERY} ${TARGET} $4/aln $4/search  && checkReturnCode "Search step died"
notExists "$4/aln_offset" && $MMSEQS offsetalignment "$QUERY_ORF" "$TARGET_ORF" "$4/aln"  "$4/aln_offset"  && checkReturnCode "Offset step died"

# post processing
mv -f "$4/aln_offset" "$3"
mv -f "$4/aln_offset.index" "$3.index"
checkReturnCode "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$4/orfs" "$4/orfs.index" "$4/orfs.dbtype"
  rm -f "$4/orfs_aa" "$4/orfs_aa.index" "$4/orfs_aa.dbtype"
fi


