#!/bin/bash -e
# Clustering workflow script
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
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

cd $(dirname $1)
QUERY_FILE=$(basename $1)
ABS_QUERY="$(pwd)/${QUERY_FILE}"
cd -

cd $4
TMP_PATH=$(pwd)
cd -
INPUT=$1
TARGET=$2
# call prefilter module
notExists "$TMP_PATH/pref" && $RUNNER $MMSEQS prefilter   "$INPUT" "$TARGET_DB_PREF" "$TMP_PATH/pref" $PREFILTER_PAR          && checkReturnCode "Prefilter died"
notExists "$TMP_PATH/pref_swaped" && $RUNNER $MMSEQS swapresults "$INPUT" "$TARGET_DB_PREF" "$TMP_PATH/pref" "$TMP_PATH/pref_swaped" && checkReturnCode "Swapresults pref died"
# call alignment module
notExists "$TMP_PATH/aln_swaped"  && $RUNNER $MMSEQS align       "$TARGET_DB_PREF" "$INPUT" "$TMP_PATH/pref_swaped" "$TMP_PATH/aln_swaped" $ALIGNMENT_PAR  && checkReturnCode "Alignment died"
notExists "$TMP_PATH/aln"         && $RUNNER $MMSEQS swapresults "$TARGET_DB_PREF" "$INPUT" "$TMP_PATH/aln_swaped"  "$TMP_PATH/aln" && checkReturnCode "Swapresults aln died"

# post processing
mv -f "$TMP_PATH/aln" "$3"
mv -f "$TMP_PATH/aln.index" "$3.index"
checkReturnCode "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -f "$TMP_PATH/pref" "$TMP_PATH/pref.index"
    rm -f "$TMP_PATH/pref_swaped" "$TMP_PATH/pref_swaped.index"
    rm -f "$TMP_PATH/aln_swaped" "$TMP_PATH/aln_swaped.index"
fi


