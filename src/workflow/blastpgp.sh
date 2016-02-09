#!/bin/bash
# Clustering workflow script
checkReturnCode () {
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () {
	[ ! -f "$1" ]
}
#pre processing
[ -z "$MMDIR" ] && echo "Please set the environment variable \$MMDIR to your MMSEQS installation directory." && exit 1;
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
$1=$(pwd)"/"$QUERY_FILE
cd -

cd $(dirname $4)
$4=$(pwd)"/"
cd -

STEP=0
QUERYDB="$1"
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
	# call prefilter module
	$RUNNER mmseqs prefilter "$QUERYDB" "$2" "$4/pref_$STEP"  $PREFILTER_PAR            && checkReturnCode "Prefilter died"

    if [ $STEP -ge 1 ]; then
        # pref -aln
        mmseqs substractresult "$4/pref_$STEP" "$4/aln_0" "$4/pref_next_$STEP" $SUBSTRACT_PAR  && checkReturnCode "Substract died"
        mv -f "$4/pref_next_$STEP" "$4/pref_$STEP"
        mv -f "$4/pref_next_$STEP.index" "$4/pref_$STEP.index"
    fi

	# call alignment module
	$RUNNER mmseqs alignment "$QUERYDB" "$2" "$4/pref_$STEP" "$4/aln_$STEP" $ALIGNMENT_PAR  && checkReturnCode "Alignment died"
# merge all accepted hits to aln_0
    if [ $STEP -gt 0 ]; then
        mmseqs mergeffindex "$QUERYDB" "$4/aln_new" "$4/aln_0" "$4/aln_$STEP"
        mv -f "$4/aln_new" "$4/aln_0"
        mv -f "$4/aln_new.index" "$4/aln_0.index"
    fi
# create profiles with all found hits
    if [ $STEP -ne $((NUM_IT  - 1)) ]; then
        mmseqs result2profile "$QUERYDB" "$2" "$4/aln_0" "$4/profile_$STEP" $PROFILE_PAR && checkReturnCode "Create profile died"
        ln -s $QUERYDB"_h" "$4/profile_$STEP""_h"
        ln -s $QUERYDB"_h.index" "$4/profile_$STEP""_h.index"
    fi
	QUERYDB="$4/profile_$STEP"
    if [ $STEP -eq 0 ]; then
        PREFILTER_PAR=$PREFILTER_PAR" --profile"
        ALIGNMENT_PAR=$ALIGNMENT_PAR" --profile"
        PROFILE_PAR=$PROFILE_PAR" --profile"
    fi
	let STEP=STEP+1
done
# post processing
let STEP=STEP-1
mv -f "$4/aln_0" "$3"
mv -f "$4/aln_0.index" "$3.index"
checkReturnCode "Could not move result to $3"

if [ -n "$KEEP_TEMP" ]; then
 echo "Keeping temporary files"
 exit 0
fi

STEP=0
#while [ $STEP -lt $NUM_IT ]; do
#    rm -f "$4/pref_$STEP" "$4/pref_$STEP.index"
#    rm -f "$4/aln_$STEP" "$4/aln_$STEP.index"
#    rm -f "$4/profile_$STEP" "$4/profile_$STEP.index" "$4/profile_$STEP""_h" "$4/profile_$STEP""_h.index"
#	let STEP=STEP+1
#done