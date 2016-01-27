#!/bin/sh
# Clustering workflow script
checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
#pre processing
[ -z "$MMDIR" ] && echo "Please set the environment variable $MMDIR to your MMSEQS installation directory." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exsists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

STEP=0
QUERYDB="$1"
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
	# call prefilter module
	notExists "$4/pref_$STEP"    && mmseqs prefilter "$QUERYDB" "$2" "$4/pref_$STEP"  $PREFILTER_PAR                         && checkReturnCode "Prefilter died"
	# call alignment module
	notExists "$4/aln_$STEP"     && mmseqs alignment "$QUERYDB" "$2" "$4/pref_$STEP" "$4/aln_$STEP" $ALIGNMENT_PAR           && checkReturnCode "Alignment died"
	# create profiles
	notExists "$4/profile_$STEP" && mmseqs clustertoprofiledb "$4/aln_$STEP" "$1" "$2" "$4/profile_$STEP" $PROFILE_PAR && checkReturnCode "Create profile died"
	QUERYDB="$4/profile_$STEP"
	PREFILTER_PAR=$PREFILTER_PAR" --profile"
	ALIGNMENT_PAR=$ALIGNMENT_PAR" --profile"
	let STEP=STEP+1
done
# post processing
let STEP=STEP-1
cp "$4/aln_$STEP" "$3"
cp "$4/aln_$STEP.index" "$3.index"
checkReturnCode "Could not copy result to $3"
rm -f "$4"/*
