#!/bin/bash
# Clustering workflow script
# helper functions
checkReturnCode () { 
	if [ $? -ne 0 ]; then
	    echo "$1"
	    exit 1
    fi
}
notExists () { 
	[ ! -f "$1" ] 
}
hasCommand () {
    command -v $1 >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}
#pre processing
#[ -z "$MMDIR" ] && echo "Please set the environment variable \$MMDIR to your MMSEQS installation directory." && exit 1;
# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;

hasCommand awk

export OMP_PROC_BIND=TRUE

INPUT="$1"
notExists "$3/aln_redundancy" && $MMSEQS clusthash "$INPUT" "$3/aln_redundancy" ${DETECTREDUNDANCY_PAR} && checkReturnCode "Fast filter step $STEP died"
notExists "$3/clu_redundancy" && $MMSEQS clust $INPUT "$3/aln_redundancy" "$3/clu_redundancy" ${CLUSTER1_PAR} && checkReturnCode "Fast Cluster filter step $STEP died"
awk '{ print $1 }' "$3/clu_redundancy.index" > "$3/order_redundancy"
notExists "$3/input_step_redundancy" && $MMSEQS createsubdb "$3/order_redundancy" $INPUT "$3/input_step_redundancy" && checkReturnCode "MMseqs order step $STEP died"

INPUT="$3/input_step_redundancy"
STEP=0
while [ $STEP -lt 4 ]; do
    PARAM=PREFILTER${STEP}_PAR
    notExists "$3/pref_step$STEP" \
        && $RUNNER $MMSEQS prefilter "$INPUT" "$INPUT" "$3/pref_step$STEP" ${!PARAM} \
        && checkReturnCode "Prefilter step $STEP died"
    PARAM=ALIGNMENT${STEP}_PAR
    notExists "$3/aln_step$STEP" \
        && $RUNNER $MMSEQS align "$INPUT" "$INPUT" "$3/pref_step$STEP" "$3/aln_step$STEP" ${!PARAM} \
        && checkReturnCode "Alignment step $STEP died"
    PARAM=CLUSTER${STEP}_PAR
    notExists "$3/clu_step$STEP" \
        && $MMSEQS clust "$INPUT" "$3/aln_step$STEP" "$3/clu_step$STEP" ${!PARAM} \
        && checkReturnCode "Clustering step $STEP died"

    NEXTINPUT="$3/input_step$((STEP+1))"
    if [ $STEP -eq 3 ]; then
        notExists "$3/clu" \
            && $MMSEQS mergeclusters "$1" "$3/clu" "$3/clu_redundancy" "$3/clu_step0" "$3/clu_step1" "$3/clu_step2" "$3/clu_step3" \
            && checkReturnCode "Merging of clusters has died"
    else
        notExists "$3/order_step$STEP" \
            && awk '{ print $1 }' "$3/clu_step$STEP.index" > "$3/order_step$STEP" \
            && checkReturnCode "Awk step $STEP died"
        notExists "$NEXTINPUT" \
            && $MMSEQS createsubdb "$3/order_step$STEP" "$INPUT" "$NEXTINPUT" \
            && checkReturnCode "Order step $STEP died"
    fi

	INPUT=$NEXTINPUT
	let STEP=STEP+1
done

# post processing
mv -f "$3/clu" "$2"
mv -f "$3/clu.index" "$2.index"
checkReturnCode "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/clu_redundancy" "$3/clu_redundancy.index"
 rm -f "$3/aln_redundancy" "$3/aln_redundancy.index"
 rm -f "$3/input_step_redundancy" "$3/input_step_redundancy.index"
 STEP=0
 while [ $STEP -lt 4 ]; do
    rm -f "$3/pref_step$STEP" "$3/pref_step$STEP.index"
    rm -f "$3/aln_step$STEP" "$3/aln_step$STEP.index"
    rm -f "$3/clu_step$STEP" "$3/clu_step$STEP.index"
    rm -f "$3/input_step$STEP" "$3/input_step$STEP.index"
    rm -f "$3/order_step$STEP"
	let STEP=STEP+1
 done
fi


