#!/bin/bash
# Clustering workflow script
checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
#pre processing
#[ -z "$MMDIR" ] && echo "Please set the environment variable \$MMDIR to your MMSEQS installation directory." && exit 1;
# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

# processing

INPUT="$1"
notExists "$3/aln_redundancy" && $MMSEQS clusthash "$INPUT" "$3/aln_redundancy" ${DETECTREDUNDANCY_PAR}  && checkReturnCode "Fast filter step $STEP died"
notExists "$3/clu_redundancy" && $MMSEQS clust $INPUT "$3/aln_redundancy" "$3/clu_redundancy" ${CLUSTER_PAR} && checkReturnCode "Fast Cluster filter step $STEP died"
awk '{ print $1 }' "$3/clu_redundancy.index" > "$3/order_redundancy"
notExists "$3/input_step_redundancy" && $MMSEQS createsubdb "$3/order_redundancy" $INPUT "$3/input_step_redundancy" && checkReturnCode "MMseqs order step $STEP died"
INPUT="$3/input_step_redundancy"
# call prefilter module
notExists "$3/pref" && $RUNNER $MMSEQS prefilter "$INPUT" "$INPUT" "$3/pref" $PREFILTER_PAR           && checkReturnCode "Prefilter died"
# call alignment module
notExists "$3/aln"  && $RUNNER $MMSEQS align "$INPUT" "$INPUT" "$3/pref" "$3/aln" $ALIGNMENT_PAR  && checkReturnCode "Alignment died"
# call cluster module
notExists "$3/clu_step0"  && $MMSEQS clust "$INPUT" "$3/aln" "$3/clu_step0" $CLUSTER_PAR          && checkReturnCode "Clustering died"
# merg clu_redundancy and clu
notExists "$3/clu" && $MMSEQS mergeclusters "$1" "$3/clu" "$3/clu_redundancy" "$3/clu_step0" && checkReturnCode "Merging of clusters has died"
# post processing
mv -f "$3/clu" "$2"
mv -f "$3/clu.index" "$2.index"
checkReturnCode "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref" "$3/pref.index"
 rm -f "$3/aln" "$3/aln.index"
 rm -f "$3/clu_step0" "$3/clu_step0.index"
 rm -f "$3/clu_redundancy" "$3/clu_redundancy.index"
 rm -f "$3/aln_redundancy" "$3/aln_redundancy.index"
 rm -f "$3/input_step_redundancy" "$3/input_step_redundancy.index"
fi