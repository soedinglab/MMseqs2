#!/bin/sh
# Clustering workflow script
# helper functions
checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
hasCommand () {
    command -v $1 >/dev/null 2>&1 || { echo >&2 "Please make sure that $1 is in $PATH."; exit 1; }
}
#pre processing
[ -z "$MMDIR" ] && echo "Please set the environment variable $MMDIR to your MMSEQS installation directory." && exit 1;
# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;
hasCommand ffindex_order
hasCommand awk

# processing
################ clustering step 1 ################
# call prefilter module
notExists "$3/pref_step1" && mmseqs prefilter "$1" "$1" "$3/pref_step1" $PREFILTER1_PAR                                        && checkReturnCode "Prefilter step 1 died"
# call alignment module
notExists "$3/aln_step1"  && mmseqs alignment "$1" "$1" "$3/pref_step1" "$3/aln_step1" $ALIGNMENT1_PAR                         && checkReturnCode "Alignment step 1 died"
# call cluster module
notExists "$3/clu_step1"  && mmseqs cluster   "$1" "$3/aln_step1" "$3/clu_step1" $CLUSTER1_PAR                                 && checkReturnCode "Clustering step 1 died"
# extract representative sequences index
notExists "$3/order_step1" && awk '{ print $1 }' "$3/clu_step1.index" > "$3/order_step1"
notExists "$3/input_step2" && ffindex_order "$3/order_step1" "$1" "$1.index" "$3/input_step2" "$3/input_step2.index"
################ clustering step 2 ################
notExists "$3/pref_step2" && mmseqs prefilter "$3/input_step2" "$3/input_step2" "$3/pref_step2" $PREFILTER2_PAR                && checkReturnCode "Prefilter step 2 died"
# call alignment module
notExists "$3/aln_step2"  && mmseqs alignment "$3/input_step2" "$3/input_step2" "$3/pref_step2" "$3/aln_step2" $ALIGNMENT2_PAR && checkReturnCode "Alignment step 2 died"
# call cluster module
notExists "$3/clu_step2"  && mmseqs cluster   "$3/input_step2" "$3/aln_step2"   "$3/clu_step2"  $CLUSTER2_PAR                  && checkReturnCode "Clustering step 2 died"
# extract representative sequences index
notExists "$3/order_step2" && awk '{ print $1 }' "$3/clu_step2.index" > "$3/order_step2"
notExists "$3/input_step3" && ffindex_order "$3/order_step2" "$1" "$1.index" "$3/input_step3" "$3/input_step3.index"
################ clustering step 3 ################
notExists "$3/pref_step3" && mmseqs prefilter "$3/input_step3" "$3/input_step3" "$3/pref_step3" $PREFILTER3_PAR                && checkReturnCode "Prefilter step 3 died"
# call alignment module
notExists "$3/aln_step3"  && mmseqs alignment "$3/input_step3" "$3/input_step3" "$3/pref_step3" "$3/aln_step3" $ALIGNMENT3_PAR && checkReturnCode "Alignment step 3 died"
# call cluster module
notExists "$3/clu_step3"  && mmseqs cluster   "$3/input_step3" "$3/aln_step3"   "$3/clu_step3"  $CLUSTER3_PAR                  && checkReturnCode "Clustering step 3 died"
# merge cluster results
notExists $3/clu && mmseqs mergecluster "$1" "$3/clu" "$3/clu_step1" "$3/clu_step2" "$3/clu_step3"
# post processing
mv -f "$3/clu" "$2"
mv -f "$3/clu.index" "$2.index"
checkReturnCode "Could not move result to $2"
rm -f "$3/pref_step*"
rm -f "$3/aln_step*"
rm -f "$3/clu_step*"
