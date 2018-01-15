#!/bin/bash
# Clustering workflow script
checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
# 1. Finding exact $k$-mer matches.
notExists "$3/pref"          && $MMSEQS kmermatcher "$INPUT" "$3/pref" ${KMERMATCHER_PAR}                    && checkReturnCode "Kmer matching step died"
# 2. Hamming distance pre-clustering
notExists "$3/pref_rescore1" && $MMSEQS rescorediagonal $INPUT $INPUT "$3/pref" "$3/pref_rescore1" ${HAMMING_PAR} && checkReturnCode "Rescore with hamming distance step died"
notExists "$3/pre_clust"     && $MMSEQS clust  $INPUT "$3/pref_rescore1" "$3/pre_clust" ${CLUSTER_PAR}            && checkReturnCode "Pre-clustering step died"
awk '{ print $1 }' "$3/pre_clust.index" > "$3/order_redundancy"
notExists "$3/input_step_redundancy" && $MMSEQS createsubdb "$3/order_redundancy" $INPUT "$3/input_step_redundancy"     && checkReturnCode "Createsubdb step died"
notExists "$3/pref_filter1"  && $MMSEQS createsubdb "$3/order_redundancy" "$3/pref" "$3/pref_filter1"                    && checkReturnCode "Createsubdb step died"
notExists "$3/pref_filter2"  && $MMSEQS filterdb "$3/pref_filter1" "$3/pref_filter2" --filter-file "$3/order_redundancy" && checkReturnCode "Filterdb step died"
INPUT="$3/input_step_redundancy"
# 3. Ungapped alignment filtering
RESULTDB="$3/pref_filter2"
if [ -n "$FILTER" ]; then
notExists "$3/pref_rescore2" && $MMSEQS rescorediagonal "$INPUT" "$INPUT" "$3/pref_filter2" "$3/pref_rescore2" ${UNGAPPED_ALN_PAR} && checkReturnCode "Ungapped alignment step died"
RESULTDB="$3/pref_rescore2"
fi
# 4. Local gapped sequence alignment.
notExists "$3/aln"           && $RUNNER $MMSEQS align "$INPUT" "$INPUT" $RESULTDB "$3/aln"  ${ALIGNMENT_PAR}            && checkReturnCode "Alignment step died"
# 5. Clustering using greedy set cover.
notExists "$3/clust"         && $MMSEQS clust "$INPUT" "$3/aln" "$3/clust" ${CLUSTER_PAR}  && checkReturnCode "Clustering step died"
notExists "$3/clu"           && $MMSEQS mergeclusters "$1" "$3/clu" "$3/pre_clust" "$3/clust"

# post processing
mv -f "$3/clu" "$2"
mv -f "$3/clu.index" "$2.index"
checkReturnCode "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref" "$3/pref.index"
 rm -f "$3/pref_rescore1" "$3/pref_rescore1.index"
 rm -f "$3/pre_clust" "$3/pre_clust.index"
 rm -f "$3/input_step_redundancy" "$3/input_step_redundancy.index"
 rm -f "$3/pref_filter1" "$3/pref_filter1.index"
 rm -f "$3/pref_filter2" "$3/pref_filter2.index"
 rm -f "$3/aln" "$3/aln.index"
 rm -f "$3/clust" "$3/clust.index"
fi
