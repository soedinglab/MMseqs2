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
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

# processing
# call prefilter module
notExists "$3/pref" && mmseqs prefilter "$1" "$1" "$3/pref" $PREFILTER_PAR           && checkReturnCode "Prefilter died"
# call alignment module
notExists "$3/aln"  && mmseqs alignment "$1" "$1" "$3/pref" "$3/aln" $ALIGNMENT_PAR  && checkReturnCode "Alignment died"
# call cluster module
notExists "$3/clu"  && mmseqs cluster   "$1" "$3/aln" "$3/clu" $CLUSTER_PAR          && checkReturnCode "Clustering died"

# post processing
mv -f "$3/clu" "$2"
mv -f "$3/clu.index" "$2.index"
checkReturnCode "Could not move result to $2"
rm -f "$3/pref*"
rm -f "$3/aln*"
rm -f "$3/clu*"
