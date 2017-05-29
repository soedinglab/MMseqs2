#!/bin/bash
# Assembler workflow script
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
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

INPUT="$1"
# 1. Finding exact $k$-mer matches.
notExists "$3/pref"          && $MMSEQS kmermatcher "$INPUT" "$3/pref" ${KMERMATCHER_PAR}                    && checkReturnCode "Kmer matching step died"
# 2. Ungapped alignment
notExists "$3/aln" && $MMSEQS rescorediagonal "$INPUT" "$INPUT" "$3/pref" "$3/aln" ${UNGAPPED_ALN_PAR} && checkReturnCode "Ungapped alignment step died"
# 3. Assemble
notExists "$3/assembly"         && $MMSEQS assembleresults "$INPUT" "$3/aln" "$3/assembly"  && checkReturnCode "Assembly step died"

# post processing
mv -f "$3/assembly" "$2"
mv -f "$3/assembly.index" "$2.index"
checkReturnCode "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref" "$3/pref.index"
 rm -f "$3/aln" "$3/aln.index"
 rm -f "$3/assembly" "$3/assembly.index"
fi