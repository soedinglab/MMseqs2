#!/bin/sh -e
tmp=$(mktemp /tmp/mmseq-hhmake.XXXXXXXX)

# Exit, HUP, INT, QUIT, PIPE, TERM
trap "rm -f $tmp; exit 1" 0 1 2 3 13 15  

cat /dev/stdin > $tmp

fasta_headers=$(grep -o '^>' $tmp | wc -l)
if [ "$fasta_headers" -ge "50" ]; then
	hhmake -i $tmp -o stdout 2>/dev/null
fi;

rm -f $tmp
trap 0
exit
