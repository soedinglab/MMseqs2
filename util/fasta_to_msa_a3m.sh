#!/bin/sh -e
in_tmp=$(mktemp /tmp/mmseq-msa.XXXXXXXX)
out_tmp=$(mktemp /tmp/mmseq-msa.XXXXXXXX)
a3m_tmp=$(mktemp /tmp/mmseq-a3m.XXXXXXXX)

# Exit, HUP, INT, QUIT, PIPE, TERM
trap "rm -f $in_tmp; rm -f $out_tmp; rm -f $a3m_tmp; exit 1" 0 1 2 3 13 15  

cat /dev/stdin > $in_tmp

# Clustal O does not handle MSAs with only one sequence well
fasta_headers=$(grep -o '^>' $in_tmp | wc -l)
if [ "$fasta_headers" -gt "1" ]; then
	clustalo --force --infile=$in_tmp --outfile=$out_tmp >/dev/null 2>&1
else
	rm -f $out_tmp
	out_tmp=$in_tmp
fi;

reformat.pl fas a3m $out_tmp $a3m_tmp >/dev/null 2>&1 
hhconsensus -maxres 65535 -i $a3m_tmp -o stdout 2>/dev/null

# Cleanup
rm -f $in_tmp
rm -f $out_tmp
rm -f $a3m_tmp
trap 0
exit
