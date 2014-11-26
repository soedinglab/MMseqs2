#!/bin/sh -e
tmp=$(mktemp /tmp/mmseq-cs.XXXXXXXX.a3m)

# Exit, HUP, INT, QUIT, PIPE, TERM
trap "rm -f $tmp; exit 1" 0 1 2 3 13 15  

cat /dev/stdin > $tmp
# parameters x and c were optimized for uniprot a long time ago
cstranslate -x 0.3 -c 4 -i $tmp -o stdout -A $HHLIB/data/cs219.lib 2>/dev/null

rm -f $tmp
trap 0
exit
