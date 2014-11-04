#!/bin/sh -e
miss_tmp=$(mktemp /tmp/mmseq-miss.XXXXXXXX)
orig_tmp=$(mktemp /tmp/mmseq-orig.XXXXXXXX)
diff_tmp=$(mktemp /tmp/mmseq-diff.XXXXXXXX)
# Exit, HUP, INT, QUIT, PIPE, TERM

trap "rm -f $miss_tmp; rm -f $orig_tmp; rm -f $diff_tmp; exit 1" 0 1 2 3 13 15  

awk '{ print $1 }' $1 > $miss_tmp
awk '{ print $1 }' $2 > $orig_tmp
diff --new-line-format="" --unchanged-line-format="" $orig_tmp $miss_tmp > $diff_tmp || true
awk 'NR==FNR{a[$0];next;}$1 in a' $diff_tmp $2

# Cleanup
rm -f $miss_tmp
rm -f $orig_tmp
rm -f $diff_tmp
trap 0
exit
