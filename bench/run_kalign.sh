#!/bin/zsh
# kalign wrapper: returns the sequence for single-sequence inputs

input=`cat /dev/stdin`
aln=$(echo $input | kalign 2>/dev/null)
if [ "$aln" = "" ]; then
    aln=$input
fi
echo $aln
