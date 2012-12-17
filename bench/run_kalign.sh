#!/bin/zsh
# kalign wrapper

input=`cat /dev/stdin`
aln=$(echo $input | kalign 2>/dev/null)
if [ "$aln" = "" ]; then
    aln=$input
fi
echo $aln
