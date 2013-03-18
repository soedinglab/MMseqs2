#!/bin/zsh
# generates a HHM for the alignment from stdin, only for alignments with sequence number > 50, and writes the HHM to the stdout

aln=""
while read LINE; do
    aln=$aln$LINE"\n"
done

seqnum=$(echo $aln | grep "^>" | wc -l)
if [ $seqnum -ge 51 ]; then
    echo $aln | /cluster/bioprogs/hh/hhmake -pca 1.4 -pcb 2.5 -pcm 2 -Blosum65 -v 0 -o stdout
fi
