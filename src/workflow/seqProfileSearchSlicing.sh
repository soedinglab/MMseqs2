#!/bin/sh

if [ "$#" -ne 6 ]; then
	echo "Usage : $0 <SeqDB> <ProfileDB> <OutDB> <tmpDir> <eValThs> <k>"
	echo "Where k is the maximum number of hit per sequence. "
	exit -1;
fi

SEQDB=$1
OUTDB=$3
TMP=$4
# Copy of the profile DB that can be reduced as the search processes
ln -s $2 $TMP/profileDB
sort -k1,1 $2.index > $TMP/profileDB.index
PROFILEDB=profileDB

eval=$5
k=$6
MEMORY_FOR_SWAPPING=$(free|grep Mem|awk '{print $4}')#1000000 #10000000
offset=0 #start with the first result

nProfiles=$(wc -l $PROFILEDB.index|cut -f1 -d' ')

while [$nProfiles -g 0]
do
 
    # correct the evalue (searching with profiles vs searching with sequences)
    currentEval=$eval
    nSequences=$(wc -l $SEQDB.index|cut -f1 -d' ')
    currentEval=`echo $currentEval | sed -e 's/[eE]+*/\\*10\\^/'`
    currentEval=`echo "$currentEval*$nProfiles/$nSequences"|bc -l`

    # compute the max number of sequence that are reasonable to swap
    # according to the number of profiles
    let MAX_SEQS=MEMORY_FOR_SWAPPING*1024/nProfiles/90 # 90 bytes/query-result line max.

    rm -f $TMP/aln_4.* $TMP/pref_4* $TMP/searchOut.notSwapped* $TMP/searchOut.current*
    mmseqs search $PROFILEDB $SEQDB $TMP/searchOut.notSwapped $TMP -e $currentEval --max-seqs $MAX_SEQS --profile --offset-result $offset
    mmseqs swapresults $PROFILEDB $SEQDB $TMP/searchOut.notSwapped $TMP/searchOut.current
    # note here : we recover the right evalue, since it is computed according to the target db which is the full profiledb   
    
    if [ -f $TMP/searchOut ]; then
        if [ $(wc -l $TMP/searchOut|cut -f1 -d ' ') -ge 1 ]; then
            echo "Merging with older results..."
            mmseqs mergeffindex $SEQDB $TMP/searchOut.new $TMP/searchOut.current $TMP/searchOut
            mv -f $TMP/searchOut.new $TMP/searchOut
            mv -f $TMP/searchOut.new.index $TMP/searchOut.index
        fi
    else
        echo "First iteration : Creating searchOut..."
        mv -f $TMP/searchOut.current $TMP/searchOut 
        mv -f $TMP/searchOut.current.index $TMP/searchOut.index
    fi
    
    let offset=offset+MAX_SEQS # keep for the prefilter only the next hits
    
    # now remove the profiles that reached their eval threshold
    mmseqs result2stats $PROFILEDB $SEQDB $TMP/searchOut.notSwapped $TMP/searchOut.count  --stat linecount
    mmseqs filterdb $TMP/searchOut.count $TMP/searchOut.toKeep --filter-column 1 --comparison-operator ge --comparison-value $MAX_SEQS

    join <(awk '$3>1{print $1}' $TMP/searchOut.toKeep.index)  $PROFILEDB.index >  $PROFILEDB.index.tmp 
    mv -f $PROFILEDB.index.tmp $PROFILEDB.index # reduce the profile DB
    nProfiles=$(wc -l $PROFILEDB.index|cut -f1 -d' ')
    
done
# Save the results
mv -f $TMP/searchOut $OUTDB
mv -f $TMP/searchOut.index $OUTDB.index

# Clean up
rm -f $TMP/searchOut.toKeep $TMP/searchOut.toKeep.index
rm -f $TMP/searchOut.count $TMP/searchOut.count.index
rm -f $TMP/aln_4.* $TMP/pref_4* $TMP/searchOut.notSwapped* $TMP/searchOut.current*
