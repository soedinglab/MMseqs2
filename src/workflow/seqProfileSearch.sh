#!/bin/sh

if [ "$#" -ne 6 ]; then
	echo "Usage : $0 <SeqDB> <ProfileDB> <OutDB> <tmpDir> <eValThs> <k>"
	echo "Where k is the maximum number of hit per sequence. "
	exit -1;
fi

SEQDB=$1
PROFILEDB=$2
OUTDB=$3
TMP=$4

eval=$5
k=$6

Kp=10000 # maximum number of hits per profile

sampling=1000 # number of profile we sample to compute the eValue at each iteration
MAX_SEQS=1000000 # maximum number of hit per profile in the sampling step

let n=$Kp*$sampling



ffindexapply=/home/clovis/Software/hh-suite/bin/ffindex_apply_mpi


    # Sample the profiles
    cut -f 1 $PROFILEDB.index|sort -R > $TMP/profiles.lst
    head -n $sampling $TMP/profiles.lst > $TMP/profilesSampled.lst
    tail -n +`expr $sampling + 1` $TMP/profiles.lst > $TMP/profilesRemaining.lst
    mmseqs order $TMP/profilesSampled.lst $PROFILEDB $TMP/profileDBsampled
    mmseqs order $TMP/profilesRemaining.lst $PROFILEDB $TMP/profileDBremaining


while [$eval -le $currentEval]
do
    # Search the samples against the seq DB
    rm $TMP/aln_4.* $TMP/pref_4* -f
    rm $TMP/SamplingSearch{,.index} -f
    mmseqs search $TMP/profileDBsampled $SEQDB $TMP/SamplingSearch $TMP -e $eval --max-seqs $MAX_SEQS --profile
    
    # Evaluate an Evalue that would fit the accepted number of hit per profile
    if [ $(wc -l $TMP/SamplingSearch|cut -f1 -d ' ') -ge $n ]; then # if there was enough hits in the sampling search
	currentEval=`cut -f4 $TMP/SamplingSearch| sort -g | sed -n '$n p'` # grab the eval threshold

	# correct the evalue (sampling profile DB vs remaining profile DB)
	nRemainingProfiles=$(wc -l $TMP/profilesRemaining.lst|cut -f1 -d' ')
	nSamplingProfiles=$(wc -l $TMP/profilesSampled.lst|cut -f1 -d' ')
	currentEval=`echo $currentEval | sed -e 's/[eE]+*/\\*10\\^/'`
	currentEval=`echo $currentEval*$nRemainingProfiles/$nSamplingProfiles|bc -l`
    else
	currentEval=$eval
    fi
    if [ $currentEval -ge $eval]
        currentEval=$eval
    

    # correct the evalue (searching with profiles vs searching with sequences)
    nProfiles=$(wc -l $PROFILEDB.index|cut -f1 -d' ')
    nSequences=$(wc -l $SEQDB.index|cut -f1 -d' ')
    currentEval=`echo $currentEval | sed -e 's/[eE]+*/\\*10\\^/'`
    currentEval=`echo "$currentEval*$nProfiles/$nSequences"|bc -l`

    # Search the remaining profiles with this e-velue threshold
    rm -f $TMP/aln_4* $TMP/pref_4*
    rm -f $TMP/RemainingSearch{,.index}
    mmseqs search $TMP/profileDBremaining $SEQDB $TMP/RemainingSearch $TMP -e $currentEval --max-seqs $Kp --profile
    
    # Trim results to $Kp seq/profile in the sampling search
    if [ -f $TMP/SamplingSearch ] && [ $(wc -l $TMP/SamplingSearch|cut -f1 -d ' ') -ge 1 ]; then
    	mmseqs filterdb $TMP/SamplingSearch $TMP/SamplingSearchTrimmed --extract-lines $Kp
    else
	cp -f $TMP/SamplingSearch $TMP/SamplingSearchTrimmed
	cp -f $TMP/SamplingSearch.index $TMP/SamplingSearchTrimmed.index
    fi

    # Merge the results with the results of the sampling
    mmseqs dbconcat $TMP/SamplingSearchTrimmed $TMP/RemainingSearch $TMP/mergedSearch.notSwapped --preserve-keys
    mmseqs swapresults $PROFILEDB $SEQDB $TMP/mergedSearch.notSwapped $TMP/mergedSearch
         # note here : we recover the right evalue, since it is computed according to the target db which is the full profiledb
    

   # if [ -f $TMP/mergedSearch && $(wc -l $TMP/mergedSearch|cut -f1 -d ' ') -ge 1]; then
   #     mmseqs mergeffindex $SEQDB $TMP/mergedSearch.new $TMP/currentSearch.swapped $TMP/mergedSearch 
   #     mv -f $TMP/mergedSearch.new $TMP/mergedSearch
   #     mv -f $TMP/mergedSearch.new.index $TMP/mergedSearch.index
   # else
   #     mv -f $TMP/currentSearch.swapped $TMP/mergedSearch 
   #     mv -f $TMP/currentSearch.swapped.index $TMP/mergedSearch.index
   # fi

    # Remove the sequences that has reached the maximum number of hits we desire per sequence
    mmseqs result2stats $SEQDB $PROFILEDB $TMP/mergedSearch  $TMP/mergedSearch.count  --stat linecount
    mmseqs filterdb $TMP/mergedSearch.count $TMP/mergedSearch.toKeep --filter-column 1 --comparison-operator ge --comparison-value $k
    #$ffindexapply $TMP/mergedSearch $TMP/mergedSearch.index -d $TMP/mergedSearch.count -i $TMP/mergedSearch.count.index -- ./testNumberOfLines.sh $k
    #mmseqs filterdb $TMP/mergedSearch.count --filter-regex NOP $TMP/mergedSearch.toKeep
    awk '$3>1{print $1}' $TMP/mergedSearch.toKeep.index > $TMP/mergedSearch.toKeep.lst
    rm $TMP/newSeqDb
    mmseqs order $TMP/mergedSearch.toKeep.lst $SEQDB $TMP/newSeqDb

    # Copy the search resuts of those sequences (that reached the max number of hits) 
    # to the final search out file.
    awk '$3==1{print $1}' $TMP/mergedSearch.toKeep.index > $TMP/mergedSearch.reachedMax.lst
    mmseqs order $TMP/mergedSearch.reachedMax.lst $TMP/mergedSearch $TMP/mergedSearch.reachedMax
    if [ -f $TMP/finalSearch ]; then
	mmseqs dbconcat $TMP/mergedSearch.reachedMax $TMP/finalSearch $TMP/finalSearch.tmp --preserve-keys
	mv -f $TMP/finalSearch.tmp  $TMP/finalSearch
	mv -f $TMP/finalSearch.tmp.index  $TMP/finalSearch.index
    else
        mv -f $TMP/mergedSearch.reachedMax $TMP/finalSearch
        mv -f $TMP/mergedSearch.reachedMax.index $TMP/finalSearch.index
    fi



    SEQDB="$TMP/newSeqDb"
    
    # Clean up
    rm -f $TMP/mergedSearch.toKeep $TMP/mergedSearch.toKeep.index
    rm -f $TMP/profileDBsampled $TMP/profileDBsampled.index
    rm -f $TMP/profileDBremaining $TMP/profileDBremaining.index
    rm -f $TMP/RemainingSearch $TMP/RemainingSearch.index
    rm -f $TMP/SamplingSearch $TMP/SamplingSearch.index
    rm -f $TMP/mergedSearch.count $TMP/mergedSearch.count.index
done

mmseqs dbconcat $TMP/mergedSearch $TMP/finalSearch $TMP/finalSearch.tmp --preserve-keys
mv -f $TMP/finalSearch.tmp  $TMP/finalSearch
mv -f $TMP/finalSearch.tmp.index  $TMP/finalSearch.index

# The results should be in $TMP/mergedSearch 
cp $TMP/finalSearch $OUTDB
cp $TMP/finalSearch.index $OUTDB.index

