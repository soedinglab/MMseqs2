#!/bin/sh -ex
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}


hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand awk
hasCommand join
hasCommand sort


# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <resultDB> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "TMP directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
RESULTS="$3"
TMP_PATH="$4"


# Copy of the profile DB that can be reduced as the search processes
ln -sf $(realpath $2) $TMP_PATH/profileDB
sort -k1,1 $2.index > $TMP_PATH/profileDB.index
ln -sf $(realpath $2).dbtype $TMP_PATH/profileDB.dbtype
PROFILEDB=$TMP_PATH/profileDB
FULLPROFILEDB=$TMP_PATH/fullProfileDB
ln -sf $(realpath $2) $FULLPROFILEDB
ln -sf $(realpath $2.index) $FULLPROFILEDB.index
ln -sf $(realpath $2.dbtype) $FULLPROFILEDB.dbtype



offset=0 #start with the first result

nProfiles=$(awk '{n=n+1}END{print n}' $PROFILEDB.index)
nSequences=$(awk '{n=n+1}END{print n}' $INPUT.index)


STEP=0
while [ $STEP -lt $MAX_STEPS ] && [ $nProfiles -gt 0 ];
do
    
    #MEMORY_FOR_SWAPPING=$(free|grep Mem|awk '{print $4}') #1000000 #10000000
    # compute the max number of sequence that are reasonable to swap
    # according to the number of profiles
    # 90 bytes/query-result line max.
    MAX_SEQS=$(awk -v mem=$AVAIL_MEM -v nProfiles=$nProfiles 'BEGIN{printf("%.0f\n", 1024*(mem/nProfiles)/90);}')
    
    
    # Max result to go in the prefilter list
    SEARCH_LIM=$(awk -v offset=$offset -v maxSeqs=$MAX_SEQS 'BEGIN{print offset+maxSeqs}')

    #echo "Memory for swapping: $MEMORY_FOR_SWAPPING" >> $TMP_PATH/log.txt
    #echo "Current iteration: searching for $MAX_SEQS sequences using $nProfiles profiles with an offset of $offset in the results" >> $TMP_PATH/log.txt

    #rm -f $TMP_PATH/aln_.* $TMP_PATH/pref_* $TMP_PATH/searchOut.notSwapped.$nProfiles* #$TMP_PATH/searchOut.current*

    if notExists "${TMP_PATH}/pref_notSwapped.$offset"; then
        "$MMSEQS" prefilter $PROFILEDB $INPUT $TMP_PATH/pref_notSwapped.$offset --max-seqs $SEARCH_LIM \
            --offset-result $offset ${PREFILTER_PAR} \
             || fail "Pref died"
    fi

    if notExists "${TMP_PATH}/pref.$offset"; then
        "$MMSEQS" swapresults $PROFILEDB $INPUT $TMP_PATH/pref_notSwapped.$offset $TMP_PATH/pref.$offset \
             || fail "Swapping died"
    fi
    # note here : we recover the right evalue, since it is computed according to the target db
    # which is the full profiledb
    if notExists "${TMP_PATH}/aln.$offset"; then
        "$MMSEQS" align $INPUT $FULLPROFILEDB $TMP_PATH/pref.$offset $TMP_PATH/aln.$offset ${ALIGNMENT_PAR} \
             || fail "Alignment died"
    fi

    if [ -f $TMP_PATH/searchOut ]; then
        mmseqs mergedbs $INPUT $TMP_PATH/searchOut.new "${TMP_PATH}/aln.$offset" $TMP_PATH/searchOut
    else
        cp "${TMP_PATH}/aln.$offset" $TMP_PATH/searchOut.new
        cp "${TMP_PATH}/aln.$offset.index" $TMP_PATH/searchOut.new.index
    fi
    "$MMSEQS" filterdb $TMP_PATH/searchOut.new $TMP_PATH/searchOut.new.sorted --sort-entries 1 --filter-column 4 $COMMONS \
        || fail "Filter (sorting) died"
    rm -f $TMP_PATH/searchOut.new{,.index}
    "$MMSEQS" filterdb $TMP_PATH/searchOut.new.sorted $TMP_PATH/searchOut.new.sorted.trunc --extract-lines $MAX_RESULTS_PER_QUERY $COMMONS \
        || fail "Filter (extract lines) died"
    rm -f $TMP_PATH/searchOut.new.sorted{,.index}
    mv -f $TMP_PATH/searchOut.new.sorted.trunc $TMP_PATH/searchOut
    mv -f $TMP_PATH/searchOut.new.sorted.trunc.index $TMP_PATH/searchOut.index


    # now remove the profiles that reached their eval threshold
    notExists $TMP_PATH/searchOut.$offset.count && "$MMSEQS" result2stats $PROFILEDB $INPUT $TMP_PATH/pref_notSwapped.$offset $TMP_PATH/searchOut.$offset.count  --stat linecount $COMMONS
    notExists $TMP_PATH/searchOut.$offset.toKeep && "$MMSEQS" filterdb $TMP_PATH/searchOut.$offset.count $TMP_PATH/searchOut.$offset.toKeep --filter-column 1 --comparison-operator ge --comparison-value $MAX_SEQS $COMMONS

    sort -k1,1 $TMP_PATH/searchOut.$offset.toKeep.index >$TMP_PATH/searchOut.$offset.toKeep.keys
    join $TMP_PATH/searchOut.$offset.toKeep.keys $PROFILEDB.index >$PROFILEDB.index.tmp 
    mv -f $PROFILEDB.index.tmp $PROFILEDB.index # reduce the profile DB
    
    offset=$SEARCH_LIM # keep for the prefilter only the next hits
    nProfiles=$(awk '{n=n+1}END{print n}' $PROFILEDB.index)
    STEP=$((STEP+1))
done

# Save the results
mv -f $TMP_PATH/searchOut $RESULTS
mv -f $TMP_PATH/searchOut.index $RESULTS.index

# Clean up
#rm -f $TMP_PATH/searchOut.toKeep $TMP_PATH/searchOut.toKeep.index
#rm -f $TMP_PATH/searchOut.count $TMP_PATH/searchOut.count.index
#rm -f $TMP_PATH/aln_4.* $TMP_PATH/pref_4* $TMP_PATH/searchOut.notSwapped* $TMP_PATH/searchOut.current*
#rm $TMP_PATH/profileDB*
