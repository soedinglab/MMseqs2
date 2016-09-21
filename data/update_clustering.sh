#!/bin/bash
checkReturnCode () {
	if [ $? -ne 0 ]; then
	    echo "$1"
	    exit 1
	fi
}

notExists () {
	[ ! -f "$1" ]
}

#pre processing
# check number of input variables
[ "$#" -ne 5 ] && echo "Please provide  <oldDB> <newDB> <oldDB_clustering> <newDB_clustering> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$3" ] &&  echo "$3 not found!" && exit 1;
[   -f "$4" ] &&  echo "$4 exists already!" && exit 1;
[ ! -d "$5" ] &&  echo "tmp directory $5 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

OLDDB=$1 #"../data/DB"
OLDCLUST=$3 #"DBclustered"
NEWDB=$2 #"../data/targetDB"
TMP=$5 #"tmp/"

notExists "$TMP/removedSeqs" && $MMSEQS diffseqdbs $OLDDB $NEWDB $TMP/removedSeqs $TMP/mappingSeqs $TMP/newSeqs && checkReturnCode "Diff died"

echo "==================================================="
echo "====== Filter out the new from old sequences ======"
echo "==================================================="
notExists "$TMP/NEWDB.newSeqs" && $MMSEQS createsubdb $TMP/newSeqs $NEWDB $TMP/NEWDB.newSeqs && checkReturnCode "Order died"

echo "==================================================="
echo "=== Update the old clustering with the new keys ==="
echo "==================================================="
notExists "$TMP/OLDCLUST.mapped" && $MMSEQS filterdb $OLDCLUST $TMP/OLDCLUST.mapped --mapping-file $TMP/mappingSeqs && checkReturnCode "Filterdb died"



echo "==================================================="
echo "======= Extract representative sequences =========="
echo "==================================================="
notExists "$TMP/OLDDB.mapped.repSeq" && $MMSEQS result2msa   $NEWDB  $NEWDB $TMP/OLDCLUST.mapped $TMP/OLDDB.mapped.repSeq --only-rep-seq  && checkReturnCode "Result2msa died"


echo "==================================================="
echo "======= Search the new sequences against =========="
echo "========= previous (rep seq of) clusters =========="
echo "==================================================="
notExists "$TMP/newSeqsHits" && $RUNNER $MMSEQS search $TMP/NEWDB.newSeqs $TMP/OLDDB.mapped.repSeq $TMP/newSeqsHits $TMP --max-seqs 1 $SEARCH_PAR && checkReturnCode "Search died"
notExists "$TMP/newSeqsHits.swapped.all" && $MMSEQS swapresults $TMP/NEWDB.newSeqs $TMP/OLDDB.mapped.repSeq $TMP/newSeqsHits $TMP/newSeqsHits.swapped.all && checkReturnCode "Swapresults died"
notExists "$TMP/newSeqsHits.swapped" && $MMSEQS filterdb $TMP/newSeqsHits.swapped.all $TMP/newSeqsHits.swapped --trim-to-one-column && checkReturnCode "Trimming died"

echo "==================================================="
echo "=  Merge found sequences with previous clustering ="
echo "==================================================="
if [ -f $TMP/newSeqsHits.swapped ]; then
    notExists "$TMP/updatedClust" && $MMSEQS mergedbs $TMP/OLDCLUST.mapped  $TMP/updatedClust $TMP/newSeqsHits.swapped $TMP/OLDCLUST.mapped && checkReturnCode "Mergeffindex died"
else
    notExists "$TMP/updatedClust" && mv $TMP/OLDCLUST.mapped $TMP/updatedClust  && checkReturnCode "Mv Oldclust to update died"
    notExists "$TMP/updatedClust.index" && mv $TMP/OLDCLUST.mapped.index $TMP/updatedClust.index  && checkReturnCode "Mv Oldclust to update died"
fi

echo "==================================================="
echo "=========== Extract unmapped sequences ============"
echo "==================================================="
notExists "$TMP/noHitSeqList" && awk '$3==1{print $1}' $TMP/newSeqsHits.index > $TMP/noHitSeqList
notExists "$TMP/toBeClusteredSeparately" && $MMSEQS createsubdb $TMP/noHitSeqList $NEWDB $TMP/toBeClusteredSeparately  && checkReturnCode "Order of no hit seq. died"


echo "==================================================="
echo "===== Cluster separately the alone sequences ======"
echo "==================================================="
if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files 1/2"
    rm -f $TMP/aln_* $TMP/pref_* $TMP/clu_* $TMP/input_*
else
    mkdir -p "$TMP/search"
    mv $TMP/aln_* $TMP/pref_* $TMP/clu_* $TMP/input_* "$TMP/search"
fi
notExists "$TMP/newClusters" && $MMSEQS cluster $TMP/toBeClusteredSeparately $TMP/newClusters $TMP $CLUST_PAR && checkReturnCode "Clustering of new seq. died"

echo "==================================================="
echo "==== Merge the updated clustering together with ==="
echo "=====         the new clusters               ======"
echo "==================================================="
if [ -f $TMP/newClusters ]; then
    notExists "$4" && $MMSEQS concatdbs $TMP/updatedClust $TMP/newClusters $4 && checkReturnCode "Dbconcat died"
else
    notExists "$4" && mv $TMP/updatedClust $4 && checkReturnCode "Mv died"
    notExists "$4.index" && mv $TMP/updatedClust.index $4.index && checkReturnCode "Mv died"
fi


if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files 2/2"
	rm -f "$TMP/newClusters" "$TMP/toBeClusteredSeparately" "$TMP/toBeClusteredSeparately.index" "$TMP/noHitSeqList" "$TMP/newSeqsHits.index" "$TMP/newSeqsHits" "$TMP/newSeqsHits.swapped" "$TMP/newSeqsHits.swapped.index"
	rm -f "$TMP/newSeqsHits.swapped.all" "$TMP/newSeqsHits.swapped.all.index" "$TMP/NEWDB.newSeqs" "$TMP/NEWDB.newSeqs.index" "$TMP/OLDCLUST.mapped" "$TMP/OLDCLUST.mapped.index"  "$TMP/mappingSeqs" "$TMP/newSeqs" "$TMP/removedSeqs"
	rm -f "$TMP/OLDDB.mapped.repSeq" "$TMP/OLDDB.mapped.repSeq.index" "$TMP/updatedClust" "$TMP/updatedClust.index"
fi



