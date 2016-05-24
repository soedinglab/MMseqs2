#!/bin/sh
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
[ -z "$MMDIR" ] && echo "Please set the environment variable \$MMDIR to your MMSEQS installation directory." && exit 1;
# check amount of input variables
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

$RUNNER mmseqs diff $OLDDB $NEWDB $TMP/removedSeqs $TMP/mappingSeqs $TMP/newSeqs && checkReturnCode "Diff died"

echo "==================================================="
echo "====== Filter out the new from old sequences ======"
echo "==================================================="
$RUNNER mmseqs order $TMP/newSeqs $NEWDB $TMP/NEWDB.newSeqs && checkReturnCode "Order died"

echo "==================================================="
echo "=== Update the old clustering with the new keys ==="
echo "==================================================="

$RUNNER mmseqs filterdb $OLDCLUST $TMP/OLDCLUST.mapped --mapping-file $TMP/mappingSeqs && checkReturnCode "Filterdb died"



echo "==================================================="
echo "======= Extract representative sequences =========="
echo "==================================================="
$RUNNER mmseqs result2msa   $NEWDB  $NEWDB $TMP/OLDCLUST.mapped $TMP/OLDDB.mapped.repSeq --only-rep-seq  && checkReturnCode "Result2msa died"


echo "==================================================="
echo "======= Search the new sequences against =========="
echo "========= previous (rep seq of) clusters =========="
echo "==================================================="
$RUNNER mmseqs search $TMP/NEWDB.newSeqs $TMP/OLDDB.mapped.repSeq $TMP/newSeqsHits $TMP --max-seqs 1 --only-key-hit  && checkReturnCode "Search died"
$RUNNER mmseqs swapresults $TMP/newSeqsHits $TMP/newSeqsHits.swapped --split 1  && checkReturnCode "Swapresults died"

echo "==================================================="
echo "= Merge the found sequence with previous clustering"
echo "==================================================="
$RUNNER mmseqs mergeffindex $TMP/OLDCLUST.mapped  $TMP/updatedClust $TMP/newSeqsHits.swapped $TMP/OLDCLUST.mapped && checkReturnCode "Mergeffindex died"

echo "==================================================="
echo "=========== Extract unmapped sequences ============"
echo "==================================================="
awk '$3==1{print $1}' $TMP/newSeqsHits.index > $TMP/noHitSeqList
$RUNNER mmseqs order $TMP/noHitSeqList $NEWDB $TMP/toBeClusteredSeparately  && checkReturnCode "Order of no hit seq. died"


echo "==================================================="
echo "===== Cluster separately the alone sequences ======"
echo "==================================================="
$RUNNER mmseqs clusteringworkflow $TMP/toBeClusteredSeparately $TMP/newClusters $TMP && checkReturnCode "Clustering of new seq. died"

#

echo "==================================================="
echo "==== Merge the updated clustering together with ==="
echo "=====         the new clusters               ======"
echo "==================================================="
$RUNNER mmseqs dbconcat $TMP/updatedClust $TMP/newClusters $4 && checkReturnCode "Dbconcat died"



if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
	rm -f "$TMP/newClusters" "$TMP/toBeClusteredSeparately" "$TMP/noHitSeqList" "$TMP/newSeqsHits.index" "$TMP/newSeqsHits" "$TMP/newSeqsHits.swapped"
	rm -f "$TMP/newSeqsHits" "$TMP/NEWDB.newSeqs" "$TMP/OLDCLUST.mapped" "$TMP/mappingSeqs" "$TMP/newSeqs" "$TMP/removedSeqs"
fi



