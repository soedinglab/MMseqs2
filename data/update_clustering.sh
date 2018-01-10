#!/bin/bash -e

function debugWait() {
    if ((${#UPDATING_DEBUG[@]})); then
        read -n1
    fi
}

function fail() {
    echo "Error: $1"
    exit 1
}

function notExists() {
	[ ! -f "$1" ]
}

function abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d $(dirname "$1") ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}

function joinAndReplace() {
    INPUT="$1"
    OUTPUT="$2"
    MAPPING="$3"
    FIELDS="$4"

    LC_ALL=C join -t $'\t' -o "$FIELDS" <(LC_ALL=C sort -T "$TMP" -k1,1 "$MAPPING") \
        <(LC_ALL=C sort -T "$TMP" -k1,1 "$INPUT") \
        | LC_ALL=C sort -T "$TMP" -k1,1 > "$OUTPUT"
}

function hasCommand () {
    command -v $1 >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand awk
hasCommand comm
hasCommand join
hasCommand sort

# pre processing
# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <i:oldSequenceDB> <i:newSequenceDB> <i:oldClusteringDB> <o:newMappedSequenceDB> <o:newClusteringDB> <o:tmpDir>" && exit 1;

# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$3" ] &&  echo "$3 not found!" && exit 1;
[   -f "$5" ] &&  echo "$5 exists already!" && exit 1;
[ ! -d "$6" ] &&  echo "tmp directory $6 not found!" && exit 1;

OLDDB="$(abspath $1)" #"../data/DB"
NEWDB="$(abspath $2)" #"../data/targetDB"
OLDCLUST="$(abspath $3)" #"DBclustered"
NEWMAPDB="$(abspath $4)"
NEWCLUST="$(abspath $5)"
TMP="$(abspath $6)" #"tmp/"

MMSEQS=${MMSEQS:-"mmseqs"}

hasCommand ${MMSEQS}

if notExists "$TMP/removedSeqs"; then
    $MMSEQS diffseqdbs "$OLDDB" "$NEWDB" "$TMP/removedSeqs" "$TMP/mappingSeqs" "$TMP/newSeqs" ${DIFF_PAR} \
        || fail "Diff died"
fi

if [ ! -s "$TMP/mappingSeqs" ]; then
    cat <<WARN
WARNING: There are no common sequences between $OLDDB and $NEWDB.
If you aim to add the sequences of $NEWDB to your previous clustering $OLDCLUST, you can run:

mmseqs concatdbs \"$OLDDB\" \"$NEWDB\" \"${OLDDB}.withNewSequences\"
mmseqs concatdbs \"${OLDDB}_h\" \"${NEWDB}_h\" \"${OLDDB}.withNewSequences_h\"
mmseqs clusterupdate \"$OLDDB\" \"${OLDDB}.withNewSequences\" \"$OLDCLUST\" \"$NEWCLUST\" \"$TMP\"
WARN
    rm -f "$TMP/removedSeqs"  "$TMP/mappingSeqs" "$TMP/newSeqs"
    exit 1
fi

if [ -n "${RECOVER_DELETED}" ] && [ -s "$TMP/removedSeqs" ]; then
    echo "==================================================="
    echo "============ Recover removed sequences ============"
    echo "==================================================="

    if notExists "$TMP/OLDDB.removedMapping"; then
        (
            HIGHESTID="$(sort -T "$TMP" -r -n -k1,1 "${NEWDB}.index"| head -n 1 | cut -f1)"
            awk -v highest="$HIGHESTID" \
                'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
                "$TMP/removedSeqs" > "$TMP/OLDDB.removedMapping"
            cat "$TMP/OLDDB.removedMapping" >> "$TMP/mappingSeqs"
        ) || fail "Could not create $TMP/OLDDB.removedMapping"
    fi

    if notExists "$TMP/NEWDB.withOld"; then
        (
            ln -sf "${OLDDB}" "${TMP}/OLDDB.removedDb"
            ln -sf "${OLDDB}_h" "${TMP}/OLDDB.removedDb_h"
            joinAndReplace "${OLDDB}.index" "$TMP/OLDDB.removedDb.index" "$TMP/OLDDB.removedMapping" "1.2 2.2 2.3"
            joinAndReplace "${OLDDB}_h.index" "$TMP/OLDDB.removedDb_h.index" "$TMP/OLDDB.removedMapping" "1.2 2.2 2.3"
            joinAndReplace "${OLDDB}.lookup" "$TMP/OLDDB.removedDb.lookup" "$TMP/OLDDB.removedMapping" "1.2 2.2"
            $MMSEQS concatdbs "$NEWDB" "$TMP/OLDDB.removedDb" "$TMP/NEWDB.withOld" --preserve-keys
            $MMSEQS concatdbs "${NEWDB}_h" "$TMP/OLDDB.removedDb_h" "$TMP/NEWDB.withOld_h" --preserve-keys
            cat "${NEWDB}.lookup" "$TMP/OLDDB.removedDb.lookup" > "$TMP/NEWDB.withOld.lookup"
        ) || fail "Could not create $TMP/NEWDB.withOld"
    fi

    NEWDB="$TMP/NEWDB.withOld"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files 1/3"
        rm -f "$TMP/OLDDB.removedMapping" "$TMP/OLDDB.removed"{Db,Db.index,Db_h,Db_h.index,Db.lookup}
    fi
fi

debugWait
echo "==================================================="
echo "=== Update the new sequences with the old keys ===="
echo "==================================================="

if notExists "$TMP/newMappingSeqs"; then
    (
        OLDHIGHESTID="$(sort -T "$TMP" -r -n -k1,1 "${OLDDB}.index"| head -n 1 | cut -f1)"
        NEWHIGHESTID="$(sort -T "$TMP" -r -n -k1,1 "${NEWDB}.index"| head -n 1 | cut -f1)"
        MAXID="$(($OLDHIGHESTID>$NEWHIGHESTID?$OLDHIGHESTID:$NEWHIGHESTID))"
        awk -v highest="$MAXID" \
            'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
            "$TMP/newSeqs" > "$TMP/newSeqs.mapped"
        awk '{ print $2"\t"$1 }' "$TMP/mappingSeqs" > "$TMP/mappingSeqs.reverse"
        cat "$TMP/mappingSeqs.reverse" "$TMP/newSeqs.mapped" > "$TMP/newMappingSeqs"
        awk '{ print $2 }' "$TMP/newSeqs.mapped" > "$TMP/newSeqs"
    ) || fail "Could not create $TMP/newMappingSeqs"
fi

if notExists "$TMP/NEWDB.index"; then
    joinAndReplace "${NEWDB}.index" "${NEWMAPDB}.index" "$TMP/newMappingSeqs" "1.2 2.2 2.3" \
        || fail "join died"
fi


if notExists "$TMP/NEWDB_h.index"; then
    joinAndReplace "${NEWDB}_h.index" "${NEWMAPDB}_h.index" "$TMP/newMappingSeqs" "1.2 2.2 2.3" \
        || fail "join died"
fi

if notExists "$TMP/NEWDB.lookup"; then
    joinAndReplace "${NEWDB}.lookup" "${NEWMAPDB}.lookup" "$TMP/newMappingSeqs" "1.2 2.2" \
        || fail "join died"
fi

ln -sf "${NEWDB}" "${NEWMAPDB}"
ln -sf "${NEWDB}_h" "${NEWMAPDB}_h"
ln -sf "${NEWDB}".dbtype "${NEWMAPDB}".dbtype
NEWDB="${NEWMAPDB}"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files 2/3"
    rm -f "$TMP/NEWDB.withOld"{,.index,.lookup,_h,_h.index}
fi

debugWait
echo "==================================================="
echo "====== Filter out the new from old sequences ======"
echo "==================================================="
if notExists "$TMP/NEWDB.newSeqs"; then
    $MMSEQS createsubdb "$TMP/newSeqs" "$NEWDB" "$TMP/NEWDB.newSeqs" && \
    ln -sf "$NEWDB".dbtype "$TMP/NEWDB.newSeqs".dbtype \
        || fail "Order died"
fi

debugWait
echo "==================================================="
echo "======= Extract representative sequences =========="
echo "==================================================="
if notExists "$TMP/OLDDB.repSeq"; then
    $MMSEQS result2repseq "$OLDDB" "$OLDCLUST" "$TMP/OLDDB.repSeq" \
    && ln -sf "$OLDDB".dbtype "$TMP/OLDDB.repSeq".dbtype \
    || fail "Result2msa died"
fi

debugWait
echo "==================================================="
echo "======== Search the new sequences against ========="
echo "========= previous (rep seq of) clusters =========="
echo "==================================================="
mkdir -p "$TMP/search"
if notExists "$TMP/newSeqsHits"; then
    $MMSEQS search "$TMP/NEWDB.newSeqs" "$TMP/OLDDB.repSeq" "$TMP/newSeqsHits" "$TMP/search" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "$TMP/newSeqsHits.swapped.all"; then
    $MMSEQS swapresults "$TMP/NEWDB.newSeqs" "$TMP/OLDDB.repSeq" "$TMP/newSeqsHits" "$TMP/newSeqsHits.swapped.all" \
        || fail "Swapresults died"
fi

if [ -s "$TMP/newSeqsHits.swapped.all.index" ]; then
    if notExists "$TMP/newSeqsHits.swapped"; then
        $MMSEQS filterdb "$TMP/newSeqsHits.swapped.all" "$TMP/newSeqsHits.swapped" --trim-to-one-column \
            || fail "Trimming died"
    fi
fi

debugWait
echo "==================================================="
echo "=  Merge found sequences with previous clustering ="
echo "==================================================="
if [ -f "$TMP/newSeqsHits.swapped" ]; then
    if notExists "$TMP/updatedClust"; then
        $MMSEQS mergedbs "$OLDCLUST" "$TMP/updatedClust" "$OLDCLUST" "$TMP/newSeqsHits.swapped" \
            || fail "Mergeffindex died"
    fi
else
    if notExists "$TMP/updatedClust"; then
        ln -sf "$OLDCLUST" "$TMP/updatedClust" \
            || fail "Mv Oldclust to update died"
    fi
    if notExists "$TMP/updatedClust.index"; then
        ln -sf "$OLDCLUST.index" "$TMP/updatedClust.index" \
            || fail "Mv Oldclust to update died"
    fi
fi

debugWait
echo "==================================================="
echo "=========== Extract unmapped sequences ============"
echo "==================================================="
if notExists "$TMP/noHitSeqList"; then
    awk '$3==1 {print $1}' "$TMP/newSeqsHits.index" > "$TMP/noHitSeqList" \
        || fail "awk died"
fi
if notExists "$TMP/toBeClusteredSeparately"; then
    $MMSEQS createsubdb "$TMP/noHitSeqList" "$NEWDB" "$TMP/toBeClusteredSeparately" \
    && ln -sf "$NEWDB".dbtype "$TMP/toBeClusteredSeparately".dbtype \
        || fail "Order of no hit seq. died"
fi

debugWait
echo "==================================================="
echo "===== Cluster separately the alone sequences ======"
echo "==================================================="

mkdir -p "$TMP/cluster"
if notExists "$TMP/newClusters"; then
    if  [ -s "$TMP/toBeClusteredSeparately" ]; then
        $MMSEQS cluster "$TMP/toBeClusteredSeparately" "$TMP/newClusters" "$TMP/cluster" ${CLUST_PAR} \
            || fail "Clustering of new seq. died"
    fi
fi

debugWait
echo "==================================================="
echo "==== Merge the updated clustering together with ==="
echo "=====         the new clusters               ======"
echo "==================================================="
if [ -f "$TMP/newClusters" ]; then
    if notExists "$NEWCLUST"; then
        $MMSEQS concatdbs "$TMP/updatedClust" "$TMP/newClusters" "$NEWCLUST" --preserve-keys \
            || fail "Dbconcat died"
    fi
else
    if notExists "$NEWCLUST"; then
        mv "$TMP/updatedClust" "$NEWCLUST" \
            || fail "Mv died"
    fi

    if notExists "${NEWCLUST}.index"; then
        mv "$TMP/updatedClust.index" "${NEWCLUST}.index" \
            || fail "Mv died"
    fi
fi

debugWait
if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files 3/3"
    rm -f "$TMP/newSeqs.mapped" "$TMP/mappingSeqs.reverse" "$TMP/newMappingSeqs"

	rm -f "$TMP/newClusters" "$TMP/newClusters.index" \
	      "$TMP/toBeClusteredSeparately" "$TMP/toBeClusteredSeparately.index" \
	      "$TMP/noHitSeqList" "$TMP/newSeqsHits.index" "$TMP/newSeqsHits" \
	      "$TMP/newSeqsHits.swapped" "$TMP/newSeqsHits.swapped.index"

	rm -f "$TMP/newSeqsHits.swapped.all" "$TMP/newSeqsHits.swapped.all.index" \
	      "$TMP/NEWDB.newSeqs" "$TMP/NEWDB.newSeqs.index" \
	      "$TMP/mappingSeqs" "$TMP/newSeqs" "$TMP/removedSeqs"

	rm -f "$TMP/OLDDB.repSeq" "$TMP/OLDDB.repSeq.index" \
	      "$TMP/updatedClust" "$TMP/updatedClust.index"

	rmdir "$TMP/search" "$TMP/cluster"

    rm -f "$TMP/update_clustering.sh"
fi
