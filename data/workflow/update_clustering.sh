#!/bin/bash -e

debugWait() {
    if [ -n "${UPDATING_DEBUG}" ]; then
        read -r -n1
    fi
}

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

joinAndReplace() {
    INPUT="$1"
    OUTPUT="$2"
    MAPPING="$3"
    FIELDS="$4"

    LC_ALL=C join -t $'\t' -o "$FIELDS" <(LC_ALL=C sort -T "${TMP_PATH}" -k1,1 "$MAPPING") \
        <(LC_ALL=C sort -T "${TMP_PATH}" -k1,1 "$INPUT") \
        | LC_ALL=C sort -T "${TMP_PATH}" -k1,1 > "$OUTPUT"
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand awk
hasCommand comm
hasCommand join
hasCommand sort

# pre processing
# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <i:oldSequenceDB> <i:newSequenceDB> <i:oldClusteringDB> <o:newMappedSequenceDB> <o:newClusteringDB> <o:tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1;
[   -f "$5.dbtype" ] && echo "$5.dbtype exists already!" && exit 1;
[ ! -d "$6" ] && echo "tmp directory $6 not found!" && exit 1;

OLDDB="$(abspath "$1")"
NEWDB="$(abspath "$2")"
OLDCLUST="$(abspath "$3")"
NEWMAPDB="$(abspath "$4")"
NEWCLUST="$(abspath "$5")"
TMP_PATH="$(abspath "$6")"

if notExists "${TMP_PATH}/removedSeqs"; then
    # shellcheck disable=SC2086
    "$MMSEQS" diffseqdbs "$OLDDB" "$NEWDB" "${TMP_PATH}/removedSeqs" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs" ${DIFF_PAR} \
        || fail "Diff died"
fi

if [ ! -s "${TMP_PATH}/mappingSeqs" ]; then
    cat <<WARN
WARNING: There are no common sequences between $OLDDB and $NEWDB.
If you aim to add the sequences of $NEWDB to your previous clustering $OLDCLUST, you can run:

mmseqs concatdbs \"$OLDDB\" \"$NEWDB\" \"${OLDDB}.withNewSequences\"
mmseqs concatdbs \"${OLDDB}_h\" \"${NEWDB}_h\" \"${OLDDB}.withNewSequences_h\"
mmseqs clusterupdate \"$OLDDB\" \"${OLDDB}.withNewSequences\" \"$OLDCLUST\" \"$NEWCLUST\" \"${TMP_PATH}\"
WARN
    rm -f "${TMP_PATH}/removedSeqs" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs"
    exit 1
fi

if [ -n "${RECOVER_DELETED}" ] && [ -s "${TMP_PATH}/removedSeqs" ]; then
    echo "==================================================="
    echo "============ Recover removed sequences ============"
    echo "==================================================="

    if notExists "${TMP_PATH}/OLDDB.removedMapping"; then
        (
            HIGHESTID="$(sort -T "${TMP_PATH}" -r -n -k1,1 "${NEWDB}.index"| head -n 1 | cut -f1)"
            awk -v highest="$HIGHESTID" \
                'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
                "${TMP_PATH}/removedSeqs" > "${TMP_PATH}/OLDDB.removedMapping"
            cat "${TMP_PATH}/OLDDB.removedMapping" >> "${TMP_PATH}/mappingSeqs"
        ) || fail "Could not create ${TMP_PATH}/OLDDB.removedMapping"
    fi

    if notExists "${TMP_PATH}/NEWDB.withOld.dbtype"; then
        (
            ln -sf "${OLDDB}" "${TMP}/OLDDB.removedDb"
            ln -sf "${OLDDB}_h" "${TMP}/OLDDB.removedDb_h"
            joinAndReplace "${OLDDB}.index" "${TMP_PATH}/OLDDB.removedDb.index" "${TMP_PATH}/OLDDB.removedMapping" "1.2 2.2 2.3"
            joinAndReplace "${OLDDB}_h.index" "${TMP_PATH}/OLDDB.removedDb_h.index" "${TMP_PATH}/OLDDB.removedMapping" "1.2 2.2 2.3"
            joinAndReplace "${OLDDB}.lookup" "${TMP_PATH}/OLDDB.removedDb.lookup" "${TMP_PATH}/OLDDB.removedMapping" "1.2 2.2"
            $MMSEQS concatdbs "$NEWDB" "${TMP_PATH}/OLDDB.removedDb" "${TMP_PATH}/NEWDB.withOld" --preserve-keys
            $MMSEQS concatdbs "${NEWDB}_h" "${TMP_PATH}/OLDDB.removedDb_h" "${TMP_PATH}/NEWDB.withOld_h" --preserve-keys
            cat "${NEWDB}.lookup" "${TMP_PATH}/OLDDB.removedDb.lookup" > "${TMP_PATH}/NEWDB.withOld.lookup"
        ) || fail "Could not create ${TMP_PATH}/NEWDB.withOld"
    fi

    NEWDB="${TMP_PATH}/NEWDB.withOld"

    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files 1/3"
        rm -f "${TMP_PATH}/OLDDB.removedMapping" "${TMP_PATH}/OLDDB.removed"{Db,Db.index,Db_h,Db_h.index,Db.lookup}
    fi
fi

debugWait
echo "==================================================="
echo "=== Update the new sequences with the old keys ===="
echo "==================================================="

if notExists "${TMP_PATH}/newMappingSeqs"; then
    (
        OLDHIGHESTID="$(sort -T "${TMP_PATH}" -r -n -k1,1 "${OLDDB}.index"| head -n 1 | cut -f1)"
        NEWHIGHESTID="$(sort -T "${TMP_PATH}" -r -n -k1,1 "${NEWDB}.index"| head -n 1 | cut -f1)"
        MAXID="$((OLDHIGHESTID>NEWHIGHESTID?OLDHIGHESTID:NEWHIGHESTID))"
        awk -v highest="$MAXID" \
            'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
            "${TMP_PATH}/newSeqs" > "${TMP_PATH}/newSeqs.mapped"
        awk '{ print $2"\t"$1 }' "${TMP_PATH}/mappingSeqs" > "${TMP_PATH}/mappingSeqs.reverse"
        cat "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newMappingSeqs"
        awk '{ print $2 }' "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newSeqs"
    ) || fail "Could not create ${TMP_PATH}/newMappingSeqs"
fi

if notExists "${TMP_PATH}/NEWDB.index"; then
    joinAndReplace "${NEWDB}.index" "${NEWMAPDB}.index" "${TMP_PATH}/newMappingSeqs" "1.2 2.2 2.3" \
        || fail "join died"
fi


if notExists "${TMP_PATH}/NEWDB_h.index"; then
    joinAndReplace "${NEWDB}_h.index" "${NEWMAPDB}_h.index" "${TMP_PATH}/newMappingSeqs" "1.2 2.2 2.3" \
        || fail "join died"
fi

if notExists "${TMP_PATH}/NEWDB.lookup"; then
    joinAndReplace "${NEWDB}.lookup" "${NEWMAPDB}.lookup" "${TMP_PATH}/newMappingSeqs" "1.2 2.2" \
        || fail "join died"
fi

ln -sf "${NEWDB}" "${NEWMAPDB}"
ln -sf "${NEWDB}_h" "${NEWMAPDB}_h"
ln -sf "${NEWDB}.dbtype" "${NEWMAPDB}.dbtype"
ln -sf "${NEWDB}_h.dbtype" "${NEWMAPDB}_h.dbtype"

NEWDB="${NEWMAPDB}"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files 2/3"
    "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld"
    "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld_h"
fi

debugWait
echo "==================================================="
echo "====== Filter out the new from old sequences ======"
echo "==================================================="
if notExists "${TMP_PATH}/NEWDB.newSeqs.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/newSeqs" "$NEWDB" "${TMP_PATH}/NEWDB.newSeqs" ${VERBOSITY} --subdb-mode 1 \
        || fail "Order died"
    ln -sf "${NEWDB}.dbtype" "${TMP_PATH}/NEWDB.newSeqs.dbtype"
fi

debugWait
echo "==================================================="
echo "======= Extract representative sequences =========="
echo "==================================================="
if notExists "${TMP_PATH}/OLDDB.repSeq.dbtype"; then
    "$MMSEQS" result2repseq "$OLDDB" "$OLDCLUST" "${TMP_PATH}/OLDDB.repSeq" \
    || fail "result2repseq died"
    ln -sf "${OLDDB}.dbtype" "${TMP_PATH}/OLDDB.repSeq.dbtype"
fi

debugWait
echo "==================================================="
echo "======== Search the new sequences against ========="
echo "========= previous (rep seq of) clusters =========="
echo "==================================================="
mkdir -p "${TMP_PATH}/search"
if notExists "${TMP_PATH}/newSeqsHits.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${TMP_PATH}/NEWDB.newSeqs" "${TMP_PATH}/OLDDB.repSeq" "${TMP_PATH}/newSeqsHits" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/newSeqsHits.swapped.all"; then
    "$MMSEQS" swapdb "${TMP_PATH}/newSeqsHits" "${TMP_PATH}/newSeqsHits.swapped.all" \
        || fail "Swapresults died"
fi

if [ -s "${TMP_PATH}/newSeqsHits.swapped.all.index" ]; then
    if notExists "${TMP_PATH}/newSeqsHits.swapped"; then
        "$MMSEQS" filterdb "${TMP_PATH}/newSeqsHits.swapped.all" "${TMP_PATH}/newSeqsHits.swapped" --trim-to-one-column \
            || fail "Trimming died"
    fi
fi

debugWait
echo "==================================================="
echo "=  Merge found sequences with previous clustering ="
echo "==================================================="
if [ -f "${TMP_PATH}/newSeqsHits.swapped.dbtype" ]; then
    if notExists "${TMP_PATH}/updatedClust.dbtype"; then
        "$MMSEQS" mergedbs "$OLDCLUST" "${TMP_PATH}/updatedClust" "$OLDCLUST" "${TMP_PATH}/newSeqsHits.swapped" \
            || fail "mergedbs died"
    fi
else
    if notExists "${TMP_PATH}/updatedClust"; then
        ln -sf "$OLDCLUST" "${TMP_PATH}/updatedClust" \
            || fail "Mv Oldclust to update died"
    fi
    if notExists "${TMP_PATH}/updatedClust.index"; then
        ln -sf "$OLDCLUST.index" "${TMP_PATH}/updatedClust.index" \
            || fail "Mv Oldclust to update died"
    fi
    if notExists "${TMP_PATH}/updatedClust.dbtype"; then
        ln -sf "$OLDCLUST.dbtype" "${TMP_PATH}/updatedClust.dbtype" \
            || fail "Mv Oldclust to update died"
    fi
fi

debugWait
echo "==================================================="
echo "=========== Extract unmapped sequences ============"
echo "==================================================="
if notExists "${TMP_PATH}/noHitSeqList.dbtype"; then
    awk '$3==1 {print $1}' "${TMP_PATH}/newSeqsHits.index" > "${TMP_PATH}/noHitSeqList" \
        || fail "awk died"
fi
if notExists "${TMP_PATH}/toBeClusteredSeparately.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/noHitSeqList" "$NEWDB" "${TMP_PATH}/toBeClusteredSeparately" ${VERBOSITY} --subdb-mode 1 \
        || fail "Order of no hit seq. died"
    ln -sf "${NEWDB}.dbtype" "${TMP_PATH}/toBeClusteredSeparately.dbtype"
fi

debugWait
echo "==================================================="
echo "===== Cluster separately the alone sequences ======"
echo "==================================================="

mkdir -p "${TMP_PATH}/cluster"
if notExists "${TMP_PATH}/newClusters.dbtype"; then
    if  [ -s "${TMP_PATH}/toBeClusteredSeparately" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" cluster "${TMP_PATH}/toBeClusteredSeparately" "${TMP_PATH}/newClusters" "${TMP_PATH}/cluster" ${CLUST_PAR} \
            || fail "Clustering of new seq. died"
    fi
fi

debugWait
echo "==================================================="
echo "==== Merge the updated clustering together with ==="
echo "=====         the new clusters               ======"
echo "==================================================="
if [ -f "${TMP_PATH}/newClusters.dbtype" ]; then
    if notExists "$NEWCLUST"; then
        "$MMSEQS" concatdbs "${TMP_PATH}/updatedClust" "${TMP_PATH}/newClusters" "$NEWCLUST" --preserve-keys \
            || fail "Dbconcat died"
    fi
else
    "$MMSEQS" mvdb "${TMP_PATH}/updatedClust" "$NEWCLUST" || fail "Mvdb died"
fi

debugWait
if [ -n "$REMOVE_TMP" ]; then
    rm -f "${TMP_PATH}/newSeqs.mapped" "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newMappingSeqs"
	rm -f  "${TMP_PATH}/noHitSeqList" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs" "${TMP_PATH}/removedSeqs"

    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits.swapped" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/newClusters" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/toBeClusteredSeparately" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/NEWDB.newSeqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/newSeqsHits.swapped.all" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/OLDDB.repSeq" ${VERBOSITY}
    # shellcheck disable=SC2086
	"$MMSEQS" rmdb "${TMP_PATH}/updatedClust" ${VERBOSITY}

	rmdir "${TMP_PATH}/search" "${TMP_PATH}/cluster"

    rm -f "${TMP_PATH}/update_clustering.sh"
fi
