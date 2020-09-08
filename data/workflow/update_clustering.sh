#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

log() {
    if [ "${VERBOSITY}" = "-v 3" ]; then
        echo "$@"
    fi
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

# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <i:oldSequenceDB> <i:newSequenceDB> <i:oldClusteringDB> <o:newMappedSequenceDB> <o:newClusteringDB> <o:tmpDir>" && exit 1
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1
[   -f "$5.dbtype" ] && echo "$5.dbtype exists already!" && exit 1
[ ! -d "$6" ] && echo "tmp directory $6 not found!" && exit 1

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

if [ -s "${TMP_PATH}/removedSeqs" ]; then
    if [ -n "${RECOVER_DELETED}" ]; then
        log "=== Recover removed sequences"
        if notExists "${TMP_PATH}/OLDDB.removedMapping"; then
            HIGHESTID="$(awk '$1 > max { max = $1 } END { print max }' "${NEWDB}.index")"
            awk -v highest="$HIGHESTID" 'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
                "${TMP_PATH}/removedSeqs" > "${TMP_PATH}/OLDDB.removedMapping"
            cat "${TMP_PATH}/OLDDB.removedMapping" >> "${TMP_PATH}/mappingSeqs"
        fi

        if notExists "${TMP_PATH}/NEWDB.withOld.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" renamedbkeys "${TMP_PATH}/OLDDB.removedMapping" "${OLDDB}" "${TMP_PATH}/OLDDB.removedDb" --subdb-mode 1 ${VERBOSITY} \
                || fail "renamedbkeys died"
            # shellcheck disable=SC2086
            "$MMSEQS" concatdbs "$NEWDB" "${TMP_PATH}/OLDDB.removedDb" "${TMP_PATH}/NEWDB.withOld" --preserve-keys --threads 1 ${VERBOSITY} \
                || fail "concatdbs died"
            # shellcheck disable=SC2086
            "$MMSEQS" concatdbs "${NEWDB}_h" "${TMP_PATH}/OLDDB.removedDb_h" "${TMP_PATH}/NEWDB.withOld_h" --preserve-keys --threads 1 ${VERBOSITY} \
                || fail "concatdbs died"
        fi
        NEWDB="${TMP_PATH}/NEWDB.withOld"

        if [ -n "$REMOVE_TMP" ]; then
            echo "Remove temporary files 1/3"
            rm -f "${TMP_PATH}/OLDDB.removedMapping"
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/OLDDB.removedDb" ${VERBOSITY}
        fi
    else
        if notExists "${TMP_PATH}/OLCLUST.withoutDeletedKeys.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/mappingSeqs" "${OLDCLUST}" "${TMP_PATH}/OLCLUST.withoutDeletedKeys" --subdb-mode 1 ${VERBOSITY} \
                || fail "createsubdb died"
        fi
        if notExists "${TMP_PATH}/OLCLUST.withoutDeleted.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" filterdb "${TMP_PATH}/OLCLUST.withoutDeletedKeys" "${TMP_PATH}/OLCLUST.withoutDeleted" --filter-file "${TMP_PATH}/removedSeqs" --positive-filter ${THREADS_PAR} \
                || fail "filterdb died"
        fi
        OLDCLUST="${TMP_PATH}/OLCLUST.withoutDeleted"
    fi
fi

if notExists "${TMP_PATH}/newMappingSeqs"; then
    log "=== Update new sequences with old keys"
    MAXID="$(awk '$1 > max { max = $1 } END { print max }' "${OLDDB}.index" "${NEWDB}.index")"
    awk -v highest="$MAXID" 'BEGIN { start=highest+1 } { print $1"\t"start; start=start+1; }' \
        "${TMP_PATH}/newSeqs" > "${TMP_PATH}/newSeqs.mapped"
    awk '{ print $2"\t"$1 }' "${TMP_PATH}/mappingSeqs" > "${TMP_PATH}/mappingSeqs.reverse"
    cat "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newMappingSeqs"
    awk '{ print $2 }' "${TMP_PATH}/newSeqs.mapped" > "${TMP_PATH}/newSeqs"
fi

if notExists "${NEWMAPDB}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/newMappingSeqs" "${NEWDB}" "${NEWMAPDB}" ${VERBOSITY} \
        || fail "renamedbkeys died"
fi
NEWDB="${NEWMAPDB}"

if notExists "${TMP_PATH}/NEWDB.newSeqs.dbtype"; then
    log "=== Filter out new from old sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/newSeqs" "$NEWDB" "${TMP_PATH}/NEWDB.newSeqs" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb died"
fi

if notExists "${TMP_PATH}/OLDDB.repSeq.dbtype"; then
    log "=== Extract representative sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "$OLDDB" "$OLDCLUST" "${TMP_PATH}/OLDDB.repSeq" ${RESULT2REPSEQ_PAR} \
        || fail "result2repseq died"
fi

if notExists "${TMP_PATH}/newSeqsHits.dbtype"; then
    log "=== Search new sequences against representatives"
    # shellcheck disable=SC2086
    "$MMSEQS" search "${TMP_PATH}/NEWDB.newSeqs" "${TMP_PATH}/OLDDB.repSeq" "${TMP_PATH}/newSeqsHits" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search died"
fi

if notExists "${TMP_PATH}/newSeqsHits.swapped.all.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${TMP_PATH}/newSeqsHits" "${TMP_PATH}/newSeqsHits.swapped.all" ${THREADS_PAR} \
        || fail "swapdb died"
    awk '$3 > 1 { print 1; exit; }' "${TMP_PATH}/newSeqsHits.swapped.all.index" > "${TMP_PATH}/newSeqsHits.swapped.hasHits"
fi

if [ -s "${TMP_PATH}/newSeqsHits.swapped.hasHits" ] && notExists "${TMP_PATH}/newSeqsHits.swapped.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/newSeqsHits.swapped.all" "${TMP_PATH}/newSeqsHits.swapped" --trim-to-one-column ${THREADS_PAR} \
        || fail "filterdb died"
fi

UPDATEDCLUST="${TMP_PATH}/updatedClust"
if [ -f "${TMP_PATH}/newSeqsHits.swapped.dbtype" ]; then
    if notExists "${TMP_PATH}/updatedClust.dbtype"; then
        log "=== Merge found sequences with previous clustering"
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "$OLDCLUST" "${TMP_PATH}/updatedClust" "$OLDCLUST" "${TMP_PATH}/newSeqsHits.swapped" ${VERBOSITY} \
            || fail "mergedbs died"
    fi
else
    UPDATEDCLUST="$OLDCLUST"
fi

if notExists "${TMP_PATH}/toBeClusteredSeparately.dbtype"; then
    log "=== Extract unmapped sequences"
    awk '$3 == 1 {print $1}' "${TMP_PATH}/newSeqsHits.index" > "${TMP_PATH}/noHitSeqList"
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/noHitSeqList" "$NEWDB" "${TMP_PATH}/toBeClusteredSeparately" ${VERBOSITY} --subdb-mode 1 \
        || fail "createsubdb of not hit seq. died"
fi

if notExists "${TMP_PATH}/newClusters.dbtype" && [ -s "${TMP_PATH}/toBeClusteredSeparately.index" ]; then
    log "=== Cluster separately the singleton sequences"
    # shellcheck disable=SC2086
    "$MMSEQS" cluster "${TMP_PATH}/toBeClusteredSeparately" "${TMP_PATH}/newClusters" "${TMP_PATH}/cluster" ${CLUST_PAR} \
        || fail "cluster of new seq. died"
fi

if [ -f "${TMP_PATH}/newClusters.dbtype" ]; then
    if notExists "$NEWCLUST"; then
        log "=== Merge updated clustering with new clusters"
        # shellcheck disable=SC2086
        "$MMSEQS" concatdbs "${UPDATEDCLUST}" "${TMP_PATH}/newClusters" "$NEWCLUST" --preserve-keys ${THREADS_PAR} \
            || fail "concatdbs died"
    fi
else
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${UPDATEDCLUST}" "$NEWCLUST" ${VERBOSITY}
fi

if [ -n "$REMOVE_TMP" ]; then
    rm -f "${TMP_PATH}/newSeqs.mapped" "${TMP_PATH}/mappingSeqs.reverse" "${TMP_PATH}/newMappingSeqs"
    rm -f "${TMP_PATH}/noHitSeqList" "${TMP_PATH}/mappingSeqs" "${TMP_PATH}/newSeqs" "${TMP_PATH}/removedSeqs"
    rm -f "${TMP_PATH}/newSeqsHits.swapped.hasHits"

    if [ -n "${RECOVER_DELETED}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/NEWDB.withOld_h" ${VERBOSITY}
    else
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/OLCLUST.withoutDeletedKeys" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/OLCLUST.withoutDeleted" ${VERBOSITY}
    fi

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

    rm -rf "${TMP_PATH}/search" "${TMP_PATH}/cluster"
    rm -f "${TMP_PATH}/update_clustering.sh"
fi
