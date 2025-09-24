#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
   [ ! -f "$1" ]
}


if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$@" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "query createdb died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${CLUSTER_MODULE}" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "linclust died"
fi

if notExists "${RESULTS}_protein_cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${RESULTS}_protein_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv protein cluster died"
fi

if notExists "${TMP_PATH}/aln_proteome.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" proteomecluster "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/aln_proteome" "${TMP_PATH}/cluster_count" "${TMP_PATH}/aln_protein" ${PROTEOMECLUSTER_PAR} \
        || fail "proteomecluster died"
fi

if notExists "${RESULTS}_cluster_count.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/cluster_count" "${RESULTS}_cluster_count.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster count report died"
fi

if [ -z "${CASCADED_PROTEOME_CLUSTERING}" ] && notExists "${RESULTS}_proteome_cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_proteome" "${RESULTS}_proteome_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster died"
fi

if [ -n "${PROTEOME_HIDDEN_REPORT}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_proteome_productionReport" "${RESULTS}_proteome_cluster_production.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster died"
fi

if notExists "${RESULTS}_protein_align.tsv" && [ -n "${WRITE_ALIGN_PROTEOME}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_protein" "${RESULTS}_protein_align.tsv" ${THREADS_PAR} \
            || fail "createtsv protein align died"
else
    rm -rf "${TMP_PATH}/aln_protein"*
fi

# cascade 
awk 'NR==FNR { sub(/^\x00/, "", $1); a[$1]; next } !($1 in a)' "${TMP_PATH}/aln_proteome" "${TMP_PATH}/input.source" > "${TMP_PATH}/source_filtered"
SOURCEtoNEXTITERATION="${TMP_PATH}/source_filtered"
STEP=2

SUBDB_LOOKUP_LIST="${TMP_PATH}/input.lookup"

# Corrected while loop condition with proper spacing and quoting
while [ -s "$SOURCEtoNEXTITERATION" ] && [ -n "${CASCADED_PROTEOME_CLUSTERING}" ]; do
    echo "Step $STEP: $(wc -l < "$SOURCEtoNEXTITERATION") sources left"
    # Make "sublookup_STEP" from lines in input.lookup whose 3rd field is in the set from source_filtered
    awk 'NR==FNR {sources[$1]; next} $3 in sources' "$SOURCEtoNEXTITERATION" "${SUBDB_LOOKUP_LIST}" > "${TMP_PATH}/sublookup_${STEP}"

    # Create a smaller DB from sublookup
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/sublookup_${STEP}" "${TMP_PATH}/input" "${TMP_PATH}/input_${STEP}" --subdb-mode 1
    NEXTINPUT="${TMP_PATH}/input_${STEP}"

    # Run linclust on the newly created sub-DB
    echo "Run linclust for iter $STEP" 
    if notExists "${TMP_PATH}/clu_${STEP}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" "${CLUSTER_MODULE}" "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${TMP_PATH}/clu_tmp_${STEP}" ${CLUSTER_PAR} \
            || fail "linclust died"
    fi

    echo "Run createtsv: protein clust result for iter $STEP"
    if notExists "${RESULTS}_protein_cluster_${STEP}.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${RESULTS}_protein_cluster_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv protein cluster died"
    fi

    if [ -n "$REMOVE_TMP" ]; then
        # shellcheck disable=SC2086
        rm -rf "${TMP_PATH}/clu_tmp_${STEP}"
    fi

    echo "Run ProteomeCluster for iter $STEP"
    # Run proteomecluster on the newly created sub-DB
    if notExists "${TMP_PATH}/aln_proteome_${STEP}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" proteomecluster "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${TMP_PATH}/aln_proteome_${STEP}" "${TMP_PATH}/cluster_count_${STEP}" "${TMP_PATH}/aln_protein_${STEP}" ${PROTEOMECLUSTER_PAR} \
            || fail "proteomecluster died"
    fi

    echo "Run createtsv: clustercount report for iter $STEP"
    if notExists "${RESULTS}_cluster_count_${STEP}.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${TMP_PATH}/cluster_count_${STEP}" "${RESULTS}_cluster_count_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv proteome cluster count report died"
    fi

    if [ -n "${PROTEOME_HIDDEN_REPORT}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_proteome_${STEP}_productionReport" "${RESULTS}_proteome_cluster_production_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv proteome cluster died"
    fi

    if notExists "${RESULTS}_protein_align_${STEP}.tsv" && [ -n "${WRITE_ALIGN_PROTEOME}" ]; then
        echo "Run createtsv: protein align result for iter $STEP"
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${NEXTINPUT}" "${TMP_PATH}/aln_protein_${STEP}" "${RESULTS}_protein_align_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv protein align died"
    else
        rm -rf "${TMP_PATH}/aln_protein_${STEP}"*
    fi

    echo "Run concatdbs of aln_proteome for iter $STEP"
    # Concatenate new proteome alignments into the master aln_proteome
    "$MMSEQS" concatdbs "${TMP_PATH}/aln_proteome" "${TMP_PATH}/aln_proteome_${STEP}" "${TMP_PATH}/aln_proteome" --preserve-keys 1

    # Repeat the AWK-based filtering to update source_filtered
    awk 'NR==FNR { sub(/^\x00/, "", $1); a[$1]; next } !($1 in a)' "${TMP_PATH}/aln_proteome" "${TMP_PATH}/input.source" > "${TMP_PATH}/source_filtered"
    SUBDB_LOOKUP_LIST="${TMP_PATH}/sublookup_${STEP}"
    # rm -f "${TMP_PATH}/sublookup_${STEP}"
    STEP=$((STEP + 1))
done

if [ -n "${CASCADED_PROTEOME_CLUSTERING}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_proteome" "${RESULTS}_proteome_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_protein" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_proteome" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/clu_tmp"
    rm -f "${TMP_PATH}/easyproteomecluster.sh"
fi