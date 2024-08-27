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
    "$MMSEQS" linclust "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" proteomecluster "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/aln_protein" "${TMP_PATH}/aln_proteome" ${PROTEOMECLUSTER_PAR} \
        || fail "Convert Alignments died"
fi

if notExists "${RESULTS}_protein_cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${RESULTS}_protein_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv protein cluster died"
fi

if notExists "${RESULTS}_protein_align.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_protein" "${RESULTS}_protein_align.tsv" ${THREADS_PAR} \
            || fail "createtsv protein align died"
fi

if notExists "${RESULTS}_proteome_cluster.tsv"; then
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
