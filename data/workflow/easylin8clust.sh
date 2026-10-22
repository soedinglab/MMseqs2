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
    $RUNNER "$MMSEQS" createdb "$@" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "query createdb died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lin8clust "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/createtsv.done"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/cluster.tsv" ${CREATETSV_PAR} \
        || fail "Convert Alignments died"
    touch "${TMP_PATH}/createtsv.done"
fi

if notExists "${TMP_PATH}/result2flat.done"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" result2flat "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/rep_seq.fasta" ${RESULT2FLAT_PAR} \
        || fail "result2flat died"
    touch "${TMP_PATH}/result2flat.done"
fi

for FILE in "${TMP_PATH}/cluster.tsv".* "${TMP_PATH}/rep_seq.fasta".*; do
    NAME="${FILE##*/}"
    mv "$FILE" "${RESULTS}_${NAME}"
done
rm -f "${TMP_PATH}/createtsv.done" "${TMP_PATH}/result2flat.done"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/clu_tmp"
    rm -f "${TMP_PATH}/easylin8clust.sh"
fi
