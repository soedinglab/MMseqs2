#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
   [ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceFASTA> <outFile> <tmp>" && exit 1;
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
RESULTS="$2"
TMP_PATH="$3"

if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "query createdb died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${CLUSTER_MODULE}" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "Search died"
fi


if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

if notExists "${TMP_PATH}/rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_rep" ${THREADS_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu_rep" "${TMP_PATH}/rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

if notExists "${TMP_PATH}/all_seqs.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createseqfiledb "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_seqs" ${THREADS_PAR} \
            || fail "Result2repseq  died"

    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu_seqs" "${TMP_PATH}/all_seqs.fasta" ${VERBOSITY_PAR} \
            || fail "result2flat died"
fi

mv "${TMP_PATH}/all_seqs.fasta"  "${RESULTS}_all_seqs.fasta"
mv "${TMP_PATH}/rep_seq.fasta"  "${RESULTS}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"  "${RESULTS}_cluster.tsv"

if [ -n "${REMOVE_TMP}" ]; then
    echo "Removing temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/input"
    "$MMSEQS" rmdb "${TMP_PATH}/clu_seqs"
    "$MMSEQS" rmdb "${TMP_PATH}/clu_rep"
    "$MMSEQS" rmdb "${TMP_PATH}/clu"
    rm -rf "${TMP_PATH}/clu_tmp"
    rm -f "${TMP_PATH}/easycluster.sh"
fi
