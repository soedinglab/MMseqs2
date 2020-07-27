#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <contigsDb> <taxSeqDB> <taxPerContigDb> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

CONTIGS_DB="$1"
TAX_SEQ_DB="$2"
RESULTS="$3"
TMP_PATH="$4"

if [ ! -e "${TMP_PATH}/orfs_aa.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${CONTIGS_DB}" "${TMP_PATH}/orfs_aa" ${EXTRACT_ORFS_PAR} \
        || fail "extractorfs died"
fi

if [ ! -e "${TMP_PATH}/orfsTax.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomy "${TMP_PATH}/orfs_aa" "${TAX_SEQ_DB}" "${TMP_PATH}/orfsTax" "${TMP_PATH}/tmpTaxonomy" ${TAXONOMY_PAR} \
        || fail "taxonomy died"
fi


if [ ! -e "${TMP_PATH}/orfs_aa_h_swapped.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${TMP_PATH}/orfs_aa_h" "${TMP_PATH}/orfs_aa_h_swapped" \
        || fail "swapdb died"
fi

if [ ! -e "${RESULTS}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" aggregatetaxweights "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_aa_h_swapped" "${TMP_PATH}/orfsTax" "${TMP_PATH}/orfsTax_aln" "${RESULTS}" ${AGGREGATETAX_PAR} \
        || fail "aggregatetaxweights died"
fi



if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_aa"
     # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfsTax"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfsTax_aln"
     # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_aa_h_swapped"

    rm -rf "${TMP_PATH}/tmpTaxonomy"

    rm -f "${TMP_PATH}/taxpercontig.sh"
fi
