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

ORFS_DB="${TMP_PATH}/orfs_aa"
if [ ! -e "${ORFS_DB}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${CONTIGS_DB}" "${ORFS_DB}" ${EXTRACT_ORFS_PAR} \
        || fail "extractorfs died"
fi

if [ -n "${ORF_FILTER}" ]; then
    if notExists "${TMP_PATH}/orfs_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" prefilter "${ORFS_DB}" "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_pref" --min-ungapped-score 3 -s 3 -k 6 --diag-score 0 --spaced-kmer-mode 0 --max-seqs 1 ${THREAD_COMP_PAR} \
            || fail "orf filter prefilter died"
    fi

    if notExists "${TMP_PATH}/orfs_aln.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "${ORFS_DB}" "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_pref" "${TMP_PATH}/orfs_aln" --rescore-mode 2 -e 100 ${THREAD_COMP_PAR} \
            || fail "orf filter rescorediagonal died"
    fi

    if notExists "${TMP_PATH}/orfs_aln.list"; then
        awk '$3 > 1 { print $1 }' "${TMP_PATH}/orfs_aln.index" > "${TMP_PATH}/orfs_aln.list"
    fi

    if notExists "${TMP_PATH}/orfs_filter.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/orfs_aln.list" "${ORFS_DB}" "${TMP_PATH}/orfs_filter" ${CREATESUBDB_PAR} \
            || fail "createsubdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/orfs_filter_h" ${VERBOSITY}
    fi

    if notExists "${TMP_PATH}/orfs_filter_h.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/orfs_aln.list" "${ORFS_DB}_h" "${TMP_PATH}/orfs_filter_h" ${CREATESUBDB_PAR} \
            || fail "createsubdb died"
    fi

    ORFS_DB="${TMP_PATH}/orfs_filter"
fi

if [ ! -e "${TMP_PATH}/orfs_tax.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomy "${ORFS_DB}" "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_tax" "${TMP_PATH}/tmp_taxonomy" ${TAXONOMY_PAR} \
        || fail "taxonomy died"
fi

if [ ! -e "${ORFS_DB}_h_swapped.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${ORFS_DB}_h" "${TMP_PATH}/orfs_h_swapped" ${SWAPDB_PAR} \
        || fail "swapdb died"
fi

if [ ! -e "${RESULTS}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" aggregatetaxweights "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_h_swapped" "${TMP_PATH}/orfs_tax" "${TMP_PATH}/orfs_tax_aln" "${RESULTS}" ${AGGREGATETAX_PAR} \
        || fail "aggregatetaxweights died"
fi



if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_aa" ${VERBOSITY}
    if [ -n "${ORF_FILTER}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/orfs_pref" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/orfs_aln" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/orfs_filter" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/orfs_filter_h" ${VERBOSITY}
        rm -f "${TMP_PATH}/orfs_aln.list"
    fi
     # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_tax" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_tax_aln" ${VERBOSITY}
     # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/orfs_h_swapped" ${VERBOSITY}

    rm -rf "${TMP_PATH}/tmp_taxonomy"
    rm -f "${TMP_PATH}/taxpercontig.sh"
fi
