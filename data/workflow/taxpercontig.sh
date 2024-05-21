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
        $RUNNER "$MMSEQS" prefilter "${ORFS_DB}" "${TARGETDB_IDX}" "${TMP_PATH}/orfs_pref" ${ORF_FILTER_PREFILTER} \
            || fail "orf filter prefilter died"
    fi

    if notExists "${TMP_PATH}/orfs_aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" rescorediagonal "${ORFS_DB}" "${TARGETDB_IDX}" "${TMP_PATH}/orfs_pref" "${TMP_PATH}/orfs_aln" ${ORF_FILTER_RESCOREDIAGONAL} \
            || fail "orf filter rescorediagonal died"
    fi

    if notExists "${TMP_PATH}/orfs_aln.list"; then
        # shellcheck disable=SC2086
        "$MMSEQS" recoverlongestorf "${ORFS_DB}" "${TMP_PATH}/orfs_aln" "${TMP_PATH}/orfs_aln_recovered.list" ${THREADS_PAR} \
            || fail "recoverlongestorf died"
        awk '$3 > 1 { print $1 }' "${TMP_PATH}/orfs_aln.index" > "${TMP_PATH}/orfs_aln_remain.list"
        cat  "${TMP_PATH}/orfs_aln_recovered.list" "${TMP_PATH}/orfs_aln_remain.list" > "${TMP_PATH}/orfs_aln.list"
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

if [ ! -e "${TMP_PATH}/orfs_h_swapped.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${ORFS_DB}_h" "${TMP_PATH}/orfs_h_swapped" ${SWAPDB_PAR} \
        || fail "swapdb died"
    awk 'BEGIN { printf("%c%c%c%c",5,0,0,0); exit; }' > "${TMP_PATH}/orfs_h_swapped.dbtype"
fi

if [ "${TAX_OUTPUT}" = "0" ] || [ "${TAX_OUTPUT}" = "2" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" aggregatetaxweights "${TAX_SEQ_DB}" "${TMP_PATH}/orfs_h_swapped" "${TMP_PATH}/orfs_tax" "${TMP_PATH}/orfs_tax_aln" "${RESULTS}" ${AGGREGATETAX_PAR} \
        || fail "aggregatetaxweights died"
fi

if [ "${TAX_OUTPUT}" = "1" ] || [ "${TAX_OUTPUT}" = "2" ]; then
    if [ ! -e "${TMP_PATH}/orfs_tax_aln_first.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/orfs_tax_aln" "${TMP_PATH}/orfs_tax_aln_first" --extract-lines 1 ${THREADS_COMP_PAR} \
            || fail "filterdb died"
    fi

    OUT_ALN="${RESULTS}_aln"
    if [ "${TAX_OUTPUT}" = "1" ]; then
        OUT_ALN="${RESULTS}"
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" mergeresultsbyset "${TMP_PATH}/orfs_h_swapped"  "${TMP_PATH}/orfs_tax_aln_first" "${OUT_ALN}" ${THREADS_COMP_PAR} \
        || fail "mvdb died"
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

    if [ "${TAX_OUTPUT}" = "1" ] || [ "${TAX_OUTPUT}" = "2" ]; then
          # shellcheck disable=SC2086
          "$MMSEQS" rmdb "${TMP_PATH}/orfs_tax_aln_first" ${VERBOSITY}
    fi

    rm -rf "${TMP_PATH}/tmp_taxonomy"
    rm -f "${TMP_PATH}/taxpercontig.sh"
fi
