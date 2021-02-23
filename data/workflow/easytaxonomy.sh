#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${TMP_PATH}/query.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$@" "${TMP_PATH}/query" ${CREATEDB_QUERY_PAR} \
        || fail "query createdb died"
fi

if notExists "${TMP_PATH}/result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomy "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/taxonomy_tmp" ${TAXONOMY_PAR} \
        || fail "Search died"
fi

if notExists "${RESULTS}_lca.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/query" "${TMP_PATH}/result" "${RESULTS}_lca.tsv" ${CREATETSV_PAR} \
        || fail "createtsv died"
fi

# shellcheck disable=SC2086
"$MMSEQS" taxonomyreport "${TARGET}" "${TMP_PATH}/result" "${RESULTS}_report" ${TAXONOMYREPORT_PAR} \
        || fail "taxonomyreport died"

#if notExists "${TMP_PATH}/result_aln.dbtype"; then
#    # shellcheck disable=SC2086
#     "$MMSEQS" filterdb "${TMP_PATH}/result" "${TMP_PATH}/result_aln" --extract-lines 1 ${THREADS_COMP_PAR} \
#        || fail "filterdb died"
#fi

if notExists "${TMP_PATH}/result_aln_swapped.dbtype"; then
    # shellcheck disable=SC2086
     "$MMSEQS" swapresults "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/result_aln" "${TMP_PATH}/result_aln_swapped" ${SWAPRESULT_PAR}  \
        || fail "filterdb died"
fi

if notExists "${TMP_PATH}/result_aln_swapped_sum.dbtype"; then
    # shellcheck disable=SC2086
     "$MMSEQS" summarizealis "${TMP_PATH}/result_aln_swapped" "${TMP_PATH}/result_aln_swapped_sum" ${THREADS_COMP_PAR}  \
        || fail "filterdb died"
fi

if notExists "${TMP_PATH}/result_aln_swapped_sum_tax.dbtype"; then
    # shellcheck disable=SC2086
     "$MMSEQS" addtaxonomy "${TARGET}" "${TMP_PATH}/result_aln_swapped_sum" "${TMP_PATH}/result_aln_swapped_sum_tax" ${ADDTAXONOMY_PAR} \
        || fail "filterdb died"
fi

# shellcheck disable=SC2086
"$MMSEQS" createtsv "${TARGET}" "${TMP_PATH}/result_aln_swapped_sum_tax" "${RESULTS}_tophit_report" ${CREATETSV_PAR} \
        || fail "filterdb died"

# shellcheck disable=SC2086
"$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/result_aln" "${RESULTS}_tophit_aln" ${CONVERT_PAR} \
        || fail "convertalis died"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result_aln" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result_aln_swapped" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result_aln_swapped_sum" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result_aln_swapped_sum_tax" ${VERBOSITY}

    rm -rf "${TMP_PATH}/taxonomy_tmp"
    rm -f "${TMP_PATH}/easytaxonomy.sh"
fi
