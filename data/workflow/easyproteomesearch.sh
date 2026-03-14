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
    "$MMSEQS" createdb "$@" "${TMP_PATH}/query" ${CREATEDB_PAR} \
        || fail "query createdb died"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

if notExists "${TMP_PATH}/search_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/search_result" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/search_result" "${RESULTS}_search.m8" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

INTERMEDIATE="${TMP_PATH}/parsealn_result"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" parseproteomealignments "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/search_result" "${INTERMEDIATE}" ${PARSEPROTEOMEALN_PAR} \
        || fail "parseproteomealignments died"
fi

if notExists "${RESULTS}.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${RESULTS}_score.tsv" ${THREADS_PAR} \
        || fail "createtsv died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    if [ -z "${LEAVE_INPUT}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY_PAR}
    fi
   
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/search_result" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/parsealn_result" ${VERBOSITY_PAR}

    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easyproteomesearch.sh"
fi