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
    "$MMSEQS" createdb "${QUERY}" "${TMP_PATH}/query" ${CREATEDB_QUERY_PAR} \
        || fail "query createdb died"
    QUERY="${TMP_PATH}/query"
fi

if [ -n "${GPU}" ]; then
    if notExists "${TMP_PATH}/query_pad"; then
        # shellcheck disable=SC2086
        "$MMSEQS" makepaddedseqdb "${TMP_PATH}/query" "${TMP_PATH}/query_pad" ${MAKEPADDEDSEQDB_PAR} \
            || fail "makepaddedseqdb died"
    fi
    QUERY="${TMP_PATH}/query_pad"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"

    if [ -n "${GPU}" ]; then
        if notExists "${TMP_PATH}/target_pad"; then
            # shellcheck disable=SC2086
            "$MMSEQS" makepaddedseqdb "${TMP_PATH}/target" "${TMP_PATH}/target_pad" ${MAKEPADDEDSEQDB_PAR} \
                || fail "makepaddedseqdb died"
        fi
        TARGET="${TMP_PATH}/target_pad"
    fi
fi

if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" rbh "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/rbh_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${RESULTS}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
        fi
        if [ -f "${TMP_PATH}/target_pad" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_pad" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/target_pad_h" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        if [ -f "${TMP_PATH}/query_pad" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_pad" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/query_pad_h" ${VERBOSITY}
        fi
    fi
    rm -rf "${TMP_PATH}/rbh_tmp"
    rm -f "${TMP_PATH}/easyrbh.sh"
fi
