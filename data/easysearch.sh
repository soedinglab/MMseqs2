#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryFASTA> <targetFASTA>|<targetDB> <outFile> <tmp>" && exit 1;
# check paths
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$2"
RESULTS="$3"
TMP_PATH="$4"

if notExists "${TMP_PATH}/query.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
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

if [ -n "${LINSEARCH}" ] && notExists "${TARGET}.linidx"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createlinindex "${TARGET}" "${TMP_PATH}/index_tmp" ${CREATELININDEX_PAR} \
        || fail "createlinindex died"
fi

INTERMEDIATE="${TMP_PATH}/result"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${SEARCH_MODULE}" "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if [ -n "${GREEDY_BEST_HITS}" ]; then
    if notExists "${TMP_PATH}/result_best.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
            || fail "Search died"
    fi
    INTERMEDIATE="${TMP_PATH}/result_best"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${TMP_PATH}/alis" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

mv -f "${TMP_PATH}/alis" "${RESULTS}" || fail "Could not move result to ${RESULTS}"
if [ -f "${TMP_PATH}/alis.index" ]; then
    mv -f "${TMP_PATH}/alis.index" "${RESULTS}.index" || fail "Could not move result index to ${RESULTS}"
fi
if [ -f "${TMP_PATH}/alis.dbtype" ]; then
    mv -f "${TMP_PATH}/alis.dbtype" "${RESULTS}.dbtype" || fail "Could not move result dbtype to ${RESULTS}"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Removing temporary files"
    if [ -n "${GREEDY_BEST_HITS}" ]; then
        "$MMSEQS" rmdb "${TMP_PATH}/result_best"
    fi
    "$MMSEQS" rmdb "${TMP_PATH}/result"
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            "$MMSEQS" rmdb "${TMP_PATH}/target"
            "$MMSEQS" rmdb "${TMP_PATH}/target_h"
        fi
        "$MMSEQS" rmdb "${TMP_PATH}/query"
        "$MMSEQS" rmdb "${TMP_PATH}/query_h"
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easysearch.sh"
fi
