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
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;

QUERY="$1"
TARGET="$2"
OUTPUT="$3"
TMP_PATH="$4"

if notExists "${TMP_PATH}/result"; then
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/result_filter"; then
    "${MMSEQS}" filterdb "${TMP_PATH}/result" "${TMP_PATH}/result_filter" --join-db "${TARGET}_orfs_orf_lookup" --column-to-take 0 ${THREADS_PAR} \
        || fail "filterdb failed"
fi

if notExists "${TMP_PATH}/aggregate"; then
    "${MMSEQS}" aggregate "${QUERY}" "${TARGET}" "${TMP_PATH}/result_filter" "${TMP_PATH}/aggregate" ${AGGREGATE_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/query_orfs_size"; then
    "${MMSEQS}" result2stats "${QUERY}_set_lookup" "${QUERY}_set_lookup" "${QUERY}" "${TMP_PATH}/query_orfs_size" ${RESULT2STATS_PAR} \
        || fail "result2stats failed"
fi

if notExists "${OUTPUT}"; then
    "${MMSEQS}" mergeorfcontigs "${TMP_PATH}/aggregate" "${TMP_PATH}/query_orfs_size" "${OUTPUT}" ${THREADS_PAR} \
        || fail "mergeclusters failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    rm -f "${TMP_PATH}/result" "${TMP_PATH}/result.index"
    rm -f "${TMP_PATH}/result_filter" "${TMP_PATH}/result_filter.index"
    rm -f "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate.index"
    rm -f "${TMP_PATH}/query_orfs_size" "${TMP_PATH}/query_orfs_size.index"
    rm -f "${TMP_PATH}/multihitsearch.sh"
fi

