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
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# TO DO??? add check if $3.dbtype already exists before entire workfolw ???

QUERY="$1"
TARGET="$2"
OUTPUT="$3"
TMP_PATH="$4"

if notExists "${TMP_PATH}/result.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/aggregate.index"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${OUTPUT}.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${OUTPUT}" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rmdir "${TMP_PATH}/search"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate" ${VERBOSITY}
    rm -f "${TMP_PATH}/multihitsearch.sh"
fi

