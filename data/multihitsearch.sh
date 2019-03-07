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
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/aggregate"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${OUTPUT}"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${OUTPUT}" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    "$MMSEQS" rmdb "${TMP_PATH}/result"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate"
    rm -f "${TMP_PATH}/multihitsearch.sh"
fi

