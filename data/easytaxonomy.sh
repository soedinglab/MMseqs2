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

if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomy "${TMP_PATH}/query" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/taxonomy_tmp" ${TAXONOMY_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/result_lca.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TARGET}"  "${TMP_PATH}/result"  "${TMP_PATH}/result_lca" ${LCA_PAR} \
        || fail "Convert Alignments died"
fi

if notExists "${RESULTS}_lca.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/query" "${TMP_PATH}/result_lca" "${RESULTS}_lca.tsv" ${CREATETSV_PAR} \
        || fail "Convert Alignments died"
fi

if notExists "${RESULTS}_report"; then
    # shellcheck disable=SC2086
    "$MMSEQS" taxonomyreport "${TARGET}" "${TMP_PATH}/result_lca" "${RESULTS}_report" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Removing temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/result"
    if [ -z "${LEAVE_INPUT}" ]; then
        "$MMSEQS" rmdb "${TMP_PATH}/query"
        "$MMSEQS" rmdb "${TMP_PATH}/query_h"
    fi
    rm -rf "${TMP_PATH}/taxonomy_tmp"
    rm -f "${TMP_PATH}/easytaxonomy.sh"
fi
