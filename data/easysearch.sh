#!/bin/bash -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [[ -d "$1" ]]; then
        echo "$(cd "$1"; pwd)"
    elif [[ -f "$1" ]]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [[ -d $(dirname "$1") ]]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryFASTA> <targetFASTA>|<targetDB> <outFile> <tmp>" && exit 1;
# check paths
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$(abspath "$1")"
TARGET="$(abspath "$2")"
RESULTS="$(abspath "$3")"
TMP_PATH="$(abspath "$4")"

if notExists "${TMP_PATH}/query"; then
   "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/query" \
        || fail "query createdb died"
fi

if notExists "${TARGET}.dbtype"; then
   if notExists "${TMP_PATH}/target"; then
       "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

INTERMEDIATE="${TMP_PATH}/result"
if notExists "${INTERMEDIATE}"; then
   "$MMSEQS" search "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if [ -n "${GREEDY_BEST_HITS}" ]; then
    if notExists "${TMP_PATH}/result_best"; then
       "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
            || fail "Search died"
    fi
    INTERMEDIATE="${TMP_PATH}/result_best"
fi

if notExists "${TMP_PATH}/alis"; then
   "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/alis" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

mv -f "${TMP_PATH}/alis" "${RESULTS}" || fail "Could not move result to ${RESULTS}"
if [[ -f "${TMP_PATH}/alis.index" ]]; then
    mv -f "${TMP_PATH}/alis.index" "${RESULTS}.index" || fail "Could not move result index to ${RESULTS}"
fi


if [[ -n "${REMOVE_TMP}" ]]; then
    echo "Removing temporary files"
    if [[ -n "${GREEDY_BEST_HITS}" ]]; then
        rm -f "${TMP_PATH}/result_best" "${TMP_PATH}/result_best.index"
    fi
    rm -f "${TMP_PATH}/result" "${TMP_PATH}/result.index"
    if [[ ! -n "${LEAVE_INPUT}" ]]; then
        if [[ -f "${TMP_PATH}/target" ]]; then
            rm -f "${TMP_PATH}/target" "${TMP_PATH}/target.index" "${TMP_PATH}/target_h" "${TMP_PATH}/target_h.index" "${TMP_PATH}/target.lookup" "${TMP_PATH}/target.dbtype"
        fi
        rm -f "${TMP_PATH}/query" "${TMP_PATH}/query.index" "${TMP_PATH}/query_h" "${TMP_PATH}/query_h.index" "${TMP_PATH}/query.lookup" "${TMP_PATH}/query.dbtype"
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easysearch.sh"
fi
