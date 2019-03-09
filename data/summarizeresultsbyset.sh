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
[ "$#" -ne 5 ] && echo "Please provide <quertDB>  <targetDB> <resultDB> <outDB> <tmpDir>" && exit 1;

QUERYDB="$1"
TARGETDB="$2"
RESULTDB="$3"
OUTDB="$4"
TMP_PATH="$5"

if notExists "${TMP_PATH}/result_pos.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" filterdb "${RESULTDB}" "${TMP_PATH}/result_pos" --compute-positions "${TARGETDB}_nucl_orf_aligned_to_contig" ${THREADS_PAR} \
        || fail "filterdb failed"
fi

if notExists "${OUTDB}.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" resultsbyset "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result_pos" "${OUTDB}" "${TMP_PATH}" ${RESULTSBYSET_PAR} \
        || fail "resultsbyset failed"
fi


if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/result_pos"
    rm -f "${TMP_PATH}/summarizeresultsbyset.sh"
fi
