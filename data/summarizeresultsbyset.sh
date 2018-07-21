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
[ "$#" -lt 2 ] && echo "Please provide <outputDB>  <fastFile1> ... <fastFileN> <tmpDir>" && exit 1;

QUERYDB="$1"
TARGETDB="$2"
RESULTDB="$3"
TMP_PATH="$4"

if notExists "${TMP_PATH}/result_pos"; then
    "${MMSEQS}" filterdb "${RESULTDB}" "${TMP_PATH}/result_pos" --compute-positions "${TARGETDB}_orf_set_lookup" \
        || fail "filterdb failed"
fi

if notExists "${OUTDB}"; then
    "${MMSEQS}" resultsbyset "${QUERYDB}" "${TARGETDB}" "${TMP_PATH}/result_pos" "${OUTDB}" "${TMP_PATH}" \
        || fail "resultsbyset failed"
fi


if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir  "${TMP_PATH}/result_pos"  "${TMP_PATH}/result_pos.index"
    rm -f "${TMP_PATH}/multihitdb.sh"
fi
