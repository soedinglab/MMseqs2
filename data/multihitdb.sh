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

OUTDB="$1"
TMPDIR="$2"

if [ -n "${NUCL}" ]; then
    mv -f "${OUTDB}" "${OUTDB}_nucl"
    mv -f "${OUTDB}.index" "${OUTDB}_nucl.index"
    mv -f "${OUTDB}_h" "${OUTDB}_nucl_h"
    mv -f "${OUTDB}_h.index" "${OUTDB}_nucl_h.index"
    mv -f "${OUTDB}_set_lookup" "${OUTDB}_nucl_set_lookup"
    mv -f "${OUTDB}_set_lookup.index" "${OUTDB}_nucl_set_lookup.index"
    mv -f "${OUTDB}.lookup" "${OUTDB}_nucl.lookup"
    mv -f "${OUTDB}.dbtype" "${OUTDB}_nucl.dbtype"


    if notExists "${OUTDB}_orf"; then
        "${MMSEQS}" extractorfs "${OUTDB}_nucl" "${OUTDB}_orf" ${EXTRACTORFS_PAR} \
            || fail "extractorfs failed"
    fi

    if notExists "${OUTDB}"; then
        "${MMSEQS}" translatenucs "${OUTDB}_orf" "${OUTDB}" ${TRANSLATENUCS_PAR} \
            || fail "translatenucs failed"
    fi

    if notExists "${TMPDIR}/nucl_element_lookup"; then
        "${MMSEQS}" swapdb "${OUTDB}_nucl_set_lookup" "${TMPDIR}/nucl_element_lookup" \
            || fail "swapdb failed"
    fi

    if notExists "${OUTDB}_nucl_set.tsv"; then
        "${MMSEQS}" prefixid "${TMPDIR}/nucl_element_lookup" "${TMPDIR}/nucl_element.tsv" --tsv \
            || fail "prefixid failed"
    fi

    if notExists "${TMPDIR}/orf_set_lookup"; then
        "${MMSEQS}" filterdb "${OUTDB}_orf_set_lookup" "${TMPDIR}/orf_set_lookup" --trim-to-one-column --filter-regex "^.*$" \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_set_lookup"; then
        "${MMSEQS}" filterdb "${TMPDIR}/orf_set_lookup" "${OUTDB}_set_lookup" --mapping-file  "${TMPDIR}/nucl_element.tsv" \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_member_lookup"; then
        "${MMSEQS}" swapdb "${OUTDB}_set_lookup" "${OUTDB}_member_lookup" \
            || fail "swapdb failed"
    fi
else
    "${MMSEQS}" swapdb "${OUTDB}_set_lookup" "${OUTDB}_member_lookup" ${SWAPDB_PAR} \
            || fail "swapdb failed"
fi

if notExists "${OUTDB}_size"; then
    "${MMSEQS}" result2stats "${OUTDB}" "${OUTDB}" "${OUTDB}_member_lookup" "${OUTDB}_set_size" ${RESULT2STATS_PAR} \
        || fail "result2stats failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    if [ -n "${NUCL}" ]; then
        rm -f "${TMPDIR}/nucl_set.tsv"
    fi
    rm -f "${TMP_PATH}/multihitdb.sh"
fi
