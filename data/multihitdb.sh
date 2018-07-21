#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;

if notExists "${OUTDB}"; then
    "${MMSEQS}" createdb "$@" "${OUTDB}" ${CREATEDB_PAR} \
        || fail "createdb failed"
fi

if [ "$("${MMSEQS}" dbtype "${OUTDB}")" = "Nucleotide" ]; then
    mv -f "${OUTDB}" "${OUTDB}_nucl"
    mv -f "${OUTDB}.index" "${OUTDB}_nucl.index"
    mv -f "${OUTDB}_h" "${OUTDB}_nucl_h"
    mv -f "${OUTDB}_h.index" "${OUTDB}_nucl_h.index"
    mv -f "${OUTDB}_member_lookup" "${OUTDB}_nucl_member_lookup"
    mv -f "${OUTDB}_member_lookup.index" "${OUTDB}_nucl_member_lookup.index"
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

    if notExists "${TMP_PATH}/nucl_member_lookup"; then
        "${MMSEQS}" swapdb "${OUTDB}_nucl_member_lookup" "${TMP_PATH}/nucl_set_lookup" ${SWAPDB_PAR} \
            || fail "swapdb failed"
    fi

    if notExists "${OUTDB}_nucl_set.tsv"; then
        "${MMSEQS}" prefixid "${TMP_PATH}/nucl_set_lookup" "${TMP_PATH}/nucl_set.tsv" --tsv ${THREADS_PAR}  \
            || fail "prefixid failed"
    fi

    if notExists "${TMP_PATH}/orf_set_lookup"; then
        "${MMSEQS}" filterdb "${OUTDB}_orf_set_lookup" "${TMP_PATH}/orf_set_lookup" --trim-to-one-column --filter-regex "^.*$" ${THREADS_PAR} \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_set_lookup"; then
        "${MMSEQS}" filterdb "${TMP_PATH}/orf_set_lookup" "${OUTDB}_set_lookup" --mapping-file  "${TMP_PATH}/nucl_set.tsv" ${THREADS_PAR} \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_member_lookup"; then
        "${MMSEQS}" swapdb "${OUTDB}_set_lookup" "${OUTDB}_member_lookup" ${SWAPDB_PAR} \
            || fail "swapdb failed"
    fi
else
    "${MMSEQS}" swapdb "${OUTDB}_member_lookup" "${OUTDB}_set_lookup" ${SWAPDB_PAR} \
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
        rm -f "${TMP_PATH}/nucl_set.tsv"
    fi
    rm -f "${TMP_PATH}/multihitdb.sh"
fi
