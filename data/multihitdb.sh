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
FASTAS="$2"

if notExists "${OUTDB}/query_nucl"; then
    "${MMSEQS}" createdb ${FASTAS} "${OUTDB}_nucl" ${CREATEDB_PAR} \
        || fail "createdb of query failed"
fi

if notExists "${OUTDB}/query_orfs"; then
    "${MMSEQS}" extractorfs "${OUTDB}_nucl" "${OUTDB}_orfs" ${EXTRACTORFS_PAR}  \
        || fail "extractorfs on Query DB failed"
fi

if notExists "${OUTDB}/query_aa"; then
    "${MMSEQS}" translatenucs "${OUTDB}_orfs" "${OUTDB}" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs on Query"
fi
