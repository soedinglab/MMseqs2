#!/bin/bash -e

# Args = TargetSetOrfs, QuerySetOrfs, SearchResults, lookups, OutputFile

# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d $(dirname "$1") ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}


#pre processing
#[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
MMSEQS="/home/clovisjr/MMseqs2/cmake-build-debug/src/mmseqs"
[ "$#" -ne 7 ] && echo "Please provide <TargetSet> <QuerySet> <SearchResults> <TargetSet_ORF_Lookup> <QuerySet_Set_lookup> <OutputFile> <TempDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$3" ] &&  echo "$3 not found!" && exit 1;
[ ! -f "$4" ] &&  echo "$4 not found!" && exit 1;
[ ! -f "$5" ] &&  echo "$5 not found!" && exit 1;
[   -f "$6" ] &&  echo "$6 exists already!" && exit 1;
[ ! -d "$7" ] &&  echo "tmp directory $7 not found!" && mkdir -p "$4";

TARGETSET="$(abspath $1)"
QUERYSET="$(abspath $2)"
SEARCHRESULTS="$(abspath $3)"
TARGETSETORFLOOKUP="$(abspath $4)"
QUERYSETSETLOOKUP="$(abspath $5)"
OUTPUTFILE="$(abspath $6)"
TMP_PATH="$(abspath $7)"


if notExists "$TMP_PATH/AddedColumn" ; then
    ${MMSEQS} filterdb "$SEARCHRESULTS" "$TMP_PATH/AddedColumn" --join-db ${TARGETSETORFLOOKUP} \
    || fail "filterdb failed"
fi


if notExists "$TMP_PATH/bestHitAggregation" ; then
    ${MMSEQS} aggregate AddedColumn bestHitAggregation --mode bestHit \
    || fail "aggregate --mode bestHit failed"
fi

if notExists "$TMP_PATH/MergedBestHitAggregation" ; then
    ${MMSEQS} mergeclusters bestHitAggregation "$TMP_PATH/MergedBestHitAggregation" ${QUERYSETSETLOOKUP} \
    || fail "mergeclusters failed"
fi

if notExists "$TMP_PATH/nbrORFsTargetSet" ; then
    # Only first arg is used
    ${MMSEQS} result2stats "${TARGETSET}" "${TARGETSET}" "${QUERYSET}" "$TMP_PATH/nbrORFsTargetSet" --stat linecount \
    || fail "result2stats for Target Set failed"
fi

if notExists "$TMP_PATH/nbrORFsQuerySet" ; then
    # Only first arg is used
    ${MMSEQS} result2stats "${QUERYSET}" "${TARGETSET}" "${QUERYSET}" "$TMP_PATH/nbrORFsQuerySet" --stat linecount \
    || fail "result2stats for Query Set failed"
fi

if notExists "${OUTPUTFILE}" ; then
    ${MMSEQS} aggregate "$TMP_PATH/MergedBestHitAggregation" "${OUTPUTFILE}" --mode pval "$TMP_PATH/nbrORFsQuerySet" "$TMP_PATH/nbrORFsTargetSet" --threads 1  \
    || fail "aggregate --mode pval failed"
fi