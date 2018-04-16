#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check amount of input variables
[ "$#" -ne 6 ] && echo "Please provide <queryDB> <targetDB> <targetTaxMap> <ncbiTaxdumpDir> <outDB> <tmpDir>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$3" ] &&  echo "$3 not found!" && exit 1;
if [ -n "${LCA_PAR}" ]; then
    if [ ! -f "$4/names.dmp" ] || [ ! -f "$4/nodes.dmp" ] || [ ! -f "$4/merged.dmp" ] || [ ! -f "$4/delnodes.dmp" ]; then
        echo "Required NCBI Taxonomy files missing!"
        exit 1;
    fi
    NCBI_TAXDUMP="$4"
fi
[   -f "$5" ] &&  echo "$5 exists already!" && exit 1;
[ ! -d "$6" ] &&  echo "tmp directory $6 not found!" && mkdir -p "$6";

INPUT="$1"
TARGET="$2"
TAXON_MAPPING="$3"
RESULTS="$5"
TMP_PATH="$6"

if [ ! -e "${TMP_PATH}/first" ]; then
    mkdir -p "${TMP_PATH}/tmp_hsp1"
    # shellcheck disable=SC2086
    "$MMSEQS" search "${INPUT}" "${TARGET}" "${TMP_PATH}/first" "${TMP_PATH}/tmp_hsp1" ${SEARCH1_PAR} \
        || fail "First search died"

fi
LCA_SOURCE="${TMP_PATH}/first"

# 2bLCA mode
if [ -n "${SEARCH2_PAR}" ]; then
    if [ ! -e "${TMP_PATH}/top1" ]; then
        "$MMSEQS" filterdb "${TMP_PATH}/first" "${TMP_PATH}/top1" --extract-lines 1 \
            || fail "Filterdb died"
    fi

    if [ ! -e "${TMP_PATH}/aligned" ]; then
        "$MMSEQS" extractalignedregion "${INPUT}" "${TARGET}" "${TMP_PATH}/top1" "${TMP_PATH}/aligned" --extract-mode 2 \
            || fail "Extractalignedregion failed"
    fi

    if [ ! -e "${TMP_PATH}/round2" ]; then
        mkdir -p "${TMP_PATH}/tmp_hsp2"
        # shellcheck disable=SC2086
        "$MMSEQS" search "${TMP_PATH}/aligned" "${TARGET}" "${TMP_PATH}/round2" "${TMP_PATH}/tmp_hsp2" ${SEARCH2_PAR} \
            || fail "Second search died"
    fi

    # Concat top hit from first search with all the results from second search
    # We can use filterdb --beats-first to filter out all entries from second search that
    # do not reach the evalue of the top 1 hit
    if [ ! -e "${TMP_PATH}/merged" ]; then
        "$MMSEQS" mergedbs "${TMP_PATH}/top1" "${TMP_PATH}/merged" "${TMP_PATH}/top1" "${TMP_PATH}/round2" \
            || fail "Mergedbs died"
    fi

    if [ ! -e "${TMP_PATH}/2b_ali" ]; then
        "$MMSEQS" filterdb "${TMP_PATH}/merged" "${TMP_PATH}/2b_ali" --beats-first --filter-column 4 --comparison-operator le \
            || fail "First filterdb died"
    fi

    LCA_SOURCE="${TMP_PATH}/2b_ali"
fi

if [ ! -e "${TMP_PATH}/mapping" ]; then
    "$MMSEQS" filterdb "${LCA_SOURCE}" "${TMP_PATH}/mapping" --filter-column 1 --mapping-file "${TAXON_MAPPING}" \
        || fail "Second filterdb died"
fi

if [ ! -e "${TMP_PATH}/taxa" ]; then
    "$MMSEQS" filterdb "${TMP_PATH}/mapping" "${TMP_PATH}/taxa" --filter-column 1 --trim-to-one-column \
        || fail "Third filterdb died"
fi

if [ -n "${LCA_PAR}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TMP_PATH}/taxa" "${NCBI_TAXDUMP}" "${RESULTS}" ${LCA_PAR} \
        || fail "Lca died"
else
    mv -f "${TMP_PATH}/taxa" "${RESULTS}"
    mv -f "${TMP_PATH}/taxa.index" "${RESULTS}.index"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}/tmp_hsp1"
    rm -rf "${TMP_PATH}/tmp_hsp2"
    rm -f "${TMP_PATH}/first" "${TMP_PATH}/first.index"

    if [ -n "${SEARCH2_PAR}" ]; then
        rm -f "${TMP_PATH}/top1" "${TMP_PATH}/top1.index"
        rm -f "${TMP_PATH}/aligned" "${TMP_PATH}/aligned.index" "${TMP_PATH}/round2" "${TMP_PATH}/round2.index"
        rm -f "${TMP_PATH}/merged" "${TMP_PATH}/merged.index" "${TMP_PATH}/2b_ali" "${TMP_PATH}/2b_ali.index"
    fi

    if [ -n "${LCA_PAR}" ]; then
        rm -f "${TMP_PATH}/mapping" "${TMP_PATH}/mapping.index" "${TMP_PATH}/taxa" "${TMP_PATH}/taxa.index"
    else
        rm -f "${TMP_PATH}/mapping" "${TMP_PATH}/mapping.index"
    fi

    rm -f "${TMP_PATH}/taxonomy.sh"
fi
