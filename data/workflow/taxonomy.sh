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
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

INPUT="$1"
TARGET="$2"
RESULTS="$3"
TMP_PATH="$4"

if [ ! -e "${TMP_PATH}/first.dbtype" ]; then
    mkdir -p "${TMP_PATH}/tmp_hsp1"
    # shellcheck disable=SC2086
    "$MMSEQS" search "${INPUT}" "${TARGET}" "${TMP_PATH}/first" "${TMP_PATH}/tmp_hsp1" ${SEARCH1_PAR} \
        || fail "First search died"

fi
LCA_SOURCE="${TMP_PATH}/first"

# top hip mode
if [ -n "${TOP_HIT}" ]; then
    if [ ! -e "${TMP_PATH}/top1.dbtype" ]; then
        "$MMSEQS" filterdb "${TMP_PATH}/first" "${TMP_PATH}/top1" --beats-first --filter-column 4 --comparison-operator le \
            || fail "First filterdb died"
    fi
    LCA_SOURCE="${TMP_PATH}/top1"
else
    # 2bLCA mode
    if [ -n "${SEARCH2_PAR}" ]; then
        if [ ! -e "${TMP_PATH}/top1.dbtype" ]; then
            "$MMSEQS" filterdb "${TMP_PATH}/first" "${TMP_PATH}/top1" --extract-lines 1 \
                || fail "Filterdb died"
        fi

        if [ ! -e "${TMP_PATH}/aligned.dbtype" ]; then
            "$MMSEQS" extractalignedregion "${INPUT}" "${TARGET}" "${TMP_PATH}/top1" "${TMP_PATH}/aligned" --extract-mode 2 \
                || fail "Extractalignedregion failed"
        fi

        if [ ! -e "${TMP_PATH}/round2.dbtype" ]; then
                if [ -n "${APPROX_2BLCA}" ]; then
                    if [ ! -e "${TMP_PATH}/first_sub" ]; then
                        # shellcheck disable=SC2086
                        "$MMSEQS" createsubdb  "${TMP_PATH}/aligned" "${TMP_PATH}/first" "${TMP_PATH}/first_sub" ${VERBOSITY} --subdb-mode 1 \
                            || fail "createsubdb"
                    fi
                    # shellcheck disable=SC2086
                    $RUNNER "$MMSEQS" align "${TMP_PATH}/aligned" "${TARGET}" "${TMP_PATH}/first_sub" "${TMP_PATH}/round2" ${SEARCH2_PAR} \
                        || fail "Second search died"
                else
                    mkdir -p "${TMP_PATH}/tmp_hsp2"
                    # shellcheck disable=SC2086
                    "$MMSEQS" search "${TMP_PATH}/aligned" "${TARGET}" "${TMP_PATH}/round2" "${TMP_PATH}/tmp_hsp2" ${SEARCH2_PAR} \
                        || fail "Second search died"
                fi
        fi

        # Concat top hit from first search with all the results from second search
        # We can use filterdb --beats-first to filter out all entries from second search that
        # do not reach the evalue of the top 1 hit
        if [ ! -e "${TMP_PATH}/merged.dbtype" ]; then
            "$MMSEQS" mergedbs "${TMP_PATH}/top1" "${TMP_PATH}/merged" "${TMP_PATH}/top1" "${TMP_PATH}/round2" \
                || fail "Mergedbs died"
        fi

        if [ ! -e "${TMP_PATH}/2b_ali.dbtype" ]; then
            "$MMSEQS" filterdb "${TMP_PATH}/merged" "${TMP_PATH}/2b_ali" --beats-first --filter-column 4 --comparison-operator le \
                || fail "First filterdb died"
        fi

        LCA_SOURCE="${TMP_PATH}/2b_ali"
    fi
fi

if [ "${TAX_OUTPUT}" -eq "0" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TARGET}" "${LCA_SOURCE}" "${RESULTS}" ${LCA_PAR} \
        || fail "Lca died"
elif [ "${TAX_OUTPUT}" -eq "2" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TARGET}" "${LCA_SOURCE}" "${RESULTS}" ${LCA_PAR} \
        || fail "Lca died"
    "$MMSEQS" mvdb "${LCA_SOURCE}" "${RESULTS}_aln"  \
        || fail "mvdb died"
else # return alignment
    "$MMSEQS" mvdb "${LCA_SOURCE}" "${RESULTS}"  \
        || fail "mvdb died"
fi


if [ -n "${REMOVE_TMP}" ]; then
    rm -rf "${TMP_PATH}/tmp_hsp1"
    rm -rf "${TMP_PATH}/tmp_hsp2"

    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/first" ${VERBOSITY}

    if [ -n "${SEARCH2_PAR}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/top1" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aligned" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/round2" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/merged" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/2b_ali" ${VERBOSITY}
        if [ -n "${APPROX_2BLCA}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/first_sub" ${VERBOSITY}
        fi
    fi
    if [ -n "${LCA_PAR}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/mapping" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/taxa" ${VERBOSITY}
    else
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/mapping" ${VERBOSITY}
    fi

    rm -f "${TMP_PATH}/taxonomy.sh"
fi
