#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
ORIGINAL="$INPUT"

mkdir -p "${TMP_PATH}/linclust"
if notExists "${TMP_PATH}/clu_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" linclust "$INPUT" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/linclust" ${LINCLUST_PAR} \
        || fail "linclust died"
fi

if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
        || faill "createsubdb died"
fi

INPUT="${TMP_PATH}/input_step_redundancy"

if notExists "$TMP_PATH/query_seqs.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractframes "$INPUT" "${TMP_PATH}/query_seqs" ${EXTRACT_FRAMES_PAR}  \
        || fail "Extractframes died"
fi
QUERY="$TMP_PATH/query_seqs"

if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "$QUERY" "$INPUT" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Prefilter step died"
fi


if [ -n "$ALIGNMENT_MODE_NOT_SET" ]; then

    if notExists "${TMP_PATH}/aln_rescore.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" rescorediagonal "$QUERY" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln_ungapped" ${RESCORE_ALN_PAR}  \
             || fail "Alignment step died"
    fi

    if notExists "${TMP_PATH}/pref_subtract.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" subtractdbs "${TMP_PATH}/pref" "${TMP_PATH}/aln_ungapped" "${TMP_PATH}/pref_subtract" ${THREADSANDCOMPRESS_PAR}  \
             || fail "Alignment step died"
    fi

    if notExists "${TMP_PATH}/aln_gapped.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" align "$QUERY" "$INPUT" "${TMP_PATH}/pref_subtract" "${TMP_PATH}/aln_gapped" ${ALIGNMENT_PAR}  \
             || fail "Alignment step died"
    fi

    if notExists "${TMP_PATH}/aln.dbtype"; then
            # shellcheck disable=SC2086
         "$MMSEQS" concatdbs "${TMP_PATH}/aln_ungapped" "${TMP_PATH}/aln_gapped" "${TMP_PATH}/aln" --preserve-keys --take-larger-entry ${THREADSANDCOMPRESS_PAR}\
             || fail "Mergedbs died"
    fi

else
    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERY" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR}  \
             || fail "Alignment step died"
    fi

fi



if notExists "${TMP_PATH}/aln_off.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "${TMP_PATH}/input_step_redundancy" "${QUERY}" \
                              "${TMP_PATH}/input_step_redundancy" "${TMP_PATH}/input_step_redundancy" \
                              "${TMP_PATH}/aln" "${TMP_PATH}/aln_off" ${OFFSETALIGNMENT_PAR} \
        || fail "Offset step died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
     # shellcheck disable=SC2086
     "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_off" "${TMP_PATH}/clu" ${CLUSTER_PAR} \
          || fail "Clustering step died"
fi

# merge clu_redundancy and clu
# shellcheck disable=SC2086
"$MMSEQS" mergeclusters "$ORIGINAL" "$2" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu" ${MERGECLU_PAR} \
        || fail "Merging of clusters has died"

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/query_seqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/perf" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
        # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_off" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu" ${VERBOSITY}

    rm -f "${TMP_PATH}/nucleotide_clustering.sh"
fi

