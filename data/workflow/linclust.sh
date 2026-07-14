#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

if [ "$LINCLUST_MODULE" = "linclust2" ]; then
    # 0. clusthash
    if [ -n "$CLUSTHASH" ]; then
        if notExists "${TMP_PATH}/input_clusthash.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" clusthash "$INPUT" "${TMP_PATH}/input_clusthash" ${CLUSTHASH_PAR} \
                || fail "clusthash died"
        fi

        if notExists "${TMP_PATH}/input_clusthash_clust.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" clust "$INPUT" "${TMP_PATH}/input_clusthash" "${TMP_PATH}/input_clusthash_clust" ${CLUSTHASH_CLUST_PAR} \
                || fail "clusthash-based clust died"
        fi

        awk '{print $1}' "${TMP_PATH}/input_clusthash_clust.index" \
            > "${TMP_PATH}/order_clusthash_redundancy"

        if notExists "${TMP_PATH}/input_clusthash_redundancy.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/order_clusthash_redundancy" "$SOURCE" \
                "${TMP_PATH}/input_clusthash_redundancy" ${VERBOSITY} --subdb-mode 1 \
                || fail "createsubdb (clusthash representatives) died"
        fi
        INPUT="${TMP_PATH}/input_clusthash_redundancy"
    fi

    # 1. k-mer matching
    if notExists "${TMP_PATH}/pref.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
            || fail "kmermatcher died"
    fi
    RESULTDB="${TMP_PATH}/pref"

    # 2. Alignment + clustering (ungapped / gapped banded-block aligner)
    if notExists "${TMP_PATH}/clu.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" align2clust "$INPUT" "$RESULTDB" "${TMP_PATH}/clu" ${ALIGN2CLUST_PAR} \
            || fail "align2clust died"
    fi

    # 2a. Merge clusthash pre-clusters with alignment clusters
    if [ -n "$CLUSTHASH" ]; then
        if notExists "${TMP_PATH}/clu_merged.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu_merged" \
                "${TMP_PATH}/input_clusthash_clust" "${TMP_PATH}/clu" $MERGECLU_PAR \
                || fail "mergeclusters (clusthash + alignment) died"
        fi
        CLUDB="${TMP_PATH}/clu_merged"
    else
        CLUDB="${TMP_PATH}/clu"
    fi

    # 3. Refinement pass: re-cluster representative sequences
    if notExists "${TMP_PATH}/input_rep.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "$CLUDB" "$INPUT" "${TMP_PATH}/input_rep" ${VERBOSITY} --subdb-mode 1 \
            || fail "createsubdb (representatives) died"
    fi

    if notExists "${TMP_PATH}/pref_rep.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcher "${TMP_PATH}/input_rep" "${TMP_PATH}/pref_rep" ${KMERMATCHER_PAR2} \
            || fail "kmermatcher (representatives) died"
    fi

    if notExists "${TMP_PATH}/clu_rep.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" align2clust "${TMP_PATH}/input_rep" "${TMP_PATH}/pref_rep" "${TMP_PATH}/clu_rep" \
            ${ALIGN2CLUST_PAR} \
            --filter-cludb-file "$CLUDB" \
            --filter-seqdb-file "$SOURCE" \
            || fail "align2clust (representatives) died"
    fi

    if notExists "$2.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergeclusters "$SOURCE" "$2" \
            "$CLUDB" "${TMP_PATH}/clu_rep" $MERGECLU_PAR \
            || fail "mergeclusters died"
    fi

    # Expose alignment results (only produced when --include-align-files is set).
    # The two align2clust passes each emit their own alignments; union them keyed by the
    # final representatives ($2) and keep only lines whose target is a member of that
    # final cluster (--merge-filter-target), so the result has exactly one entry per
    # cluster containing exactly that cluster's rep->member alignments.
    if [ -f "${TMP_PATH}/clu_aln.dbtype" ]; then
        if [ -f "${TMP_PATH}/clu_rep_aln.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${2}" "${2}_aln" \
                "${TMP_PATH}/clu_aln" "${TMP_PATH}/clu_rep_aln" \
                --merge-filter-target 1 ${VERBOSITY} \
                || fail "mergedbs clu_aln died"
        else
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "${TMP_PATH}/clu_aln" "${2}_aln" ${VERBOSITY} \
                || fail "mvdb clu_aln died"
        fi
    fi

    # Optionally replace representatives by the most profile-consistent observed member,
    # reusing the alignments in ${2}_aln (no profile-vs-member realignment).
    if [ -n "$SWITCH_CONSENSUS_REP" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" pickconsensusrepfast "$1" "$2" "${TMP_PATH}/clu_switched" "${TMP_PATH}/switch_tmp" ${PICKREP_PAR} \
            || fail "pickconsensusrepfast (switch representatives) died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "$2" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" mvdb "${TMP_PATH}/clu_switched" "$2" ${VERBOSITY} \
            || fail "mvdb switched clustering died"
        if [ -z "$KEEP_SWITCH_ALN" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${2}_aln" ${VERBOSITY}
        fi
        rm -rf "${TMP_PATH}/switch_tmp"
    fi

elif [ "$LINCLUST_MODULE" = "linclust1" ]; then
    # 0. clusthash
    if [ -n "$CLUSTHASH" ]; then
        if notExists "${TMP_PATH}/input_clusthash.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" clusthash "$INPUT" "${TMP_PATH}/input_clusthash" ${CLUSTHASH_PAR} \
                    || fail "clust hash died"
        fi

        if notExists "${TMP_PATH}/input_clusthash_clust.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" clust "$INPUT" "${TMP_PATH}/input_clusthash" "${TMP_PATH}/input_clusthash_clust" ${CLUSTHASH_CLUST_PAR} \
                    || fail "clust hash based clust died"
        fi

        awk '{print $1}' "${TMP_PATH}/input_clusthash_clust.index" > "${TMP_PATH}/order_clusthash_redundancy"

        if notExists "${TMP_PATH}/input_clusthash_redundancy.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/order_clusthash_redundancy" "$SOURCE" "${TMP_PATH}/input_clusthash_redundancy" ${VERBOSITY} --subdb-mode 1 \
                || fail "Createsubdb step died"
        fi
        INPUT="${TMP_PATH}/input_clusthash_redundancy"
    fi

    # 1. Finding exact k-mer matches.
    if notExists "${TMP_PATH}/pref.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
            || fail "kmermatcher died"
    fi

    # 2. Hamming distance pre-clustering
    if notExists "${TMP_PATH}/pref_rescore1.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/pref_rescore1" ${HAMMING_PAR} \
            || fail "Rescore with hamming distance step died"
    fi

    if notExists "${TMP_PATH}/pre_clust.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" clust "$INPUT" "${TMP_PATH}/pref_rescore1" "${TMP_PATH}/pre_clust" ${CLUSTER_PAR} \
            || fail "Pre-clustering step died"
    fi

    PRECLUST="${TMP_PATH}/pre_clust"
    # 0-2. clusthash (merge with preclust)
    if [ -n "$CLUSTHASH" ]; then
        if notExists "${TMP_PATH}/pre_clust_clusthash.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/pre_clust_clusthash" "${TMP_PATH}/input_clusthash_clust" "${TMP_PATH}/pre_clust" $MERGECLU_PAR \
                    || fail "mergeclusters clusthash based clust and preclust died"
        fi
        PRECLUST="${TMP_PATH}/pre_clust_clusthash"
    fi

    awk '{ print $1 }' "${PRECLUST}.index" > "${TMP_PATH}/order_redundancy"

    if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} --subdb-mode 1 \
            || fail "Createsubdb step died"
    fi

    if notExists "${TMP_PATH}/pref_filter1.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/order_redundancy" "${TMP_PATH}/pref" "${TMP_PATH}/pref_filter1" ${VERBOSITY} --subdb-mode 1 \
            || fail "Createsubdb step died"
    fi

    if notExists "${TMP_PATH}/pref_filter2.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/pref_filter1" "${TMP_PATH}/pref_filter2" --filter-file "${TMP_PATH}/order_redundancy" ${VERBOSITYANDCOMPRESS} \
            || fail "Filterdb step died"
    fi

    INPUT="${TMP_PATH}/input_step_redundancy"
    # 3. Ungapped alignment filtering
    RESULTDB="${TMP_PATH}/pref_filter2"
    if [ -n "$FILTER" ]; then
        if notExists "${TMP_PATH}/pref_rescore2.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "$RESULTDB" "${TMP_PATH}/pref_rescore2" ${UNGAPPED_ALN_PAR} \
                || fail "Ungapped alignment step died"
        fi
        RESULTDB="${TMP_PATH}/pref_rescore2"
    fi

    # 4. Local gapped sequence alignment.
    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$INPUT" "$RESULTDB" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
            || fail "Alignment step died"
    fi
    RESULTDB="${TMP_PATH}/aln"

    # 5. Clustering using greedy set cover.
    if notExists "${TMP_PATH}/clust.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" clust "$INPUT" "$RESULTDB" "${TMP_PATH}/clust" ${CLUSTER_PAR} \
            || fail "Clustering step died"
    fi

    if notExists "$2.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergeclusters "$SOURCE" "$2" "$PRECLUST" "${TMP_PATH}/clust" $MERGECLU_PAR \
            || fail "mergeclusters died"
    fi

fi

if [ -n "$REMOVE_TMP" ]; then
    if [ "$LINCLUST_MODULE" = "linclust2" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/clu" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_rep" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_rep_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_rep" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/clu_rep" ${VERBOSITY}
        # align intermediates (only present with --include-align-files)
        if [ -f "${TMP_PATH}/clu_aln.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_aln" ${VERBOSITY}
        fi
        if [ -f "${TMP_PATH}/clu_rep_aln.dbtype" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_rep_aln" ${VERBOSITY}
        fi
        if [ -n "$CLUSTHASH" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_clust" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_redundancy" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/clu_merged" ${VERBOSITY}
            rm -f "${TMP_PATH}/order_clusthash_redundancy"
        fi
        rm -f "${TMP_PATH}/linclust.sh"
    elif [ "$LINCLUST_MODULE" = "linclust1" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_filter1" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore1" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pre_clust" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy_h" ${VERBOSITY}
        rm -f "${TMP_PATH}/order_redundancy"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_filter2" ${VERBOSITY}
        if [ -n "$FILTER" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore2" ${VERBOSITY}
        fi
        if [ -n "$CLUSTHASH" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_clust" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/input_clusthash_redundancy" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/pre_clust_clusthash" ${VERBOSITY}
            rm -f "${TMP_PATH}/order_clusthash_redundancy"
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/clust" ${VERBOSITY}
        rm -f "${TMP_PATH}/linclust.sh"
    fi
fi