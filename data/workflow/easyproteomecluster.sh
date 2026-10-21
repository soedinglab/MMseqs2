#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
   [ ! -f "$1" ]
}


if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$@" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "query createdb died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${CLUSTER_MODULE}" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "linclust died"
fi

if notExists "${RESULTS}_protein_cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${RESULTS}_protein_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv protein cluster died"
fi

if notExists "${TMP_PATH}/aln_proteome.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" proteomecluster "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/aln_proteome" "${TMP_PATH}/cluster_count" "${TMP_PATH}/aln_protein" ${PROTEOMECLUSTER_PAR} \
        || fail "proteomecluster died"
fi

if notExists "${RESULTS}_cluster_count.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/cluster_count" "${RESULTS}_cluster_count.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster count report died"
fi

if notExists "${RESULTS}_protein_align.tsv" && [ -n "${WRITE_ALIGN_PROTEOME}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/aln_protein" "${RESULTS}_protein_align.tsv" ${THREADS_PAR} \
            || fail "createtsv protein align died"
else
    rm -rf "${TMP_PATH}/aln_protein"*
fi

# cascade 
awk 'NR==FNR { sub(/^\x00/, "", $1); a[$1]; next } !($1 in a)' "${TMP_PATH}/aln_proteome" "${TMP_PATH}/input.source" > "${TMP_PATH}/source_filtered"
SOURCEtoNEXTITERATION="${TMP_PATH}/source_filtered"
STEP=2

SUBDB_LOOKUP_LIST="${TMP_PATH}/input.lookup"

if [ -n "${REASSIGN}" ] && [ -n "${CASCADED_PROTEOME_CLUSTERING}" ] && notExists "${TMP_PATH}/aln_proteome_1.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cpdb "${TMP_PATH}/aln_proteome" "${TMP_PATH}/aln_proteome_1" ${VERBOSITY_PAR} \
        || fail "cpdb died"
fi

# Corrected while loop condition with proper spacing and quoting
while [ -s "$SOURCEtoNEXTITERATION" ] && [ -n "${CASCADED_PROTEOME_CLUSTERING}" ]; do
    echo "Step $STEP: $(wc -l < "$SOURCEtoNEXTITERATION") sources left"
    # Make "sublookup_STEP" from lines in input.lookup whose 3rd field is in the set from source_filtered
    awk 'NR==FNR {sources[$1]; next} $3 in sources' "$SOURCEtoNEXTITERATION" "${SUBDB_LOOKUP_LIST}" > "${TMP_PATH}/sublookup_${STEP}"

    # Create a smaller DB from sublookup
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/sublookup_${STEP}" "${TMP_PATH}/input" "${TMP_PATH}/input_${STEP}" --subdb-mode 1
    NEXTINPUT="${TMP_PATH}/input_${STEP}"

    # Run linclust on the newly created sub-DB
    echo "Run linclust for iter $STEP" 
    if notExists "${TMP_PATH}/clu_${STEP}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" "${CLUSTER_MODULE}" "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${TMP_PATH}/clu_tmp_${STEP}" ${CLUSTER_PAR} \
            || fail "linclust died"
    fi

    echo "Run createtsv: protein clust result for iter $STEP"
    if notExists "${RESULTS}_protein_cluster_${STEP}.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${RESULTS}_protein_cluster_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv protein cluster died"
    fi

    if [ -n "$REMOVE_TMP" ]; then
        # shellcheck disable=SC2086
        rm -rf "${TMP_PATH}/clu_tmp_${STEP}"
    fi

    echo "Run ProteomeCluster for iter $STEP"
    # Run proteomecluster on the newly created sub-DB
    if notExists "${TMP_PATH}/aln_proteome_${STEP}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" proteomecluster "${NEXTINPUT}" "${TMP_PATH}/clu_${STEP}" "${TMP_PATH}/aln_proteome_${STEP}" "${TMP_PATH}/cluster_count_${STEP}" "${TMP_PATH}/aln_protein_${STEP}" ${PROTEOMECLUSTER_PAR} \
            || fail "proteomecluster died"
    fi

    echo "Run createtsv: clustercount report for iter $STEP"
    if notExists "${RESULTS}_cluster_count_${STEP}.tsv"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${TMP_PATH}/cluster_count_${STEP}" "${RESULTS}_cluster_count_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv proteome cluster count report died"
    fi

    if notExists "${RESULTS}_protein_align_${STEP}.tsv" && [ -n "${WRITE_ALIGN_PROTEOME}" ]; then
        echo "Run createtsv: protein align result for iter $STEP"
        # shellcheck disable=SC2086
        "$MMSEQS" createtsv "${NEXTINPUT}" "${NEXTINPUT}" "${TMP_PATH}/aln_protein_${STEP}" "${RESULTS}_protein_align_${STEP}.tsv" ${THREADS_PAR} \
                || fail "createtsv protein align died"
    else
        rm -rf "${TMP_PATH}/aln_protein_${STEP}"*
    fi

    echo "Run concatdbs of aln_proteome for iter $STEP"
    # Concatenate new proteome alignments into the master aln_proteome
    "$MMSEQS" concatdbs "${TMP_PATH}/aln_proteome" "${TMP_PATH}/aln_proteome_${STEP}" "${TMP_PATH}/aln_proteome" --preserve-keys 1

    # Repeat the AWK-based filtering to update source_filtered
    awk 'NR==FNR { sub(/^\x00/, "", $1); a[$1]; next } !($1 in a)' "${TMP_PATH}/aln_proteome" "${TMP_PATH}/input.source" > "${TMP_PATH}/source_filtered"
    SUBDB_LOOKUP_LIST="${TMP_PATH}/sublookup_${STEP}"
    # rm -f "${TMP_PATH}/sublookup_${STEP}"
    STEP=$((STEP + 1))
done

ALN_PROTEOME="${TMP_PATH}/aln_proteome"
if [ -n "${REASSIGN}" ]; then
    CLUSTER_STR="${TMP_PATH}/aln_proteome"
    if [ -n "${CASCADED_PROTEOME_CLUSTERING}" ]; then
        CLUSTER_STR="${TMP_PATH}/aln_proteome_1"
        STEP=2
        while [ -f "${TMP_PATH}/aln_proteome_${STEP}.dbtype" ]; do
            CLUSTER_STR="${CLUSTER_STR} ${TMP_PATH}/aln_proteome_${STEP}"
            STEP=$((STEP + 1))
        done
    fi
    if notExists "${TMP_PATH}/input_reference.dbtype"; then
        awk 'NR == FNR { r[$1]; next } $3 in r { print $1 }' "${TMP_PATH}/aln_proteome.index" "${TMP_PATH}/input.lookup" > "${TMP_PATH}/reference_proteins"
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/reference_proteins" "${TMP_PATH}/input" "${TMP_PATH}/input_reference" ${VERBOSITY_PAR} --subdb-mode 1 \
            || fail "createsubdb reference proteins died"
    fi
    REFERENCE_SEARCH_PAR="${SEARCH_PAR}"
    CANDIDATE_SEARCH_PAR="${SEARCH_PAR}"
    if [ -z "${MAX_SEQS_SET}" ]; then
        REFERENCE_SEARCH_PAR="${SEARCH_PAR} --max-seqs $(awk 'END { print (NR > 300 ? NR : 300) }' "${TMP_PATH}/aln_proteome.index")"
        CANDIDATE_SEARCH_PAR="${SEARCH_PAR} --max-seqs $(awk 'FILENAME == ARGV[1] { n++; next } END { m = FNR - n; print (m > 300 ? m : 300) }' "${TMP_PATH}/aln_proteome.index" "${TMP_PATH}/input.source")"
    fi
    if notExists "${TMP_PATH}/aln_reference.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" search "${TMP_PATH}/input_reference" "${TMP_PATH}/input_reference" "${TMP_PATH}/aln_reference" "${TMP_PATH}/search_reference_tmp" ${REFERENCE_SEARCH_PAR} \
            || fail "search reference proteomes died"
    fi
    if notExists "${TMP_PATH}/aln_reference_swap.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/aln_reference" "${TMP_PATH}/aln_reference_swap" ${THREADS_PAR} \
            || fail "swapdb reference proteomes died"
    fi
    if notExists "${TMP_PATH}/score_reference.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" parseproteomealignments "${TMP_PATH}/input_reference" "${TMP_PATH}/input_reference" "${TMP_PATH}/aln_reference_swap" "${TMP_PATH}/score_reference" ${THREADS_PAR} \
            || fail "parseproteomealignments reference proteomes died"
    fi
    if notExists "${TMP_PATH}/score_reference_pass.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/score_reference" "${TMP_PATH}/score_reference_pass" --filter-column 2 --comparison-operator ge --comparison-value "${PROTEOME_SIMILARITY}" ${THREADS_PAR} \
            || fail "filterdb reference proteomes died"
    fi
    # proteins of the references a passing group may join, and of the members of such groups
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "${TMP_PATH}/score_reference_pass" "${TMP_PATH}/score_reference_pass.tsv" --tsv ${THREADS_PAR} \
        || fail "prefixid reference scores died"
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "${TMP_PATH}/aln_proteome" "${TMP_PATH}/aln_proteome.tsv" --tsv ${THREADS_PAR} \
        || fail "prefixid proteome clusters died"
    awk -v REFERENCES="${TMP_PATH}/candidate_reference_proteins" -v MEMBERS="${TMP_PATH}/candidate_member_proteins" \
        'FILENAME == ARGV[1] { if ($1 != $2) { group[$1]; reference[$2] } next }
         FILENAME == ARGV[2] { if (($1 in group) && $2 != $1) { member[$2] } next }
         $3 in reference { print $1 > REFERENCES } $3 in member { print $1 > MEMBERS }' \
        "${TMP_PATH}/score_reference_pass.tsv" "${TMP_PATH}/aln_proteome.tsv" "${TMP_PATH}/input.lookup"
    if [ -s "${TMP_PATH}/candidate_reference_proteins" ]; then
        PROTEOME_SCORE="${TMP_PATH}/score_reference_pass"
        if [ -s "${TMP_PATH}/candidate_member_proteins" ]; then
            if notExists "${TMP_PATH}/input_candidate_member.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "${TMP_PATH}/candidate_reference_proteins" "${TMP_PATH}/input" "${TMP_PATH}/input_candidate_reference" ${VERBOSITY_PAR} --subdb-mode 1 \
                    || fail "createsubdb candidate references died"
                # shellcheck disable=SC2086
                "$MMSEQS" createsubdb "${TMP_PATH}/candidate_member_proteins" "${TMP_PATH}/input" "${TMP_PATH}/input_candidate_member" ${VERBOSITY_PAR} --subdb-mode 1 \
                    || fail "createsubdb candidate members died"
            fi
            if notExists "${TMP_PATH}/aln_candidate.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" search "${TMP_PATH}/input_candidate_reference" "${TMP_PATH}/input_candidate_member" "${TMP_PATH}/aln_candidate" "${TMP_PATH}/search_candidate_tmp" ${CANDIDATE_SEARCH_PAR} \
                    || fail "search candidate members died"
            fi
            if notExists "${TMP_PATH}/aln_candidate_swap.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" swapdb "${TMP_PATH}/aln_candidate" "${TMP_PATH}/aln_candidate_swap" ${THREADS_PAR} \
                    || fail "swapdb candidate members died"
            fi
            if notExists "${TMP_PATH}/score_candidate.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" parseproteomealignments "${TMP_PATH}/input_candidate_member" "${TMP_PATH}/input_candidate_reference" "${TMP_PATH}/aln_candidate_swap" "${TMP_PATH}/score_candidate" ${THREADS_PAR} \
                    || fail "parseproteomealignments candidate members died"
            fi
            if notExists "${TMP_PATH}/score_proteome.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" concatdbs "${TMP_PATH}/score_reference_pass" "${TMP_PATH}/score_candidate" "${TMP_PATH}/score_proteome" --preserve-keys 1 ${THREADS_PAR} \
                    || fail "concatdbs proteome scores died"
            fi
            PROTEOME_SCORE="${TMP_PATH}/score_proteome"
        fi
        if notExists "${TMP_PATH}/aln_proteome_merged.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeproteomeclusters "${TMP_PATH}/input" "${PROTEOME_SCORE}" "${TMP_PATH}/aln_proteome_merged" ${CLUSTER_STR} ${MERGEPROTEOMECLUSTERS_PAR} \
                || fail "mergeproteomeclusters died"
        fi
        ALN_PROTEOME="${TMP_PATH}/aln_proteome_merged"
    fi
fi

if notExists "${RESULTS}_proteome_cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${ALN_PROTEOME}" "${RESULTS}_proteome_cluster.tsv" ${THREADS_PAR} \
            || fail "createtsv proteome cluster died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    for DB in aln_proteome_1 input_reference input_reference_h aln_reference aln_reference_swap score_reference score_reference_pass \
              input_candidate_reference input_candidate_reference_h input_candidate_member input_candidate_member_h \
              aln_candidate aln_candidate_swap score_candidate score_proteome aln_proteome_merged; do
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/${DB}" ${VERBOSITY_PAR}
    done
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_protein" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_proteome" ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/clu_tmp" "${TMP_PATH}/search_reference_tmp" "${TMP_PATH}/search_candidate_tmp"
    rm -f "${TMP_PATH}/reference_proteins" "${TMP_PATH}/candidate_reference_proteins" "${TMP_PATH}/candidate_member_proteins" \
          "${TMP_PATH}/score_reference_pass.tsv" "${TMP_PATH}/aln_proteome.tsv"
    rm -f "${TMP_PATH}/easyproteomecluster.sh"
fi