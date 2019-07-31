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
SOURCE="$INPUT"

mkdir -p "${TMP_PATH}/linclust"
if notExists "${TMP_PATH}/clu_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" linclust "$INPUT" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/linclust" ${LINCLUST_PAR} \
        || fail "linclust died"
fi

if notExists "${TMP_PATH}/input_step_redundancy.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" ${VERBOSITY} \
        || faill "createsubdb died"
fi

INPUT="${TMP_PATH}/input_step_redundancy"
STEP=0
STEPS=${STEPS:-1}
CLUSTER_STR=""
while [ "$STEP" -lt "$STEPS" ]; do
    PARAM=PREFILTER${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/pref_step$STEP.dbtype"; then
         # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "$INPUT" "$INPUT" "${TMP_PATH}/pref_step$STEP" ${TMP} \
            || fail "Prefilter step $STEP died"
    fi
    PARAM=ALIGNMENT${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/aln_step$STEP.dbtype"; then
         # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$INPUT" "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/aln_step$STEP" ${TMP} \
            || fail "Alignment step $STEP died"
    fi
    PARAM=CLUSTER${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/clu_step$STEP.dbtype"; then
         # shellcheck disable=SC2086
        "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/clu_step$STEP" ${TMP} \
            || fail "Clustering step $STEP died"
    fi

    # FIXME: This won't work if paths contain spaces
    CLUSTER_STR="${CLUSTER_STR} ${TMP_PATH}/clu_step$STEP"
    NEXTINPUT="${TMP_PATH}/input_step$((STEP+1))"
    if [ "$STEP" -eq "$((STEPS-1))" ]; then
       if [ -n "$REASSIGN" ]; then
          if notExists "${TMP_PATH}/clu.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} \
            || fail "Merging of clusters has died"
          fi
       else
            # shellcheck disable=SC2086
            "$MMSEQS" mergeclusters "$SOURCE" "$2" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} $MERGECLU_PAR \
            || fail "Merging of clusters has died"
       fi
    else
        if notExists "$NEXTINPUT.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "$INPUT" "$NEXTINPUT" ${VERBOSITY} \
                || fail "Order step $STEP died"
        fi
    fi

	INPUT="$NEXTINPUT"
	STEP=$((STEP+1))
done

if [ -n "$REASSIGN" ]; then
    STEP=$((STEP-1))
    PARAM=ALIGNMENT${STEP}_PAR
    eval ALIGNMENT_PAR="\$$PARAM"
    # align to cluster sequences
    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$SOURCE" "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/aln" ${ALIGNMENT_REASSIGN_PAR} \
                 || fail "align1 reassign died"
    fi
    # create file of cluster that do not align based on given criteria
    if notExists "${TMP_PATH}/clu_not_accepted.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/aln" "${TMP_PATH}/clu_not_accepted" --e-profile 100000000 ${THREADSANDCOMPRESS} \
                 || fail "subtractdbs1 reassign died"
    fi
    # create file of cluster that do align based on given criteria
    if notExists "${TMP_PATH}/clu_accepted.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_accepted" --e-profile 100000000 ${THREADSANDCOMPRESS} \
                 || fail "subtractdbs2 reassign died"
    fi
    if notExists "${TMP_PATH}/clu_not_accepted_swap.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_not_accepted_swap" ${THREADSANDCOMPRESS} \
                 || fail "swapdb1 reassign died"
    fi
    # create sequences database that were wrong assigned
    if notExists "${TMP_PATH}/seq_wrong_assigned.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/clu_not_accepted_swap" "$SOURCE" "${TMP_PATH}/seq_wrong_assigned" ${VERBOSITY} \
                 || fail "createsubdb1 reassign died"
    fi
    # build seed sequences
    if notExists "${TMP_PATH}/seq_seeds.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/clu" "$SOURCE" "${TMP_PATH}/seq_seeds" ${VERBOSITY} \
                || fail "createsubdb2 reassign died"
    fi
    PARAM=PREFILTER${STEP}_PAR
    eval PREFILTER_PAR="\$$PARAM"
    # try to find best matching centroid sequences for prev. wrong assigned sequences
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${TMP_PATH}/seq_wrong_assigned" "${TMP_PATH}/seq_seeds" "${TMP_PATH}/seq_wrong_assigned_pref" ${PREFILTER_PAR} \
                 || fail "Prefilter reassign died"
    fi
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref" "${TMP_PATH}/seq_wrong_assigned_pref_swaped" ${THREADSANDCOMPRESS} \
                 || fail "swapdb2 reassign died"
    fi
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${TMP_PATH}/seq_seeds" "${TMP_PATH}/seq_wrong_assigned" \
                                            "${TMP_PATH}/seq_wrong_assigned_pref_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${ALIGNMENT_REASSIGN_PAR} \
                 || fail "align2 reassign died"
    fi
    # extract only top hit and trim to frist column
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped" ${THREADSANDCOMPRESS} \
                 || fail "swapdb3 reassign died"
    fi
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1" --extract-lines 1 ${THREADSANDCOMPRESS} \
                    || fail "filterdb1 reassign died"
    fi
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol" --trim-to-one-column ${THREADSANDCOMPRESS} \
                    || fail "filterdb2 reassign died"
    fi
    if notExists "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped" ${THREADSANDCOMPRESS} \
                        || fail "swapdb2 reassign died"
    fi
    # combine clusters
    # shellcheck disable=SC2086
    "$MMSEQS" mergedbs "$SOURCE" "$2" "${TMP_PATH}/clu_accepted" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped" ${VERBCOMPRESS} \
                             || fail "mergedbs reassign died"
    if [ -n "$REMOVE_TMP" ]; then
        echo "Remove temporary files"
        "$MMSEQS" rmdb "${TMP_PATH}/aln"
        "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted"
        "$MMSEQS" rmdb "${TMP_PATH}/clu_accepted"
        "$MMSEQS" rmdb "${TMP_PATH}/clu_not_accepted_swap"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_seeds"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol"
        "$MMSEQS" rmdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped"
    fi
fi


if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/clu_redundancy"
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy"
    STEP=0
    while [ "$STEP" -lt "$STEPS" ]; do
        "$MMSEQS" rmdb "${TMP_PATH}/pref_step$STEP"
        "$MMSEQS" rmdb "${TMP_PATH}/aln_step$STEP"
        "$MMSEQS" rmdb "${TMP_PATH}/clu_step$STEP"
        STEP=$((STEP+1))
    done

    STEP=1
    while [ "$STEP" -lt "$STEPS" ]; do
        "$MMSEQS" rmdb "${TMP_PATH}/input_step$STEP"
        STEP=$((STEP+1))
    done

    rm -f "${TMP_PATH}/cascaded_clustering.sh"
fi

