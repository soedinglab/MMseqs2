#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
SOURCE="$INPUT"

mkdir -p "${TMP_PATH}/linclust"
if notExists "${TMP_PATH}/clu_redundancy"; then
    # shellcheck disable=SC2086
    "$MMSEQS" linclust "$INPUT" "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/linclust" ${LINCLUST_PAR} \
        || fail "linclust died"
fi

if notExists "${TMP_PATH}/input_step_redundancy"; then
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_redundancy" "$INPUT" "${TMP_PATH}/input_step_redundancy" \
        || faill "createsubdb died"
fi

INPUT="${TMP_PATH}/input_step_redundancy"
STEP=0
STEPS=${STEPS:-1}
CLUSTER_STR=""
while [ "$STEP" -lt "$STEPS" ]; do
    PARAM=PREFILTER${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/pref_step$STEP"; then
         # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "$INPUT" "$INPUT" "${TMP_PATH}/pref_step$STEP" ${TMP} \
            || fail "Prefilter step $STEP died"
    fi
    PARAM=ALIGNMENT${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/aln_step$STEP"; then
         # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$INPUT" "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/aln_step$STEP" ${TMP} \
            || fail "Alignment step $STEP died"
    fi
    PARAM=CLUSTER${STEP}_PAR
    eval TMP="\$$PARAM"
    if notExists "${TMP_PATH}/clu_step$STEP"; then
         # shellcheck disable=SC2086
        "$MMSEQS" clust "$INPUT" "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/clu_step$STEP" ${TMP} \
            || fail "Clustering step $STEP died"
    fi

    # FIXME: This won't work if paths contain spaces
    CLUSTER_STR="${CLUSTER_STR} ${TMP_PATH}/clu_step$STEP"
    NEXTINPUT="${TMP_PATH}/input_step$((STEP+1))"
    if [ "$STEP" -eq "$((STEPS-1))" ]; then
        if notExists "${TMP_PATH}/clu"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergeclusters "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/clu_redundancy" ${CLUSTER_STR} \
                || fail "Merging of clusters has died"
        fi
    else
        if notExists "$NEXTINPUT"; then
            "$MMSEQS" createsubdb "${TMP_PATH}/clu_step$STEP" "$INPUT" "$NEXTINPUT" \
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
    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$SOURCE" "$SOURCE" "${TMP_PATH}/clu" "${TMP_PATH}/aln" ${ALIGNMENT_REASSIGN_PAR} \
             || fail "align1 reassign died"
    "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/aln" "${TMP_PATH}/clu_not_accepted" --e-profile 100000 \
             || fail "subtractdbs1 reassign died"
    "$MMSEQS" subtractdbs "${TMP_PATH}/clu" "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_accepted" --e-profile 100000 \
             || fail "subtractdbs2 reassign died"
    "$MMSEQS" swapdb "${TMP_PATH}/clu_not_accepted" "${TMP_PATH}/clu_not_accepted_swap" \
             || fail "swapdb1 reassign died"
    "$MMSEQS" createsubdb "${TMP_PATH}/clu_not_accepted_swap" "$SOURCE" "${TMP_PATH}/seq_wrong_assigned" \
             || fail "createsubdb1 reassign died"
    "$MMSEQS" createsubdb "${TMP_PATH}/clu" "$SOURCE" "${TMP_PATH}/seq_seeds"  \
             || fail "createsubdb2 reassign died"
    PARAM=PREFILTER${STEP}_PAR
    eval PREFILTER_PAR="\$$PARAM"
    $RUNNER "$MMSEQS" prefilter  "${TMP_PATH}/seq_wrong_assigned" "${TMP_PATH}/seq_seeds" "${TMP_PATH}/seq_wrong_assigned_pref" ${PREFILTER_PAR} \
             || fail "Prefilter reassign died"
    "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref" "${TMP_PATH}/seq_wrong_assigned_pref_swaped" \
             || fail "swapdb2 reassign died"

    $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${TMP_PATH}/seq_seeds" "${TMP_PATH}/seq_wrong_assigned" \
                                        "${TMP_PATH}/seq_wrong_assigned_pref_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" ${ALIGNMENT_REASSIGN_PAR} \
             || fail "align2 reassign died"
    "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped"  \
             || fail "swapdb3 reassign died"

    "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1" --extract-lines 1 \
                || fail "filterdb1 reassign died"
    "$MMSEQS" filterdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol" --trim-to-one-column \
                    || fail "filterdb2 reassign died"
    "$MMSEQS" swapdb "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped" \
                        || fail "swapdb2 reassign died"
    "$MMSEQS" mergedbs "$SOURCE" "${TMP_PATH}/clu_reassign" "${TMP_PATH}/clu_accepted" "${TMP_PATH}/seq_wrong_assigned_pref_swaped_aln_swaped_top1_ocol_swaped" \
                             || fail "mergedbs reassign died"
            # post processing
    mv -f "${TMP_PATH}/clu_reassign" "$2" || fail "Could not move result to $2"
    mv -f "${TMP_PATH}/clu_reassign.index" "$2.index" || fail "Could not move result to $2"
else
    # post processing
    mv -f "${TMP_PATH}/clu" "$2" || fail "Could not move result to $2"
    mv -f "${TMP_PATH}/clu.index" "$2.index" || fail "Could not move result to $2"
fi


if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "${TMP_PATH}/order_redundancy"
 rm -f "${TMP_PATH}/clu_redundancy" "${TMP_PATH}/clu_redundancy.index"
 rm -f "${TMP_PATH}/aln_redundancy" "${TMP_PATH}/aln_redundancy.index"
 rm -f "${TMP_PATH}/input_step_redundancy" "${TMP_PATH}/input_step_redundancy.index"
 STEP=0
 while [ "$STEP" -lt "$STEPS" ]; do
    rm -f "${TMP_PATH}/pref_step$STEP" "${TMP_PATH}/pref_step$STEP.index"
    rm -f "${TMP_PATH}/aln_step$STEP" "${TMP_PATH}/aln_step$STEP.index"
    rm -f "${TMP_PATH}/clu_step$STEP" "${TMP_PATH}/clu_step$STEP.index"
    rm -f "${TMP_PATH}/input_step$STEP" "${TMP_PATH}/input_step$STEP.index"
    rm -f "${TMP_PATH}/order_step$STEP"
	STEP=$((STEP+1))
 done

 rm -f "${TMP_PATH}/cascaded_clustering.sh"
fi
