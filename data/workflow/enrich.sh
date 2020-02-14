#!/bin/sh -ex
# Iterative sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "${MMSEQS}" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 6 ] && echo "Please provide <queryDB> <targetDB> <targetProf> <targetRes> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1;
[ ! -f "$4.dbtype" ] && echo "$4.dbtype not found!" && exit 1;
[   -f "$5.dbtype" ] && echo "$5.dbtype exists already!" && exit 1;
[ ! -d "$6" ] && echo "tmp directory $6 not found!" && mkdir -p "$6";

QUERYDB="$1"
PROFTARGETSEQ="$2"
TARGETPROF="$3"
PROFRESULT="$4"
RESULT="$5"
TMP_PATH="$6"

if notExists "${TMP_PATH}/search_slice"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERYDB}" "${TARGETPROF}" "${TMP_PATH}/search_slice" "${TMP_PATH}/slice_tmp" ${PROF_SEARCH_PAR} \
        || fail "search died"
fi

if notExists "${TMP_PATH}/prof_slice"; then
    # shellcheck disable=SC2086
    ${RUNNER} "${MMSEQS}" result2profile "${QUERYDB}" "${TARGETPROF}" "${TMP_PATH}/search_slice" "${TMP_PATH}/prof_slice" ${PROF_PROF_PAR} \
        || fail "result2profile died"
fi

INPUT="${TMP_PATH}/prof_slice"
STEP=0
while [ "${STEP}" -lt "${NUM_IT}" ]; do
    # call prefilter module
    if notExists "${TMP_PATH}/pref_${STEP}"; then
        PARAM="PREFILTER_PAR_${STEP}"
        eval TMP="\$$PARAM"
        # shellcheck disable=SC2086
        ${RUNNER} "${MMSEQS}" prefilter "${INPUT}" "${TARGETPROF}_consensus" "${TMP_PATH}/pref_${STEP}" ${TMP} \
            || fail "prefilter died"
    fi

    if [ ${STEP} -ge 1 ]; then
        if notExists "${TMP_PATH}/pref_${STEP}.hasnext"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" subtractdbs "${TMP_PATH}/pref_${STEP}" "${TMP_PATH}/aln_0" "${TMP_PATH}/pref_next_${STEP}" ${SUBSTRACT_PAR} \
                || fail "subtractdbs died"
            mv -f "${TMP_PATH}/pref_next_${STEP}" "${TMP_PATH}/pref_${STEP}"
            mv -f "${TMP_PATH}/pref_next_${STEP}.index" "${TMP_PATH}/pref_${STEP}.index"
            mv -f "${TMP_PATH}/pref_next_${STEP}.dbtype" "${TMP_PATH}/pref_${STEP}.dbtype"
            touch "${TMP_PATH}/pref_${STEP}.hasnext"
        fi
    fi

	# call alignment module
	if notExists "${TMP_PATH}/aln_${STEP}"; then
	    PARAM="ALIGNMENT_PAR_${STEP}"
        eval TMP="\$$PARAM"
        # shellcheck disable=SC2086
        ${RUNNER} "${MMSEQS}" "${ALIGN_MODULE}" "${INPUT}" "${TARGETPROF}_consensus" "${TMP_PATH}/pref_${STEP}" "${TMP_PATH}/aln_${STEP}" ${TMP} \
            || fail "${ALIGN_MODULE} died"
    fi


    if notExists "${TMP_PATH}/aln_${STEP}.hasexpand"; then
        PARAM="EXPAND_PAR_${STEP}"
        eval TMP="\$$PARAM"
        # shellcheck disable=SC2086
        "${MMSEQS}" expandaln "${INPUT}" "${PROFTARGETSEQ}" "${TMP_PATH}/aln_${STEP}" "${PROFRESULT}" "${TMP_PATH}/aln_exp_${STEP}" ${TMP} \
            || fail "expandaln died"
        mv -f "${TMP_PATH}/aln_exp_${STEP}" "${TMP_PATH}/aln_${STEP}"
        mv -f "${TMP_PATH}/aln_exp_${STEP}.index" "${TMP_PATH}/aln_${STEP}.index"
        mv -f "${TMP_PATH}/aln_exp_${STEP}.dbtype" "${TMP_PATH}/aln_${STEP}.dbtype"
        touch "${TMP_PATH}/aln_exp_${STEP}.hasexpand"
    fi

    if [ ${STEP} -gt 0 ]; then
        if notExists "${TMP_PATH}/aln_${STEP}.hasmerge"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" mergedbs "${INPUT}" "${TMP_PATH}/aln_new" "${TMP_PATH}/aln_0" "${TMP_PATH}/aln_${STEP}" ${VERBOSITY_PAR} \
                || fail "mergedbs died"
            mv -f "${TMP_PATH}/aln_new" "${TMP_PATH}/aln_0"
            mv -f "${TMP_PATH}/aln_new.index" "${TMP_PATH}/aln_0.index"
            mv -f "${TMP_PATH}/aln_new.dbtype" "${TMP_PATH}/aln_0.dbtype"
            touch "${TMP_PATH}/aln_${STEP}.hasmerge"
        fi
    fi

    # create profiles
    if [ "$((STEP-1))" != "${NUM_IT}" ]; then
        if notExists "${TMP_PATH}/profile_${STEP}"; then
            PARAM="PROFILE_PAR_${STEP}"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            ${RUNNER} "${MMSEQS}" result2profile "${QUERYDB}" "${PROFTARGETSEQ}" "${TMP_PATH}/aln_0" "${TMP_PATH}/profile_${STEP}" ${TMP} \
                || fail "result2profile died"
        fi
        INPUT="${TMP_PATH}/profile_${STEP}"
	fi

	STEP="$((STEP+1))"
done

mv -f "${TMP_PATH}/aln_0" "${RESULT}"
mv -f "${TMP_PATH}/aln_0.index" "${RESULT}.index"
mv -f "${TMP_PATH}/aln_0.dbtype" "${RESULT}.dbtype"

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "${STEP}" -lt "${NUM_IT}" ]; do
        rm -f "${TMP_PATH}/pref_${STEP}" "${TMP_PATH}/pref_${STEP}.index" "${TMP_PATH}/pref_${STEP}.dbtype"
        rm -f "${TMP_PATH}/aln_${STEP}" "${TMP_PATH}/aln_${STEP}.index" "${TMP_PATH}/aln_${STEP}.dbtype"
        rm -f "${TMP_PATH}/profile_${STEP}" "${TMP_PATH}/profile_${STEP}.index" "${TMP_PATH}/profile_${STEP}.dbtype"
        rm -f "${TMP_PATH}/profile_${STEP}_h" "${TMP_PATH}/profile_${STEP}_h.index" "${TMP_PATH}/profile_${STEP}_h.dbtype"
        rm -f "${TMP_PATH}/profile_${STEP}_consensus" "${TMP_PATH}/profile_${STEP}_consensus.index" "${TMP_PATH}/profile_${STEP}_consensus.dbtype"
        rm -f "${TMP_PATH}/profile_${STEP}_consensus_h" "${TMP_PATH}/profile_${STEP}_consensus_h.index" "${TMP_PATH}/profile_${STEP}_consensus_h.dbtype"
        rm -f "${TMP_PATH}/aln_${STEP}.hasmerge" "${TMP_PATH}/aln_exp_${STEP}.hasexpand" "${TMP_PATH}/pref_${STEP}.hasnext"
        STEP="$((STEP+1))"
    done
    rm -f "${TMP_PATH}/prof_slice" "${TMP_PATH}/prof_slice.index" "${TMP_PATH}/prof_slice.dbtype"
    rm -f "${TMP_PATH}/prof_slice_h" "${TMP_PATH}/prof_slice_h.index"
    rm -f "${TMP_PATH}/prof_slice_consensus" "${TMP_PATH}/prof_slice_consensus.index" "${TMP_PATH}/prof_slice_consensus.dbtype"
    rm -f "${TMP_PATH}/prof_slice_consensus_h" "${TMP_PATH}/prof_slice_consensus_h.index" "${TMP_PATH}/prof_slice_consensus_h.dbtype"
    rm -f "${TMP_PATH}/search_slice" "${TMP_PATH}/search_slice.index" "${TMP_PATH}/search_slice.dbtype"
    rm -f "${TMP_PATH}/enrich.sh"
fi

