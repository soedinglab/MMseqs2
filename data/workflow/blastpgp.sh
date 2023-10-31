#!/bin/sh -e
# Iterative sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
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

QUERYDB="$1"
TMP_PATH="$4"

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ "$STEP" -lt "$NUM_IT" ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_tmp_${STEP}.done"; then
        if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
            TMP=""
            PREF="fake_pref"
        elif [ "$PREFMODE" = "UNGAPPED" ]; then
            PARAM="UNGAPPEDPREFILTER_PAR_$STEP"
            eval TMP="\$$PARAM"
            PREF="${MMSEQS} ungappedprefilter"
        else
            PARAM="PREFILTER_PAR_$STEP"
            eval TMP="\$$PARAM"
            PREF="${MMSEQS} prefilter"
        fi

        if [ "$STEP" -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER $PREF "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" ${TMP} \
                || fail "Prefilter died"
        else
            # shellcheck disable=SC2086
            $RUNNER $PREF "$QUERYDB" "$2" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
                || fail "Prefilter died"
        fi
        touch "$TMP_PATH/pref_tmp_${STEP}.done"
    fi

    if [ "$STEP" -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/pref_$STEP" $SUBSTRACT_PAR \
                || fail "Substract died"
            "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
        fi
        touch "$TMP_PATH/pref_$STEP.done"
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_tmp_$STEP.done"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"

        if [ "$STEP" -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" ${TMP} \
                || fail "Alignment died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_tmp_$STEP" ${TMP} \
                || fail "Alignment died"
        fi
        touch "$TMP_PATH/aln_tmp_$STEP.done"
    fi

    if [ "$STEP" -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.done"; then
            STEPONE=$((STEP-1))
            if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
                "$MMSEQS" mergedbs "$QUERYDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                    || fail "Alignment died"
            else
                "$MMSEQS" mergedbs "$QUERYDB" "$3" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                    || fail "Alignment died"
            fi
            "$MMSEQS" rmdb "$TMP_PATH/aln_$STEPONE"
            "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
            touch "$TMP_PATH/aln_$STEP.done"
        fi
    fi

# create profiles
    if [ "$STEP" -ne "$((NUM_IT  - 1))" ]; then
        if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" result2profile "$QUERYDB" "$2" "$TMP_PATH/aln_$STEP" "$TMP_PATH/profile_$STEP" ${TMP} \
                || fail "Create profile died"
        fi
    fi
	QUERYDB="$TMP_PATH/profile_$STEP"
	STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "$STEP" -lt "$NUM_IT" ]; do
        if [ "$STEP" -gt 0 ]; then
            rm -f -- "$TMP_PATH/aln_$STEP.done" "$TMP_PATH/pref_$STEP.done"
        fi
        rm -f -- "$TMP_PATH/aln_tmp_$STEP.done" "$TMP_PATH/pref_tmp_${STEP}.done"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/blastpgp.sh"
fi

