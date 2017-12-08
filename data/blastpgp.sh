#!/bin/bash -e
# Iterative sequence search workflow script
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
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERYDB="$(abspath $1)"
TMP_PATH="$(abspath $4)"

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_$STEP"; then
        PARAM="PREFILTER_PAR_$STEP"
        eval TMP="\$$PARAM"
        $RUNNER $MMSEQS prefilter "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" ${TMP} \
            || fail "Prefilter died"
    fi

    if [ $STEP -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.hasnext"; then
            # pref -aln
            $MMSEQS subtractdbs "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_0" "$TMP_PATH/pref_next_$STEP" $SUBSTRACT_PAR \
                || fail "Substract died"
            mv -f "$TMP_PATH/pref_next_$STEP" "$TMP_PATH/pref_$STEP"
            mv -f "$TMP_PATH/pref_next_$STEP.index" "$TMP_PATH/pref_$STEP.index"
            touch "$TMP_PATH/pref_$STEP.hasnext"
        fi
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_$STEP"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"
        $RUNNER $MMSEQS align "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" ${TMP} \
            || fail "Alignment died"
    fi

    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.hasmerge"; then
            $MMSEQS mergedbs "$QUERYDB" "$TMP_PATH/aln_new" "$TMP_PATH/aln_0" "$TMP_PATH/aln_$STEP" \
                || fail "Merge died"
            mv -f "$TMP_PATH/aln_new" "$TMP_PATH/aln_0"
            mv -f "$TMP_PATH/aln_new.index" "$TMP_PATH/aln_0.index"
            touch "$TMP_PATH/aln_$STEP.hasmerge"
        fi
    fi

# create profiles
    if [ $STEP -ne $((NUM_IT  - 1)) ]; then
        if notExists "$TMP_PATH/profile_$STEP"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"
            $RUNNER $MMSEQS result2profile "$QUERYDB" "$2" "$TMP_PATH/aln_0" "$TMP_PATH/profile_$STEP" ${TMP} \
                || fail "Create profile died"
            ln -sf "${QUERYDB}_h" "$TMP_PATH/profile_${STEP}_h"
            ln -sf "${QUERYDB}_h.index" "$TMP_PATH/profile_${STEP}_h.index"
        fi
    fi
	QUERYDB="$TMP_PATH/profile_$STEP"

	STEP=$(($STEP+1))
done
# post processing
STEP=$(($STEP-1))
(mv -f "$TMP_PATH/aln_0" "$3" && mv -f "$TMP_PATH/aln_0.index" "$3.index") || fail "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 STEP=0
 while [ $STEP -lt $NUM_IT ]; do
    rm -f "$TMP_PATH/pref_$STEP" "$TMP_PATH/pref_$STEP.index"
    rm -f "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEP.index"
    rm -f "$TMP_PATH/profile_$STEP" "$TMP_PATH/profile_$STEP.index" "$TMP_PATH/profile_${STEP}_h" "$TMP_PATH/profile_${STEP}_h.index"
    STEP=$(($STEP+1))
 done

 rm -f "$TMP_PATH/blastpgp.sh"
fi

