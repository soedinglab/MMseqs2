#!/bin/bash
# Clustering workflow script
checkReturnCode () {
	if [ $? -ne 0 ]; then
	    echo "$1"
	    exit 1
	fi
}
notExists () {
	[ ! -f "$1" ]
}
#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[ ! -f "$TARGET_DB_PREF" ] &&  echo "$TARGET_DB_PREF not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $TMP_PATH not found!" && exit 1;

export OMP_PROC_BIND=TRUE

cd $(dirname "$1")
QUERY_FILE="$(basename $1)"
ABS_QUERY="$(pwd)/${QUERY_FILE}"
cd -

cd "$4"
TMP_PATH="$(pwd)"
cd -

STEP=0
QUERYDB="$ABS_QUERY"
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_$STEP"; then
        $RUNNER $MMSEQS prefilter "$QUERYDB" "$TARGET_DB_PREF" "$TMP_PATH/pref_$STEP" $PREFILTER_PAR
        checkReturnCode "Prefilter died"
    fi

    if [ $STEP -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.hasnext"; then
            # pref -aln
            $MMSEQS subtractdbs "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_0" "$TMP_PATH/pref_next_$STEP" $SUBSTRACT_PAR
            checkReturnCode "Substract died"
            mv -f "$TMP_PATH/pref_next_$STEP" "$TMP_PATH/pref_$STEP"
            mv -f "$TMP_PATH/pref_next_$STEP.index" "$TMP_PATH/pref_$STEP.index"
            touch "$TMP_PATH/pref_$STEP.hasnext"
        fi
    fi
    echo "RUN alignmment"

	REALIGN=""
	if [ $STEP -eq 0 ] && [ $PROFILE -eq 0 ]; then
	    REALIGN="--realign"
	fi
	# call alignment module
	if notExists "$TMP_PATH/aln_$STEP"; then
        $RUNNER $MMSEQS align "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" $ALIGNMENT_PAR $REALIGN -a
        checkReturnCode "Alignment died"
    fi

    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.hasmerge"; then
            $MMSEQS mergedbs "$QUERYDB" "$TMP_PATH/aln_new" "$TMP_PATH/aln_0" "$TMP_PATH/aln_$STEP"
            checkReturnCode "Merge died"
            mv -f "$TMP_PATH/aln_new" "$TMP_PATH/aln_0"
            mv -f "$TMP_PATH/aln_new.index" "$TMP_PATH/aln_0.index"
            touch "$TMP_PATH/aln_$STEP.hasmerge"
        fi
    fi

# create profiles
    if [ $STEP -ne $((NUM_IT  - 1)) ]; then
        if notExists "$TMP_PATH/profile_$STEP"; then
            $RUNNER $MMSEQS result2profile "$QUERYDB" "$2" "$TMP_PATH/aln_0" "$TMP_PATH/profile_$STEP" $PROFILE_PAR
            checkReturnCode "Create profile died"
            ln -sf "${QUERYDB}_h" "$TMP_PATH/profile_${STEP}_h"
            ln -sf "${QUERYDB}_h.index" "$TMP_PATH/profile_${STEP}_h.index"
        fi
    fi
	QUERYDB="$TMP_PATH/profile_$STEP"
    if [ $STEP -eq 0 ] && [ $PROFILE -eq 0 ]; then
        PREFILTER_PAR="$PREFILTER_PAR --profile"
        ALIGNMENT_PAR="$ALIGNMENT_PAR --profile"
        PROFILE_PAR="$PROFILE_PAR --profile"
    fi
	let STEP=STEP+1
done
# post processing
let STEP=STEP-1
mv -f "$TMP_PATH/aln_0" "$3"
mv -f "$TMP_PATH/aln_0.index" "$3.index"
checkReturnCode "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 STEP=0
 while [ $STEP -lt $NUM_IT ]; do
    rm -f "$TMP_PATH/pref_$STEP" "$TMP_PATH/pref_$STEP.index"
    rm -f "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEP.index"
    rm -f "$TMP_PATH/profile_$STEP" "$TMP_PATH/profile_$STEP.index" "$TMP_PATH/profile_${STEP}_h" "$TMP_PATH/profile_${STEP}_h.index"
    let STEP=STEP+1
 done
fi

