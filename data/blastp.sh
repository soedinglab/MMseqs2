#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
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


INPUT="$1"
TARGET="$2"
TMP_PATH="$4"


STEP=0
STEPS=${STEPS:-1}
while [ "$STEP" -lt "$STEPS" ]; do
    SENS_PARAM=SENSE_${STEP}
    eval SENS="\$$SENS_PARAM"
    # call prefilter module
    if notExists "$TMP_PATH/pref_$SENS"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "$INPUT" "$TARGET" "$TMP_PATH/pref_$SENS" $PREFILTER_PAR -s "$SENS" \
            || fail "Prefilter died"
    fi

    # call alignment module
    if notExists "$TMP_PATH/aln_$SENS"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" align "$INPUT" "$TARGET${ALIGNMENT_DB_EXT}" "$TMP_PATH/pref_$SENS" "$TMP_PATH/aln_$SENS" $ALIGNMENT_PAR  \
            || fail "Alignment died"
    fi

    # only merge results after first step
    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_${SENS}.hasmerge"; then
            "$MMSEQS" mergedbs "$1" "$TMP_PATH/aln_new" "$TMP_PATH/aln_${SENSE_0}" "$TMP_PATH/aln_$SENS" \
                || fail "Alignment died"
            mv -f "$TMP_PATH/aln_new" "$TMP_PATH/aln_${SENSE_0}"
            mv -f "$TMP_PATH/aln_new.index" "$TMP_PATH/aln_${SENSE_0}.index"
            touch "$TMP_PATH/aln_${SENS}.hasmerge"
        fi
    fi

    NEXTINPUT="$TMP_PATH/input_step$SENS"
    #do not create subdb at last step
    if [ $STEP -lt $((STEPS-1)) ]; then
        if notExists "$TMP_PATH/order_step$SENS"; then
            awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$SENS.index" > "$TMP_PATH/order_step$SENS" \
                || fail "Awk step $SENS died"
        fi

        if [ ! -s "$TMP_PATH/order_step$SENS" ]; then break; fi

        if notExists "$NEXTINPUT"; then
            "$MMSEQS" createsubdb "$TMP_PATH/order_step$SENS" "$INPUT" "$NEXTINPUT" \
                || fail "Order step $SENS died"
        fi
    fi
    INPUT="$NEXTINPUT"
    STEP=$((STEP+1))
done

# post processing
(mv -f "$TMP_PATH/aln_${SENSE_0}" "$3" && mv -f "$TMP_PATH/aln_${SENSE_0}.index" "$3.index" ) \
    || fail "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    STEP=0
    while [ "$STEP" -lt "$STEPS" ]; do
        SENS_PARAM=SENSE_${STEP}
        eval SENS="\$$SENS_PARAM"
        rm -f "$TMP_PATH/pref_$SENS" "$TMP_PATH/pref_$SENS.index"
        rm -f "$TMP_PATH/aln_$SENS" "$TMP_PATH/aln_$SENS.index"
        NEXTINPUT="$TMP_PATH/input_step$SENS"
        rm -f "$TMP_PATH/input_step$SENS" "$TMP_PATH/input_step$SENS.index"
        STEP=$((STEP+1))
    done

    rm -f "$TMP_PATH/blastp.sh"
fi


