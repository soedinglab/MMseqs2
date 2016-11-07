#!/bin/bash -e
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
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && exit 1;

export OMP_PROC_BIND=TRUE

cd $(dirname $1)
QUERY_FILE=$(basename $1)
ABS_QUERY="$(pwd)/${QUERY_FILE}"
cd -

cd $4
TMP_PATH=$(pwd)
cd -
INPUT=$1
TARGET=$2
SENS=$START_SENS
while [ $SENS -le $TARGET_SENS ]; do
    # call prefilter module
    notExists "$TMP_PATH/pref_$SENS" && $RUNNER $MMSEQS prefilter "$INPUT" "$TARGET_DB_PREF" "$TMP_PATH/pref_$SENS" $PREFILTER_PAR -s $SENS && checkReturnCode "Prefilter died"
    # call alignment module
    notExists "$TMP_PATH/aln_$SENS"  && $RUNNER $MMSEQS align "$INPUT" "$TARGET" "$TMP_PATH/pref_$SENS" "$TMP_PATH/aln_$SENS" $ALIGNMENT_PAR  && checkReturnCode "Alignment died"

    if [ $SENS -gt $START_SENS ]; then
        $MMSEQS mergedbs "$1" "$TMP_PATH/aln_new" "$TMP_PATH/aln_${START_SENS}" "$TMP_PATH/aln_$SENS" \
            && checkReturnCode "Alignment died"
        mv -f "$TMP_PATH/aln_new" "$TMP_PATH/aln_${START_SENS}"
        mv -f "$TMP_PATH/aln_new.index" "$TMP_PATH/aln_${START_SENS}.index"
    fi

    NEXTINPUT="$TMP_PATH/input_step$SENS"
    if [  $SENS -lt $TARGET_SENS ]; then
        notExists "$TMP_PATH/order_step$SENS" \
            && awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$SENS.index" > "$TMP_PATH/order_step$SENS" \
            && checkReturnCode "Awk step $SENS died"
        notExists "$NEXTINPUT" \
            && $MMSEQS createsubdb "$TMP_PATH/order_step$SENS" "$INPUT" "$NEXTINPUT" \
            && checkReturnCode "Order step $SENS died"
    fi
    let SENS=SENS+SENS_STEP_SIZE

    INPUT=$NEXTINPUT
done
# post processing
mv -f "$TMP_PATH/aln_${START_SENS}" "$3"
mv -f "$TMP_PATH/aln_${START_SENS}.index" "$3.index"
checkReturnCode "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    SENS=$START_SENS
    while [ $SENS -lt $TARGET_SENS ]; do
        rm -f "$TMP_PATH/pref_$SENS" "$TMP_PATH/pref_$SENS.index"
        rm -f "$TMP_PATH/aln_$SENS" "$TMP_PATH/aln_$SENS.index"
        let SENS=SENS+SENS_STEP_SIZE
        NEXTINPUT="$TMP_PATH/input_step$SENS"
        rm -f "$TMP_PATH/input_step$SENS" "$TMP_PATH/input_step$SENS.index"
    done
fi


