#!/bin/sh -e
# Iterative sequence search workflow script
fail() {
  echo "Error: $1"
  exit 1
}

notExists() {
  [ ! -f "$1" ]
}

# pre-processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -lt 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;

# check if input files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# TODO: ask for $2_aln and that this contains backtrace
#[ ! -f "$2_aln" ] && echo "$2_aln not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERYDB="$1"
TARGETDB="$2"
ALNDB="${2}_aln"
SEQDB="${2}_seq"
TMP_PATH="$4"
STEP=0

while [ "$STEP" -lt "$NUM_IT" ]; do
  # call slice search for the first iteration
  if [ "$STEP" -eq 0 ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "$QUERYDB" "$TARGETDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH" ${SEARCH_PAR} \
      || fail "Slicesearch died"
    # shellcheck disable=SC2086
    "$MMSEQS" profile2consensus "$TARGETDB" "$2_consensus" ${CONSENSUS_PAR} \
      || fail "Profile2Consensus died"
    TARGETDB="$2_consensus"
  fi
  # call prefilter module
  if [ "$STEP" -gt 0 ]; then
    PARAM="PREFILTER_PAR_$STEP"
    eval TMP="\$$PARAM"
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "$QUERYDB" "$TARGETDB" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
      || fail "Prefilter died"
    STEPPREV=$((STEP-1))
    # shellcheck disable=SC2086
    "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPPREV" "$TMP_PATH/pref_$STEP" $SUBTRACT_PAR \
      || fail "Subtract died"
    "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
    # call alignment module
    PARAM="ALIGNMENT_PAR_$STEP"
    eval TMP="\$$PARAM"
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$QUERYDB" "$2" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_tmp_$STEP" ${TMP} \
      || "Alignment died"
    # merge alignment dbs
    STEPPREV=$((STEP-1))
    "$MMSEQS" mergedbs "$QUERYDB" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPPREV" "$TMP_PATH/aln_tmp_$STEP" \
      || fail "Mergedbs died"
    #"$MMSEQS" rmdb "$TMP_PATH/aln_$STEPPREV"
    #"$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
  fi
  # expand alignment dbs
  if [ "$STEP" -ne $((NUM_IT - 1)) ]; then
    # shellcheck disable=SC2086
      "$MMSEQS" expand2profile "$QUERYDB" "$SEQDB" "$TMP_PATH/aln_$STEP" "$ALNDB" "$TMP_PATH/profile_$STEP" $EXPANDPROFILE_PAR \
      || fail 'Expand2Profile died'
  else
#    PARAM="EXPANDALN_PAR"
    # shellcheck disable=SC2086
    "$MMSEQS" expandaln "$QUERYDB" "$SEQDB" "$TMP_PATH/aln_$STEP" "$ALNDB" "$3" $EXPANDALN_PAR \
      || fail "Expandaln died"
  fi
  QUERYDB="$TMP_PATH/profile_$STEP"
  STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
  STEP=0
  while [ "$STEP" -lt "$NUM_IT" ]; do
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
    STEP=$((STEP+1))
  done
  rm -f "$TMP_PATH/iterativepp.sh"
fi
