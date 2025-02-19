#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 4 ] && echo "Please provide <seqDB> <clusterDB> <outClusterDB> <tmp>" && exit 1;
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";
TMP_PATH="$4"

if notExists "${TMP_PATH}/msa.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" result2msa "$1" "$1" "$2" "$TMP_PATH/msa" ${RESULT2MSA_PAR} \
        || fail "result2msa failed"
fi

if notExists "${TMP_PATH}/profile.dbtype"; then
  # shellcheck disable=SC2086
  "$MMSEQS" msa2profile "$TMP_PATH/msa" "$TMP_PATH/profile" ${MSA2PROFILE_PAR} \
        || fail "result2msa failed"
fi

if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" align "$TMP_PATH/profile" "$1" "$2"  "${TMP_PATH}/aln" || fail "align failed"
fi

if notExists "${TMP_PATH}/aln.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "${TMP_PATH}/aln" "${TMP_PATH}/aln.tsv" --tsv ${VERBOSITY} || fail "prefixid1 -tsv failed"
fi

awk 'FNR == NR{ best[$1]=1; rep[$1] = $1; next}
{
    cluster = $1; member = $2; score = $3;
    if (!(cluster in best) || score > best[cluster]) {
        best[cluster] = score;
        rep[cluster] = member;
    }
}
END {
    for (cluster in rep) {
        print cluster "\t" rep[cluster];
    }
}' "${TMP_PATH}/aln.index" "${TMP_PATH}/aln.tsv" > "${TMP_PATH}/rep_mapping.txt"

if notExists "${TMP_PATH}/clu.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "$2" "${TMP_PATH}/clu.tsv" --tsv ${VERBOSITY} || fail "prefixid2 -tsv failed"
fi

if notExists "${TMP_PATH}/updated_clu.tsv"; then
    awk 'FNR == NR{f[$1] = $2; next}
         $1 != prev { print f[$1] "\t" f[$1]; prev = $1; }
         $1 in f && $2 != f[$1]{print f[$1]"\t"$2}' "${TMP_PATH}/rep_mapping.txt" "${TMP_PATH}/clu.tsv" > "${TMP_PATH}/updated_clu.tsv"
fi

"$MMSEQS" tsv2db "${TMP_PATH}/updated_clu.tsv" "${3}" --output-dbtype 6 || fail "tsv2db failed"

if [ -n "$REMOVE_TMP" ]; then
    rm -f "${TMP_PATH}/updated_clu.tsv"
    rm -f "${TMP_PATH}/aln.tsv"
    rm -f "${TMP_PATH}/rep_mapping.txt"
    rm -f "${TMP_PATH}/clu.tsv"
    rm -f "${TMP_PATH}/updated_clu.tsv"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/msa" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/profile" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}

    rm -rf "${TMP_PATH}/pickconsensusrep.sh"
fi
