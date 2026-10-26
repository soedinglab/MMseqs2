#!/bin/sh -e
# Fast profile-guided representative selection.
# Reuses the representative-to-member alignments produced by clustering (run with
# --include-align-files) instead of realigning a profile against every member, then
# rewrites the cluster DB with the best-scoring observed member as the new rep.
# Output is 1:1 with pickconsensusrep (a new cluster DB), only the core is faster.

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 4 ] && echo "Please provide <seqDB> <clusterDB> <outClusterDB> <tmpDir>" && exit 1;
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";
TMP_PATH="$4"

# Prerequisite: the clustering alignment DB (${clusterDB}_aln), produced by linclust/cluster
# when run with '--include-align-files 1 -a'. Without it there is nothing to reuse.
if notExists "${2}_aln.dbtype"; then
    fail "${2}_aln not found. Re-run linclust/cluster with '--include-align-files 1 -a' so the representative-to-member alignments are available for reuse."
fi

# 1. Score every observed member against its cluster profile and pick the best one.
if notExists "${TMP_PATH}/rep_map.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" pickrepprofile "$1" "${2}_aln" "${TMP_PATH}/rep_map" ${PICKREP_PAR} \
        || fail "pickrepprofile failed"
fi

# 2. Extract the oldRep -> newRep mapping (key space).
if notExists "${TMP_PATH}/rep_map.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "${TMP_PATH}/rep_map" "${TMP_PATH}/rep_map.tsv" --tsv ${VERBOSITY} \
        || fail "prefixid rep_map failed"
fi
awk '{ print $1 "\t" $2 }' "${TMP_PATH}/rep_map.tsv" > "${TMP_PATH}/rep_mapping.txt"

# 3. Old clustering as a flat oldRep -> member table.
if notExists "${TMP_PATH}/clu.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "$2" "${TMP_PATH}/clu.tsv" --tsv ${VERBOSITY} \
        || fail "prefixid clu failed"
fi

# 4. Rewrite the clustering with the new representatives.
awk 'FNR == NR{ f[$1] = $2; next }
     $1 != prev { print f[$1] "\t" f[$1]; prev = $1; }
     $1 in f && $2 != f[$1]{ print f[$1] "\t" $2 }' \
     "${TMP_PATH}/rep_mapping.txt" "${TMP_PATH}/clu.tsv" > "${TMP_PATH}/updated_clu.tsv"

# shellcheck disable=SC2086
"$MMSEQS" tsv2db "${TMP_PATH}/updated_clu.tsv" "$3" --output-dbtype 6 ${VERBOSITY} \
    || fail "tsv2db failed"

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/rep_map" ${VERBOSITY}
    rm -f "${TMP_PATH}/rep_map.tsv" "${TMP_PATH}/rep_mapping.txt" "${TMP_PATH}/clu.tsv" "${TMP_PATH}/updated_clu.tsv"
    rm -f "${TMP_PATH}/pickconsensusrepfast.sh"
fi
