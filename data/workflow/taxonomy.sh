#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
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

INPUT="$1"
TARGET="$2"
RESULTS="$3"
TMP_PATH="$4"

if [ ! -e "${TMP_PATH}/first.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${INPUT}" "${TARGET}" "${TMP_PATH}/first" "${TMP_PATH}/tmp_hsp1" ${SEARCH_PAR} \
        || fail "First search died"
fi
LCAIN="${TMP_PATH}/first"

if [ -n "${TOPHIT_MODE}" ]; then
  if [ ! -e "${TMP_PATH}/top1.dbtype" ]; then
      # shellcheck disable=SC2086
      "$MMSEQS" filterdb "${TMP_PATH}/first" "${TMP_PATH}/top1" --beats-first --filter-column 4 --comparison-operator le ${THREADS_COMP_PAR} \
          || fail "First filterdb died"
  fi
  LCAIN="${TMP_PATH}/top1"
fi

if [ "${TAX_OUTPUT}" -eq "0" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TARGET}" "${LCAIN}" "${RESULTS}" ${LCA_PAR} \
        || fail "Lca died"
elif [ "${TAX_OUTPUT}" -eq "2" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" lca "${TARGET}" "${LCAIN}" "${RESULTS}" ${LCA_PAR} \
        || fail "Lca died"
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${LCAIN}" "${RESULTS}_aln" ${VERBOSITY} \
        || fail "mvdb died"
else
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${LCAIN}" "${RESULTS}" ${VERBOSITY} \
        || fail "mvdb died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/first" ${VERBOSITY}
    if [ -n "${TOPHIT_MODE}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/top1" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/tmp_hsp1"
    rm -f "${TMP_PATH}/taxonomy.sh"
fi
