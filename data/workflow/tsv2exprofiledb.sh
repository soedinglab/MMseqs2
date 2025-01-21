#!/bin/sh -e
# shellcheck disable=SC2086
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 2 ] && echo "Please provide <inputTSV> <outDB>" && exit 1

notExists() {
	[ ! -f "$1" ]
}

IN="$1"
OUT="$2"

[ ! -f "${IN}.tsv" ] && echo "${IN}.tsv not found!" && exit 1;
[ ! -f "${IN}_h.tsv" ] && echo "${IN}_h.tsv not found!" && exit 1;
[ ! -f "${IN}_seq.tsv" ] && echo "${IN}_seq.tsv not found!" && exit 1;
[ ! -f "${IN}_aln.tsv" ] && echo "${IN}_aln.tsv not found!" && exit 1;
[ -d "${OUT}.tsv" ] && echo "${OUT} is a directory!" && exit 1;

if notExists "${OUT}_seq.dbtype"; then
  if [ -n "${COMPRESSED}" ]; then
    "$MMSEQS" tsv2db "${IN}_seq.tsv" "${OUT}_seq_tmp" --output-dbtype 0 ${VERBOSITY}
    MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_seq_tmp" "${OUT}_seq" ${THREADS}
    "$MMSEQS" rmdb "${OUT}_seq_tmp" ${VERBOSITY}
  else
    "$MMSEQS" tsv2db "${IN}_seq.tsv" "${OUT}_seq" --output-dbtype 0 ${VERBOSITY}
  fi
fi

if notExists "${OUT}_seq_h.dbtype"; then
  MMSEQS_FORCE_MERGE=1 "$MMSEQS" tsv2db "${IN}_h.tsv" "${OUT}_seq_h" --output-dbtype 12 ${VERBOSITY}
fi

if notExists "${OUT}.dbtype"; then
  if [ -n "${COMPRESSED}" ]; then
    "$MMSEQS" tsv2db "${IN}.tsv" "${OUT}_tmp" --output-dbtype 0 ${VERBOSITY}
    MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_tmp" "${OUT}" ${THREADS}
    "$MMSEQS" rmdb "${OUT}_tmp" ${VERBOSITY}
  else
    MMSEQS_FORCE_MERGE=1 "$MMSEQS" tsv2db "${IN}.tsv" "${OUT}" --output-dbtype 0 ${VERBOSITY}
  fi
fi

if [ -n "${GPU}" ]; then
  if notExists "${OUT}.GPU_READY"; then
    "$MMSEQS" aliasdb "${OUT}_seq_h" "${OUT}_h" ${VERBOSITY}
    "$MMSEQS" makepaddedseqdb "${OUT}" "${OUT}_pad" ${THREADS}
    "$MMSEQS" rmdb "${OUT}" ${VERBOSITY}
    "$MMSEQS" rmdb "${OUT}_h" ${VERBOSITY}
    "$MMSEQS" mvdb "${OUT}_pad" "${OUT}" ${VERBOSITY}
    "$MMSEQS" mvdb "${OUT}_pad_h" "${OUT}_h" ${VERBOSITY}
    touch "${OUT}.GPU_READY"
  fi
else
  if notExists "${OUT}_h.dbtype"; then
    "$MMSEQS" aliasdb "${OUT}_seq_h" "${OUT}_h" ${VERBOSITY}
  fi
fi

if notExists "${OUT}_aln.dbtype"; then
  TSVPATH="${IN}"
  if [ -n "${GPU}" ]; then
    if notExists "${OUT}_mapped_aln.tsv"; then
      awk 'BEGIN { FS = "\t"; OFS = "\t"; sh = ""; entry = ""; } NR == 1 { last = $1; } ($1 != last) { if (sh != "") print sh; if (entry != "") print entry; last = $1; sh = ""; entry = ""; } ($1 != "" && $1 == $2) { sh = $0; next; } { if (entry != "") { entry = entry "\n" $0; } else { entry = $0; } } END { if (sh != "") print sh; if (entry != "") print entry; }' \
        "${IN}_aln.tsv" > "${OUT}_reorder_aln.tsv"
      awk 'NR == FNR { f[$3] = $1; next; } { $1 = f[$1]; print }' \
        "${OUT}.lookup" "${OUT}_reorder_aln.tsv" | sort -s -k1,1n > "${OUT}_mapped_aln.tsv"
      rm -f -- "${OUT}_reorder_aln.tsv"
    fi
    TSVPATH="${OUT}_mapped"
  fi

  if [ -n "${COMPRESSED}" ]; then
    "$MMSEQS" tsv2db "${TSVPATH}_aln.tsv" "${OUT}_aln_tmp" --output-dbtype 5 ${VERBOSITY}
    MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_aln_tmp" "${OUT}_aln" ${VERBOSITY}
    "$MMSEQS" rmdb "${OUT}_aln_tmp" ${VERBOSITY}
  else
    MMSEQS_FORCE_MERGE=1 "$MMSEQS" tsv2db "${TSVPATH}_aln.tsv" "${OUT}_aln" --output-dbtype 5 ${VERBOSITY}
  fi
fi

if [ -e "${OUT}.sh" ]; then
  if [ -n "${GPU}" ]; then
    rm -rf -- "${OUT}_mapped_aln.tsv"
  fi
  rm -f -- "${OUT}.sh"
fi
