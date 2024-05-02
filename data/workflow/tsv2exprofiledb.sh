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

if notExists "${OUT}_h.dbtype"; then
  "$MMSEQS" tsv2db "${IN}_h.tsv" "${OUT}_h" --output-dbtype 12 ${VERBOSITY}
fi

if notExists "${OUT}.dbtype"; then
  "$MMSEQS" tsv2db "${IN}.tsv" "${OUT}_tmp" --output-dbtype 0 ${VERBOSITY}
  MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_tmp" "${OUT}" ${VERBOSITY}
  "$MMSEQS" rmdb "${OUT}_tmp" ${VERBOSITY}
fi

if notExists "${OUT}_seq.dbtype"; then
  "$MMSEQS" tsv2db "${IN}_seq.tsv" "${OUT}_seq_tmp" --output-dbtype 0 ${VERBOSITY}
  MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_seq_tmp" "${OUT}_seq" ${VERBOSITY}
  "$MMSEQS" rmdb "${OUT}_seq_tmp" ${VERBOSITY}
fi

if notExists "${OUT}_aln.dbtype"; then
  "$MMSEQS" tsv2db "${IN}_aln.tsv" "${OUT}_aln_tmp" --output-dbtype 5 ${VERBOSITY}
  MMSEQS_FORCE_MERGE=1 "$MMSEQS" compress "${OUT}_aln_tmp" "${OUT}_aln" ${VERBOSITY}
  "$MMSEQS" rmdb "${OUT}_aln_tmp" ${VERBOSITY}
fi

if notExists "${OUT}_seq_h.dbtype"; then
  "$MMSEQS" aliasdb "${OUT}_h" "${OUT}_seq_h" ${VERBOSITY}
fi

if [ -e "${OUT}.sh" ]; then
  rm -f -- "${OUT}.sh"
fi
